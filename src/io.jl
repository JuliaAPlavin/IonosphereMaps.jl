# gzip also reads unix compress (.Z)
readlines_decompressed(path) = occursin(r"\.(gz|z)$"i, path) ? readlines(`gzip -dc $path`) : readlines(path)
headerlabel(l) = strip(l[min(end + 1, 61):end])

# fields of width `w`: IONEX numbers need not be separated by spaces
fixedwidth(f, s, w) = let s = rstrip(s)
    map(i -> f(s[i:min(i + w - 1, end)]), 1:w:length(s))
end
fixedwidth(T::Type, s, w) = fixedwidth(Base.Fix1(parse, T), s, w)

# a map value scaled by 10^expo; no value (9999, or `nan` from some producers) as NaN
mapvalue(s, expo) = lowercase(strip(s)) == "nan" ? NaN : let v = parse(Int, s)
    v == 9999 ? NaN : v * 10.0^expo
end

# map block starting at line `i` as M[lat, lon], scaled by the EXPONENT in force
function readmap(lines, i, expos, lats, lons)
    endlabel = replace(headerlabel(lines[i]), "START" => "END")
    stop = findnext(l -> headerlabel(l) == endlabel, lines, i)
    rows = filter(j -> headerlabel(lines[j]) == "LAT/LON1/LON2/DLON/H", i:stop)
    @assert length(rows) == length(lats)
    stack(zip(rows, lats); dims=1) do (r, lat)
        rlat, lon1, lon2, dlon = fixedwidth(Float64, lines[r][3:26], 6)
        @assert rlat ≈ lat && lon1:dlon:lon2 == lons
        vals = collect(Iterators.flatmap(l -> fixedwidth(v -> mapvalue(v, expos[r]), l, 5), lines[r+1:r+cld(length(lons), 16)]))
        @assert length(vals) == length(lons)
        vals
    end
end

# TEC maps (TECU) as tec[lat, lon, epoch] with degree and UTC keys, and their RMS maps (or `nothing`) as rms
function read_ionex(path)
    lines = readlines_decompressed(path)
    hdr(lbl) = lines[findfirst(l -> headerlabel(l) == lbl, lines)][1:60]
    # an EXPONENT record (default -1) holds until the next one, which may come inside the data section
    expos = accumulate((e, l) -> headerlabel(l) == "EXPONENT" ? parse(Int, l[1:6]) : e, lines; init=-1)
    lat1, lat2, dlat = fixedwidth(Float64, hdr("LAT1 / LAT2 / DLAT")[3:20], 6)
    lon1, lon2, dlon = fixedwidth(Float64, hdr("LON1 / LON2 / DLON")[3:20], 6)
    lats = lat1:dlat:lat2; lons = lon1:dlon:lon2
    starts = findall(l -> headerlabel(l) == "START OF TEC MAP", lines)
    # the map count is what reveals a file cut between maps: .Z has no checksum, END OF FILE is often absent
    nmaps = parse(Int, hdr("# OF MAPS IN FILE")[1:6])
    length(starts) == nmaps || throw(ArgumentError("$path has $(length(starts)) TEC maps, its header says $nmaps"))
    epochs = map(starts) do i
        e = lines[findnext(l -> headerlabel(l) == "EPOCH OF CURRENT MAP", lines, i)]
        DateTime(fixedwidth(Int, e[1:36], 6)...)   # hour 24 rolls over to next day
    end
    rmsstarts = findall(l -> headerlabel(l) == "START OF RMS MAP", lines)
    # RMS maps pair with TEC maps by index
    mapindex(i) = parse(Int, lines[i][1:6])
    @assert isempty(rmsstarts) || map(mapindex, rmsstarts) == map(mapindex, starts)
    keyed(M) = KeyedArray(M; lat=lats, lon=lons, epoch=epochs)
    (; tec=keyed(stack(map(i -> readmap(lines, i, expos, lats, lons), starts))),
       rms=isempty(rmsstarts) ? nothing : keyed(stack(map(i -> readmap(lines, i, expos, lats, lons), rmsstarts))),
       height=fixedwidth(Float64, hdr("HGT1 / HGT2 / DHGT")[3:20], 6)[1] * u"km", radius=parse(Float64, hdr("BASE RADIUS")[1:8]) * u"km",
       mapping=strip(hdr("MAPPING FUNCTION")))
end
