# gzip also reads unix compress (.Z)
readlines_decompressed(path) = occursin(r"\.(gz|z)$"i, path) ? readlines(`gzip -dc $path`) : readlines(path)
# label from column 61 on, counted in characters: a non-IONEX file need not be ASCII
headerlabel(l) = strip(chop(l; head=60, tail=0))

# index of the first `label` record from line `from` on
findrecord(path, lines, label, from=1) =
    @something(findnext(l -> headerlabel(l) == label, lines, from), throw(ArgumentError("$path has no $label record from line $from on")))

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
function readmap(path, lines, i, expos, lats, lons)
    stop = findrecord(path, lines, replace(headerlabel(lines[i]), "START" => "END"), i)
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
    hdr(lbl) = lines[findrecord(path, lines, lbl)][1:60]
    dim = parse(Int, hdr("MAP DIMENSION")[1:6])
    dim == 2 || throw(ArgumentError("$path has $dim-D maps, only 2-D maps are supported"))
    # an EXPONENT record (default -1) holds until the next one, which may come inside the data section
    expos = accumulate((e, l) -> headerlabel(l) == "EXPONENT" ? parse(Int, l[1:6]) : e, lines; init=-1)
    lat1, lat2, dlat = fixedwidth(Float64, hdr("LAT1 / LAT2 / DLAT")[3:20], 6)
    lon1, lon2, dlon = fixedwidth(Float64, hdr("LON1 / LON2 / DLON")[3:20], 6)
    lats = lat1:dlat:lat2; lons = lon1:dlon:lon2
    starts = findall(l -> headerlabel(l) == "START OF TEC MAP", lines)
    # the map count is what reveals a file cut between maps: .Z has no checksum, END OF FILE is often absent
    nmaps = parse(Int, hdr("# OF MAPS IN FILE")[1:6])
    length(starts) == nmaps || throw(ArgumentError("$path has $(length(starts)) TEC maps, its header says $nmaps"))
    # hour 24 rolls over to the next day
    epoch(i) = DateTime(fixedwidth(Int, lines[findrecord(path, lines, "EPOCH OF CURRENT MAP", i)][1:36], 6)...)
    epochs = map(epoch, starts)
    # interpolation and joining files look maps up by epoch
    issorted(epochs; lt=<=) || throw(ArgumentError("$path has TEC map epochs not strictly increasing: $epochs"))
    rmsstarts = findall(l -> headerlabel(l) == "START OF RMS MAP", lines)
    # RMS maps pair with TEC maps by position
    rmsepochs = map(epoch, rmsstarts)
    isempty(rmsstarts) || rmsepochs == epochs || throw(ArgumentError("$path has RMS maps at $rmsepochs, its TEC maps at $epochs"))
    keyed(M) = KeyedArray(M; lat=lats * u"°", lon=lons * u"°", epoch=epochs)
    readmaps(blocks) = keyed(stack(map(i -> readmap(path, lines, i, expos, lats, lons), blocks)))
    (; tec=readmaps(starts), rms=isempty(rmsstarts) ? nothing : readmaps(rmsstarts),
       height=fixedwidth(Float64, hdr("HGT1 / HGT2 / DHGT")[3:20], 6)[1] * u"km", radius=parse(Float64, hdr("BASE RADIUS")[1:8]) * u"km",
       mapping=strip(hdr("MAPPING FUNCTION")))
end
