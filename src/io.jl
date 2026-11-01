openmaybegz(f, path) = endswith(lowercase(path), ".gz") ? open(io -> f(GzipDecompressorStream(io)), path) : open(f, path)
headerlabel(l) = strip(l[min(end + 1, 61):end])

parsefloats(s) = parse.(Float64, split(s))

# TEC maps (TECU) as tec[lat, lon, epoch]; epochs in UTC
function read_ionex(path)
    lines = openmaybegz(readlines, path)
    hdr(lbl) = lines[findfirst(l -> headerlabel(l) == lbl, lines)][1:60]
    expo = parse(Int, split(hdr("EXPONENT"))[1])
    la = parsefloats(hdr("LAT1 / LAT2 / DLAT")); lo = parsefloats(hdr("LON1 / LON2 / DLON"))
    lats = la[1]:la[3]:la[2]; lons = lo[1]:lo[3]:lo[2]
    starts = findall(l -> headerlabel(l) == "START OF TEC MAP", lines)
    maps = map(starts) do i
        epoch = DateTime(parse.(Int, split(lines[i+1])[1:6])...)   # hour 24 rolls over to next day
        rows = findall(l -> headerlabel(l) == "LAT/LON1/LON2/DLON/H", lines[i+2:findnext(l -> headerlabel(l) == "END OF TEC MAP", lines, i)])
        M = stack(rows; dims=1) do r
            vals = Int[]; j = i + 1 + r
            while length(vals) < length(lons)
                j += 1; append!(vals, parse.(Int, split(lines[j])))
            end
            @. ifelse(vals == 9999, NaN, vals * 10.0^expo)
        end
        @assert size(M) == (length(lats), length(lons))
        (; epoch, M)
    end
    (; epochs=map(m -> m.epoch, maps), lats, lons, tec=stack(map(m -> m.M, maps)),
       height=1e3 * parsefloats(hdr("HGT1 / HGT2 / DHGT"))[1], radius=1e3 * parse(Float64, split(hdr("BASE RADIUS"))[1]),
       mapping=strip(hdr("MAPPING FUNCTION")))
end

