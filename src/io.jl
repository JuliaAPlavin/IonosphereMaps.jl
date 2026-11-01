openmaybegz(f, path) = endswith(lowercase(path), ".gz") ? open(io -> f(GzipDecompressorStream(io)), path) : open(f, path)
firstline(path) = openmaybegz(readline, path)
headerlabel(l) = strip(l[min(end + 1, 61):end])

isionex(path) = isfile(path) && occursin("IONEX VERSION / TYPE", firstline(path))
isrinexobs(path) = isfile(path) && (l = firstline(path); occursin("RINEX VERSION / TYPE", l) && l[21] == 'O')

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

# GPS differential code biases (ns) from the IONEX auxiliary block
function read_ionex_dcb(path)
    sat = Dict{Int, Float64}(); sta = Dict{String, Float64}(); ref = nothing
    openmaybegz(path) do io
        for l in eachline(io)
            lab = headerlabel(l)
            lab == "END OF HEADER" && break
            if lab == "PRN / BIAS / RMS" && l[4] == 'G'
                sat[parse(Int, l[5:6])] = parse(Float64, split(l[1:60])[end-1])
            elseif lab == "STATION / BIAS / RMS" && l[4] == 'G'
                sta[uppercase(l[7:10])] = parse(Float64, split(l[1:60])[end-1])
            elseif occursin(r"Reference observables for GPS\s*:", l)
                ref = strip(split(l[1:60], ':')[2])
            end
        end
    end
    isempty(sat) && error("no GPS DCB block in $path")
    (; sat, sta, ref)
end

# satellite observable-specific code biases (ns) from a Bias-SINEX OSB file, keyed by (prn, code)
function read_osb(path)
    openmaybegz(path) do io
        filtermap(eachline(io)) do l
            v = split(l)
            startswith(l, " OSB ") && length(v) >= 9 && occursin(r"^G\d\d$", v[3]) && !occursin(r"^\d", v[4]) || return nothing
            (parse(Int, v[3][2:3]), String(v[4])) => parse(Float64, v[end-1])
        end |> Dict
    end
end

const OBS_CODES = (:C1C, :C1W, :C2W, :L1C, :L2W)

# GPS records of a RINEX 3 observation file; times in GPS time; NaN for absent observables
function read_rinex3_gps(path)
    openmaybegz(path) do io
        types = String[]; ng = 0; xyz = nothing; dh = 0.0; marker = ""; leap = nothing
        while true
            l = readline(io)
            lab = headerlabel(l)
            lab == "END OF HEADER" && break
            lab == "MARKER NAME" && (marker = strip(l[1:60]))
            lab == "APPROX POSITION XYZ" && (xyz = SVector{3}(parsefloats(l[1:60])))
            lab == "ANTENNA: DELTA H/E/N" && (dh = parsefloats(l[1:60])[1])
            lab == "LEAP SECONDS" && (leap = parse(Int, split(l[1:60])[1]))
            if lab == "SYS / # / OBS TYPES" && (l[1] == 'G' || (l[1] == ' ' && length(types) < ng))
                l[1] == 'G' && (ng = parse(Int, l[4:6]))
                append!(types, split(l[8:60]))
            end
        end
        idx = map(c -> something(findfirst(==(String(c)), types), 0), OBS_CODES)
        @assert idx[3] > 0 && (idx[1] > 0 || idx[2] > 0) && idx[4] > 0 && idx[5] > 0 "$path lacks GPS C1C/C1W, C2W, L1C, L2W"
        field(l, k) = (r = 4+16(k-1):3+16(k-1)+14; k == 0 || length(l) < last(r) || isempty(strip(l[r])) ? NaN : parse(Float64, l[r]))
        lli(l, k) = (j = 3 + 16(k-1) + 15; length(l) >= j && l[j] != ' ' && isodd(parse(Int, l[j])))
        lines = readlines(io)
        obs = @p findall(startswith('>'), lines) |> filter(parse(Int, lines[_][32]) <= 1) |> flatmap() do i
            h = lines[i]
            t = DateTime(parse.(Int, (h[3:6], h[8:9], h[11:12], h[14:15], h[17:18]))...) +
                Millisecond(round(Int, 1000 * parse(Float64, h[19:29])))
            @p lines[i+1:i+parse(Int, h[33:35])] |> filter(startswith('G')) |> map() do l
                (; time=t, prn=parse(Int, l[2:3]), NamedTuple{OBS_CODES}(map(k -> field(l, k), idx))...,
                   lli=lli(l, idx[4]) || lli(l, idx[5]))
            end
        end |> StructArray
        (; marker, xyz=xyz + dh * normalize(xyz), leap, obs)
    end
end

# SP3 precise orbits: ECEF positions (m) of GPS satellites; times in GPS time
function read_sp3(path)
    lines = openmaybegz(readlines, path)
    any(l -> startswith(l, "%c") && occursin(" GPS", l), lines) || error("$path: SP3 time system is not GPS")
    eps = findall(l -> startswith(l, "* "), lines)
    times = map(eps) do i
        v = split(lines[i])
        DateTime(parse.(Int, v[2:6])...) + Millisecond(round(Int, 1000 * parse(Float64, v[7])))
    end
    nanpos = SVector(NaN, NaN, NaN)
    prns = sort(unique(filtermap(l -> startswith(l, "PG") ? parse(Int, l[3:4]) : nothing, lines)))
    pos = Dict(map(prns) do prn
        prn => map(eachindex(eps)) do e
            stop = e < length(eps) ? eps[e+1] : length(lines) + 1
            i = findfirst(l -> startswith(l, "PG") && parse(Int, l[3:4]) == prn, lines[eps[e]+1:stop-1])
            isnothing(i) ? nanpos : 1e3 * SVector{3}(parsefloats(lines[eps[e]+i][5:46]))
        end
    end)
    (; times, pos)
end

# 10-point Lagrange interpolation of an SP3 orbit at time t (GPS time)
function satpos(sp3, prn, t; n=10)
    ts = sp3.times; p = sp3.pos[prn]
    dt = Dates.value(ts[2] - ts[1]); x = Dates.value(t - ts[1]) / dt
    i0 = clamp(round(Int, x) - n ÷ 2 + 1, 1, length(ts) - n + 1)
    nodes = i0:i0+n-1
    sum(nodes) do j
        prod(m -> m == j ? 1.0 : (x - (m - 1)) / (j - m), nodes) * p[j]
    end
end
