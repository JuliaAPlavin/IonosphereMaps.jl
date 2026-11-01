const C_LIGHT = 299792458.0
const F1, F2 = 1575.42e6, 1227.60e6
const TECU_PER_METRE = 1 / (40.3082e16 * (1 / F2^2 - 1 / F1^2))   # slant TEC per metre of P2 − P1

"""
    MeanVTEC(; height=300e3, elmin=deg2rad(10))

ALBUS `RI_G01`: VTEC = mean over all satellites above `elmin` at the epoch nearest `t` of STEC·cos z′
(thin shell at `height` over a 6378.137 km sphere), mapped to the target line of sight with the same shell;
B∥ from IGRF at the target's shell pierce point. Equal weights (ALBUS weights by 1/σ of the levelling).
"""
Base.@kwdef struct MeanVTEC
    height::Float64 = 300e3
    elmin::Float64 = deg2rad(10)
end

"""
    PlanarVTEC(; height=350e3, elmin=deg2rad(15), window=Minute(10))

VTEC(pierce point) fitted as a plane in latitude/longitude plus a linear time term, to all observations above
`elmin` within ±`window` of `t` (weights sin²el); evaluated at the target pierce point (shell over 6371 km).
"""
Base.@kwdef struct PlanarVTEC
    height::Float64 = 350e3
    elmin::Float64 = deg2rad(15)
    window::Millisecond = Minute(10)
end

"""
    SingleStation(rinex_obs; sp3, ionex, osb, model=MeanVTEC(), receiver_dcb)

Slant TEC from one dual-frequency GPS receiver: geometry-free code (C2W − C1W, or C1C with OSB-derived
satellite C1C−C1W biases) levelled to the carrier phase per continuous arc, corrected for satellite and receiver
DCBs from the IONEX auxiliary block (CODE: C1W−C2W). By default, `sp3`, `ionex` and `osb` are the unique matching
files in the directory of `rinex_obs` (`osb` only needed for receivers without C1W).
"""
struct SingleStation{M, T}
    xyz::SVector{3, Float64}
    obs::T
    model::M
end

sibling(path, re) = (fs = filter(f -> occursin(re, basename(f)), readdir(dirname(abspath(path)); join=true));
                     length(fs) == 1 ? only(fs) : error("expected exactly one file matching $re next to $path, found $(basename.(fs))"))
optsibling(path, re) = (fs = filter(f -> occursin(re, basename(f)), readdir(dirname(abspath(path)); join=true)); isempty(fs) ? nothing : only(fs))

function SingleStation(path::AbstractString;
        sp3=sibling(path, r"\.sp3(\.gz)?$"i), ionex=sibling(path, r"(\.inx|\.\d\di)(\.gz)?$"i),
        osb=optsibling(path, r"osb\.bia(\.gz)?$"i), model=MeanVTEC(), receiver_dcb=nothing)
    R = read_rinex3_gps(path)
    orbits = read_sp3(sp3)
    dcb = read_ionex_dcb(ionex)
    something(dcb.ref, "C1W-C2W") == "C1W-C2W" || error("$ionex: GPS DCB reference $(dcb.ref), expected C1W-C2W")
    isnothing(dcb.ref) && @warn "$ionex does not state the DCB reference observables; assuming C1W-C2W (P1-P2)"
    rdcb = @something receiver_dcb get(dcb.sta, uppercase(R.marker[1:4]), nothing) error("no receiver DCB for $(R.marker) in $ionex; pass receiver_dcb")
    leap = @something R.leap (@warn("$path has no LEAP SECONDS header; using 18 s"); 18)
    usec1c = all(isnan, R.obs.C1W)
    p1c1 = usec1c ? (isnothing(osb) ? error("$path has no C1W: OSB file needed for satellite C1C−C1W biases") : read_osb(osb)) : nothing
    satbias(prn) = dcb.sat[prn] + (usec1c ? p1c1[(prn, "C1C")] - p1c1[(prn, "C1W")] : 0.0)
    T = geometry_free(R, orbits, usec1c, leap)
    L = level(T)
    nodcb = setdiff(unique(map(r -> r.prn, L)), keys(dcb.sat))
    isempty(nodcb) || @warn "no satellite DCB in $ionex for PRNs $nodcb; their data are not used"
    obs = @p L |> filterview(haskey(dcb.sat, _.prn)) |>
        map((; _.time, _.prn, _.el, _.k, stec=(_.P4 + C_LIGHT * 1e-9 * (satbias(_.prn) + rdcb)) * TECU_PER_METRE)) |>
        StructArray |> sort(by=_.time)
    SingleStation(R.xyz, obs, model)
end

function geometry_free(R, orbits, usec1c, leap)
    u = up(geodetic(R.xyz))
    @p R.obs |> filterview(!isnan(usec1c ? _.C1C : _.C1W) && !isnan(_.C2W) && !isnan(_.L1C) && !isnan(_.L2W) && haskey(orbits.pos, _.prn)) |>
        map() do o
            C1 = usec1c ? o.C1C : o.C1W
            k = normalize(satpos(orbits, o.prn, o.time) - R.xyz)
            L4 = C_LIGHT / F1 * o.L1C - C_LIGHT / F2 * o.L2W
            (; time=o.time - Second(leap), o.prn, k, el=asin(dot(k, u)), o.lli,
               P4=o.C2W - C1, L4, MW=(C_LIGHT * (o.L1C - o.L2W)) / (F1 - F2) - (F1 * C1 + F2 * o.C2W) / (F1 + F2))
        end |> StructArray
end

# continuous phase arcs (break on gap, loss of lock, geometry-free or Melbourne-Wübbena jump);
# phase levelled to code by the sin²el-weighted mean of P4 − L4 above `elmin`
function level(T; maxgap=Second(90), gfjump=0.15, mwjump=3.0, minlen=20, elmin=deg2rad(10))
    arcs = flatmap(unique(T.prn)) do prn
        ix = @p findall(==(prn), T.prn) |> sort(by=T.time[_])
        brk = map(zip(ix[1:end-1], ix[2:end])) do (a, b)
            T.time[b] - T.time[a] > maxgap || T.lli[b] || abs(T.L4[b] - T.L4[a]) > gfjump || abs(T.MW[b] - T.MW[a]) > mwjump
        end
        starts = [1; findall(brk) .+ 1]; stops = [starts[2:end] .- 1; length(ix)]
        map((s, e) -> ix[s:e], starts, stops)
    end
    @p arcs |> filter(length(_) >= minlen && any(i -> T.el[i] >= elmin, _)) |> flatmap() do a
        sel = filter(i -> T.el[i] >= elmin, a)
        w = sin.(T.el[sel]) .^ 2
        offset = sum(w .* (T.P4[sel] .- T.L4[sel])) / sum(w)
        map(i -> (; time=T.time[i], prn=T.prn[i], k=T.k[i], el=T.el[i], P4=T.L4[i] + offset), a)
    end
end

const ALBUS_RADIUS = 6378137.0
cosz(el, hobs, H) = sqrt(1 - ((ALBUS_RADIUS + hobs) * cos(el) / (ALBUS_RADIUS + H))^2)

function slant_tec(m::SingleStation{MeanVTEC}, o, site, k, t)
    (; height, elmin) = m.model
    (; obs) = m
    i = searchsortedfirst(obs.time, t)
    cand = filter(in(eachindex(obs.time)), (i - 1, i))
    isempty(cand) && return (; stec=NaN, pp=o)
    ep = obs.time[argmin(j -> abs(obs.time[j] - t), cand)]
    abs(ep - t) <= Minute(1) || return (; stec=NaN, pp=o)
    hrx = geodetic(m.xyz).h
    sel = filterview(s -> s.el >= elmin, view(obs, searchsorted(obs.time, ep)))
    isempty(sel) && return (; stec=NaN, pp=o)
    vtec = sum(s -> s.stec * cosz(s.el, hrx, height), sel) / length(sel)
    (; stec=vtec / cosz(asin(dot(k, up(site))), site.h, height), pp=pierce(o, k, ALBUS_RADIUS + height))
end

function slant_tec(m::SingleStation{PlanarVTEC}, o, site, k, t)
    (; height, elmin, window) = m.model
    r = 6371e3 + height
    pp = pierce(o, k, r)
    c = geocentric_latlon(pp)
    sel = filterview(s -> s.el >= elmin, view(m.obs, searchsortedfirst(m.obs.time, t - window):searchsortedlast(m.obs.time, t + window)))
    length(sel) >= 10 || return (; stec=NaN, pp)
    rows = @p sel |> map() do s
        q = pierce(m.xyz, s.k, r)
        g = geocentric_latlon(q)
        (; x=SVector(1.0, g.lat - c.lat, (g.lon - c.lon) * cos(c.lat), Dates.value(s.time - t) / 3.6e6),
           y=s.stec * dot(q, s.k) / norm(q), w=sin(s.el)^2)
    end |> StructArray
    X = stack(rows.x; dims=1); W = Diagonal(rows.w)
    β = (X' * W * X) \ (X' * W * rows.y)
    (; stec=β[1] * norm(pp) / dot(pp, k), pp)
end
