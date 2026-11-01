"""Single-layer mapping at the IONEX shell height: exact 1/cos z′ at the pierce point."""
struct SLM end
"""CODE modified single-layer mapping (H = 506.7 km, α = 0.9782, R = 6371 km) at the geodetic zenith distance."""
struct MSLM end

"""
    GlobalMap(paths; mapping)
    GlobalMap(path; mapping)

Thin-shell model from IONEX global ionosphere maps: VTEC bilinearly interpolated at the pierce point,
maps blended in time after rotating each with the Sun (IONEX method 3), also across file boundaries.
`mapping` (`SLM()` or `MSLM()`) converts VTEC to slant TEC; use the one the maps were produced with.

`paths` are files of consecutive periods in time order, e.g. the daily files of an observation that
crosses midnight. They must share the grid, base radius, shell height and mapping function, and each
must start no later than one map interval after the previous one ends.
Where consecutive files overlap — normally the midnight map that closes one daily file and opens the
next, from two separate daily solutions — the later file's maps are used. IONEX assigns a time to the
interval `Tᵢ ≤ t < Tᵢ₊₁`, so a day file covers `[00:00, 24:00)` and its 24:00 map only closes the last
interval; the epoch itself belongs to the next day. AIPS TECOR (NFILES > 1) and spinifex
(`remove_midnight_jumps`) treat the shared epoch the same way, which keeps TEC continuous in time.
"""
struct GlobalMap{I, M}
    ionex::I
    mapping::M
end
GlobalMap(paths; mapping) = GlobalMap(foldl(append_ionex, map(read_ionex, paths)), mapping)
GlobalMap(path::AbstractString; mapping) = GlobalMap((path,); mapping)

# `b`'s maps from its first epoch on replace `a`'s
function append_ionex(a, b)
    for k in (:lats, :lons, :radius, :height, :mapping)
        a[k] == b[k] || throw(ArgumentError("IONEX files differ in $k: $(a[k]) vs $(b[k])"))
    end
    first(a.epochs) < first(b.epochs) && last(a.epochs) < last(b.epochs) ||
        throw(ArgumentError("IONEX files not in time order: $(first(a.epochs))–$(last(a.epochs)) followed by $(first(b.epochs))–$(last(b.epochs))"))
    # interpolating across the junction is no coarser than within either file
    spacing = maximum(I -> maximum(diff(I.epochs); init=Millisecond(0)), (a, b))
    first(b.epochs) - last(a.epochs) <= spacing ||
        throw(ArgumentError("gap between IONEX files: $(last(a.epochs)) to $(first(b.epochs)), maps are $spacing apart"))
    keep = a.epochs .< first(b.epochs)
    (; a..., epochs=[a.epochs[keep]; b.epochs], tec=cat(a.tec[:, :, keep], b.tec; dims=3))
end

function bilinear(I, M, lat, lon)
    lon = mod(lon + 180, 360) - 180
    fi = (lat - first(I.lats)) / step(I.lats) + 1; fj = (lon - first(I.lons)) / step(I.lons) + 1
    i0 = clamp(floor(Int, fi), 1, length(I.lats) - 1); j0 = clamp(floor(Int, fj), 1, length(I.lons) - 1)
    p = fi - i0; q = fj - j0
    (1-p)*(1-q)*M[i0,j0] + p*(1-q)*M[i0+1,j0] + (1-p)*q*M[i0,j0+1] + p*q*M[i0+1,j0+1]
end

# lat, lon in degrees (geocentric)
function vtec(I, t, lat, lon)
    I.epochs[1] <= t <= I.epochs[end] || error("$t outside IONEX map epochs $(I.epochs[1])–$(I.epochs[end])")
    k = clamp(searchsortedlast(I.epochs, t), 1, length(I.epochs) - 1)
    dt0 = Dates.value(t - I.epochs[k]) / 1e3; dt1 = Dates.value(t - I.epochs[k+1]) / 1e3
    Δ = Dates.value(I.epochs[k+1] - I.epochs[k]) / 1e3
    e0 = bilinear(I, view(I.tec, :, :, k), lat, lon + dt0 * 360 / 86400)
    e1 = bilinear(I, view(I.tec, :, :, k+1), lat, lon + dt1 * 360 / 86400)
    ((Δ - dt0) * e0 + dt0 * e1) / Δ
end

mapping_factor(::SLM, pp, k, site) = norm(pp) / dot(pp, k)
mapping_factor(::MSLM, pp, k, site) = 1 / sqrt(1 - (6371 / (6371 + 506.7) * sin(0.9782 * acos(dot(k, up(site)))))^2)

function slant_tec(m::GlobalMap, o, site, k, t)
    I = m.ionex
    pp = pierce(o, k, I.radius + I.height)
    (; lat, lon) = geocentric_latlon(pp)
    (; stec=vtec(I, t, rad2deg(lat), rad2deg(lon)) * mapping_factor(m.mapping, pp, k, site), pp)
end
