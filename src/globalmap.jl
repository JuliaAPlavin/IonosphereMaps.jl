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
Pierce points beyond the grid (e.g. past ±87.5° latitude, or outside a regional map) are an error.
A missing map node (9999 in IONEX) gives NaN throughout its grid cells.

`paths` are files of consecutive periods in time order, e.g. the daily files of an observation crossing
midnight. They must share grid, base radius, shell height and mapping function, with gaps of at most one map interval.
Where files overlap (normally the midnight map of two daily solutions), the later file's maps are used:
IONEX maps cover `Tᵢ ≤ t < Tᵢ₊₁`, so a day's 24:00 map only closes its last interval. This keeps TEC continuous.
"""
struct GlobalMap{I, M}
    ionex::I
    mapping::M
end
GlobalMap(paths; mapping) = GlobalMap(foldl(append_ionex, map(read_ionex, paths)), mapping)
GlobalMap(path::AbstractString; mapping) = GlobalMap((path,); mapping)

# `b`'s maps from its first epoch on replace `a`'s
function append_ionex(a, b)
    for (k, x, y) in [(:lat, axiskeys(a.tec, :lat), axiskeys(b.tec, :lat)), (:lon, axiskeys(a.tec, :lon), axiskeys(b.tec, :lon)),
                      (:radius, a.radius, b.radius), (:height, a.height, b.height), (:mapping, a.mapping, b.mapping)]
        x == y || throw(ArgumentError("IONEX files differ in $k: $x vs $y"))
    end
    ea, eb = axiskeys(a.tec, :epoch), axiskeys(b.tec, :epoch)
    first(ea) < first(eb) && last(ea) < last(eb) ||
        throw(ArgumentError("IONEX files not in time order: $(first(ea))–$(last(ea)) followed by $(first(eb))–$(last(eb))"))
    # interpolating across the junction is no coarser than within either file
    spacing = maximum(e -> maximum(diff(e); init=Millisecond(0)), (ea, eb))
    first(eb) - last(ea) <= spacing ||
        throw(ArgumentError("gap between IONEX files: $(last(ea)) to $(first(eb)), maps are $spacing apart"))
    (; a..., tec=catmaps(a.tec, b.tec), rms=catmaps(a.rms, b.rms))
end
catmaps(A, B) = cat(A(epoch=<(first(axiskeys(B, :epoch)))), B; dims=3)
catmaps(::Nothing, ::Nothing) = nothing

# grid cell of x in the range r: x = r[i] + p * step(r), 0 ≤ p ≤ 1
function gridcell(r, x)
    f = (x - first(r)) / step(r)
    0 <= f <= length(r) - 1 || throw(DomainError(x, "outside the IONEX grid $(first(r)) to $(last(r))"))
    i = min(floor(Int, f), length(r) - 2)
    (i + 1, f - i)
end

# map M[lat, lon] at lat, lon
function bilinear(M, lat, lon)
    (i, p) = gridcell(axiskeys(M, :lat), lat)
    # a full circle without a repeated end column closes on its first column;
    # either side of a repeated (e.g. ±180) column interpolates from its own copy
    mlons = axiskeys(M, :lon)
    closed = length(mlons) * step(mlons) ≈ 360
    lons = range(first(mlons); step=step(mlons), length=length(mlons) + closed)
    (j, q) = gridcell(lons, first(lons) + mod(lon - first(lons), 360))
    j1 = mod1(j + 1, length(mlons))
    (1-p)*(1-q)*M[i,j] + p*(1-q)*M[i+1,j] + (1-p)*q*M[i,j1] + p*q*M[i+1,j1]
end

# maps A[lat, lon, epoch] at t; lat, lon in degrees (geocentric)
function interpolate(A, t, lat, lon)
    epochs = axiskeys(A, :epoch)
    k = searchsortedlast(epochs, t)
    t == epochs[k] && return bilinear(view(A, :, :, k), lat, lon)
    w = (t - epochs[k]) /ₜ (epochs[k+1] - epochs[k])
    e0 = bilinear(view(A, :, :, k), lat, lon + 360 * ((t - epochs[k]) /ₜ Day))
    e1 = bilinear(view(A, :, :, k+1), lat, lon + 360 * ((t - epochs[k+1]) /ₜ Day))
    (1 - w) * e0 + w * e1
end
interpolate(::Nothing, t, lat, lon) = nothing

# VTEC (TECU), ± its RMS where the file has RMS maps. RMS is interpolated with the TEC weights:
# errors at neighbouring nodes of a smooth map are close to fully correlated.
vtec(I, t, lat, lon) = with_rms(interpolate(I.tec, t, lat, lon), interpolate(I.rms, t, lat, lon))
with_rms(v, σ) = v ±ᵤ σ
with_rms(v, ::Nothing) = v

# below the horizon: NaN of the type a computed RM has
nan_rm(m::GlobalMap) = with_rms(NaN, isnothing(m.ionex.rms) ? nothing : NaN)

mapping_factor(::SLM, pp, k, site) = norm(pp) / dot(pp, k)
mapping_factor(::MSLM, pp, k, site) = 1 / sqrt(1 - (6371 / (6371 + 506.7) * sin(0.9782 * acos(dot(k, up(site)))))^2)

function slant_tec(m::GlobalMap, site, k, t)
    I = m.ionex
    pp = pierce(ecef(site), k, I.radius + I.height)
    (; lat, lon) = geocentric_latlon(pp)
    (; stec=vtec(I, t, rad2deg(lat), rad2deg(lon)) * mapping_factor(m.mapping, pp, k, site), pp)
end
