"""Single-layer mapping at the IONEX shell height: exact 1/cos z′ at the pierce point."""
struct SLM end
"""
CODE modified single-layer mapping `1 / √(1 − (R/(R + H))² sin²(α z))`, H = 506.7 km, α = 0.9782, R = 6371 km,
z the geodetic zenith distance at the ray origin.
"""
struct MSLM end
"""
Petrov (2023, AJ 165, 183, eq. 12) mapping `k / √(1 − (R/(R + H + ΔH))² cos²(α e))`, k = 0.85, ΔH = 56.7 km,
α = 0.9782; R and H the map's base radius and shell height, e the geocentric elevation at the ray origin.
"""
struct Petrov2023 end

"""
    GlobalMap(paths; mapping)
    GlobalMap(path; mapping)

Thin-shell model from IONEX global ionosphere maps: VTEC bilinearly interpolated at the pierce point,
maps blended in time after rotating each with the Sun (IONEX method 3), also across file boundaries.
`mapping` (`SLM()`, `MSLM()` or `Petrov2023()`) converts VTEC to slant TEC. The files' MAPPING FUNCTION record
(COSZ or NONE in practice) is used only to check that joined files agree. CODE produces its maps with MSLM;
Petrov2023 is fitted to VLBI astrometry and includes its scale factor k = 0.85.
Pierce points beyond the grid (e.g. past ±87.5° latitude, or outside a regional map), and times outside the maps, are an error.
For a regional map, the rotation with the Sun can move a point inside the grid outside it near the east and west edges
between map epochs: also an error.
A missing map node (9999 or `nan` in IONEX) gives NaN throughout its grid cells.

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
        throw(ArgumentError("gap between IONEX files: $(last(ea)) to $(first(eb)), longer than any map interval within them"))
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
    closed = length(mlons) * step(mlons) ≈ 360u"°"
    lons = range(first(mlons); step=step(mlons), length=length(mlons) + closed)
    (j, q) = gridcell(lons, first(lons) + mod(lon - first(lons), 360u"°"))
    j1 = mod1(j + 1, length(mlons))
    (1-p)*(1-q)*M[i,j] + p*(1-q)*M[i+1,j] + (1-p)*q*M[i,j1] + p*q*M[i+1,j1]
end

# index k of the map interval Tₖ ≤ t < Tₖ₊₁ of the maps I, or of the last map at its epoch
function mapindex(I, t)
    epochs = axiskeys(I.tec, :epoch)
    first(epochs) <= t <= last(epochs) || throw(DomainError(t, "outside the IONEX maps $(first(epochs)) to $(last(epochs))"))
    searchsortedlast(epochs, t)
end

# maps A[lat, lon, epoch] at t in the interval k; lat, lon in degrees (geocentric)
function interpolate(A, t, k, lat, lon)
    epochs = axiskeys(A, :epoch)
    t == epochs[k] && return bilinear(view(A, :, :, k), lat, lon)
    w = (t - epochs[k]) /ₜ (epochs[k+1] - epochs[k])
    e0 = bilinear(view(A, :, :, k), lat, lon + 360u"°" * ((t - epochs[k]) /ₜ Day))
    e1 = bilinear(view(A, :, :, k+1), lat, lon + 360u"°" * ((t - epochs[k+1]) /ₜ Day))
    (1 - w) * e0 + w * e1
end
interpolate(::Nothing, t, k, lat, lon) = nothing

# VTEC (TECU), ± its RMS where the file has RMS maps, at geocentric lat, lon (Unitful, or plain radians).
# RMS is interpolated with the TEC weights: errors at neighbouring nodes of a smooth map are close to fully correlated.
vtec(I, t, lat, lon) = vtec(I, t, mapindex(I, t), lat, lon)
function vtec(I, t, k, lat, lon)
    lat, lon = uconvert(u"°", lat), uconvert(u"°", lon)
    with_rms(interpolate(I.tec, t, k, lat, lon), interpolate(I.rms, t, k, lat, lon))
end
with_rms(v, σ) = v ±ᵤ σ
with_rms(v, ::Nothing) = v

# below the horizon: NaN of the type VTEC has with these RMS maps
nan_like(::Nothing) = NaN
nan_like(rms::AbstractArray) = NaN ±ᵤ NaN

ionospheric_tec(m::GlobalMap, v::Vertical, t::DateTime) = vtec(m.ionex, t, v.lat, v.lon) * TECU

# mapping from VTEC to slant TEC for `ray`, piercing the shell of the maps `I` at pp
mapping_factor(::SLM, I, pp, ray) = norm(pp) / dot(pp, ray.direction)
function mapping_factor(::MSLM, I, pp, ray)
    R, H, α = 6371u"km", 506.7u"km", 0.9782
    z = acos(dot(ray.direction, geodetic_up(ray.origin)))
    1 / sqrt(1 - (R / (R + H) * sin(α * z))^2)
end
function mapping_factor(::Petrov2023, I, pp, ray)
    k, ΔH, α = 0.85, 56.7u"km", 0.9782
    R, H = I.radius, I.height
    e = asin(dot(ray.direction, ray.origin) / norm(ray.origin))
    k / sqrt(1 - (R / (R + H + ΔH) * cos(α * e))^2)
end

# slant TEC (TECU) along `ray` and its pierce point pp; the time and the ray origin are checked also below the horizon
function slant_tec(m::GlobalMap, ray, t)
    I = m.ionex
    k = mapindex(I, t)
    pp = pierce(ray.origin, ray.direction, I.radius + I.height)
    stec = above_horizon(ray) ? vtec(I, t, k, geocentric_latlon(pp)...) * mapping_factor(m.mapping, I, pp, ray) : nan_like(I.rms)
    (; stec, pp)
end
