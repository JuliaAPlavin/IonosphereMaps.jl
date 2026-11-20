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
`mapping` (`SLM()`, `MSLM()` or `Petrov2023()`) converts VTEC to slant TEC; [`vtec`](@ref) looks up the VTEC itself.
The files' MAPPING FUNCTION record (COSZ or NONE in practice) is used only to check that joined files agree.
Supported are 2-D TEC maps with ascending longitudes and their RMS maps, if any, paired by map number;
height maps, 3-D maps or descending longitudes are an error.
CODE produces its maps with MSLM; Petrov2023 is fitted to VLBI astrometry and includes its scale factor k = 0.85.
Points beyond the grid (e.g. outside a regional map) and times outside the maps are an error, except the polar caps
of a map covering all longitudes whose outermost row lies one grid step from the pole (e.g. 87.5° for a 2.5° step):
there the pole is a node valued the mean of that row, and a point in the cap blends linearly in latitude between the row,
interpolated in longitude, and the pole.
For a regional map, the rotation with the Sun can move a point inside the grid outside it near the east and west edges
between map epochs: also an error.
A missing map node (9999 or `nan` in IONEX) gives NaN wherever it enters the interpolation.

`paths` are files of consecutive periods in time order, e.g. the daily files of an observation crossing
midnight. They must share grid, base radius, shell height and mapping function, all or none have RMS maps,
and gaps are at most one map interval.
Where files overlap (normally the midnight map of two daily solutions), the later file's maps are used:
IONEX maps cover `Tᵢ ≤ t < Tᵢ₊₁`, so a day's 24:00 map only closes its last interval. This keeps TEC continuous.
"""
struct GlobalMap{I, M}
    ionex::I
    mapping::M
end
GlobalMap(paths; mapping) = GlobalMap(foldl(append_ionex, map(read_ionex, paths)), mapping)
GlobalMap(path::AbstractString; mapping) = GlobalMap((path,); mapping)

# `b`'s maps from its first epoch on replace `a`'s: at midnight, the next day's 00:00 map replaces the 24:00 one, as in AIPS TECOR, CASA and spinifex
function append_ionex(a, b)
    for (k, x, y) in [(:lat, axiskeys(a.tec, :lat), axiskeys(b.tec, :lat)), (:lon, axiskeys(a.tec, :lon), axiskeys(b.tec, :lon)),
                      (:radius, a.radius, b.radius), (:height, a.height, b.height), (:mapping_function, a.mapping_function, b.mapping_function)]
        x == y || throw(ArgumentError("IONEX files differ in $k: $x vs $y"))
    end
    isnothing(a.rms) == isnothing(b.rms) || throw(ArgumentError("IONEX files differ in RMS maps: present in one, absent in the other"))
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
    mlats, mlons = axiskeys(M, :lat), axiskeys(M, :lon)
    # a full circle without a repeated end column closes on its first column;
    # either side of a repeated (e.g. ±180) column interpolates from its own copy
    closed = length(mlons) * step(mlons) ≈ 360u"°"
    lons = range(first(mlons); step=step(mlons), length=length(mlons) + closed)
    (j, q) = gridcell(lons, first(lons) + mod(lon - first(lons), 360u"°"))
    j1 = mod1(j + 1, length(mlons))
    # the columns of distinct longitudes
    ring = 1:length(mlons) - !closed
    # in a map covering all longitudes, a pole one step beyond an outermost row is a node too: the mean of that row
    fullring = length(ring) * step(mlons) ≈ 360u"°"
    before, after = map(x -> fullring && abs(x) ≈ 90u"°", (first(mlats) - step(mlats), last(mlats) + step(mlats)))
    lats = range(first(mlats) - before * step(mlats); step=step(mlats), length=length(mlats) + before + after)
    (i, p) = gridcell(lats, lat)
    n = size(M, 1)
    # row r of M interpolated in longitude; r = 0 or n + 1 is the pole next to row 1 or n
    row(r) = 1 <= r <= n ? (1-q)*M[r,j] + q*M[r,j1] : sum(view(M, clamp(r, 1, n), ring)) / length(ring)
    (1-p)*row(i - before) + p*row(i + 1 - before)
end

# index k of the map interval Tₖ ≤ t < Tₖ₊₁ of the maps I, or of the last map at its epoch
function mapindex(I, t)
    epochs = axiskeys(I.tec, :epoch)
    first(epochs) <= t <= last(epochs) || throw(DomainError(t, "outside the IONEX maps $(first(epochs)) to $(last(epochs))"))
    searchsortedlast(epochs, t)
end

# maps A[lat, lon, epoch] at lat, lon in degrees (geocentric) and t in the map interval k
function interpolate(A, lat, lon, t, k)
    epochs = axiskeys(A, :epoch)
    t == epochs[k] && return bilinear(view(A, :, :, k), lat, lon)
    w = (t - epochs[k]) / (epochs[k+1] - epochs[k])
    e0 = bilinear(view(A, :, :, k), lat, lon + 360u"°" * ((t - epochs[k]) / Day(1)))
    e1 = bilinear(view(A, :, :, k+1), lat, lon + 360u"°" * ((t - epochs[k+1]) / Day(1)))
    (1 - w) * e0 + w * e1
end
interpolate(::Nothing, lat, lon, t, k) = nothing

"""
    vtec(m::GlobalMap, p, t::DateTime)

Vertical total electron content (m⁻²) of the maps of `m` at the point `p` of their shell, geocentric `(; lat, lon)`
or any object with these properties (angles Unitful, or plain radians), at UTC `t`, interpolated as described in [`GlobalMap`](@ref). For maps with RMS maps, returns
`tec ±ᵤ σ` (`Uncertain.Value`), σ from the map RMS. This is the maps' own data, not a model prediction:
the mapping of `m` plays no role (e.g. no scale factor k of `Petrov2023`).
"""
vtec(m::GlobalMap, p, t::DateTime) = vtec_in_interval(m, p, t, mapindex(m.ionex, t))
# VTEC at t in the map interval k (from mapindex)
# RMS is interpolated with the TEC weights: errors at neighbouring nodes of a smooth map are close to fully correlated.
function vtec_in_interval(m::GlobalMap, p, t, k)
    (; tec, rms) = m.ionex
    lat, lon = uconvert(u"°", p.lat), uconvert(u"°", p.lon)
    with_rms(interpolate(tec, lat, lon, t, k), interpolate(rms, lat, lon, t, k)) * TECU
end
with_rms(v, σ) = v ±ᵤ σ
with_rms(v, ::Nothing) = v

# below the horizon: NaN of the type VTEC (TECU) has with these RMS maps
nan_like(::Nothing) = NaN
nan_like(rms::AbstractArray) = NaN ±ᵤ NaN

# mapping from VTEC to slant TEC for `ray`, piercing the shell of the maps `I` at pp, `up` the geodetic up at its origin;
# angles from the vertical as atan(sin, cos): accurate near the zenith, not limited to an exactly unit direction
mapping_factor(::SLM, I, pp, ray, up) = norm(pp) / dot(pp, ray.direction)
function mapping_factor(::MSLM, I, pp, ray, up)
    R, H, α = 6371u"km", 506.7u"km", 0.9782
    z = atan(norm(cross(ray.direction, up)), dot(ray.direction, up))
    1 / sqrt(1 - (R / (R + H) * sin(α * z))^2)
end
function mapping_factor(::Petrov2023, I, pp, ray, up)
    k, ΔH, α = 0.85, 56.7u"km", 0.9782
    R, H = I.radius, I.height
    e = atan(dot(ray.direction, ray.origin), norm(cross(ray.direction, ray.origin)))
    k / sqrt(1 - (R / (R + H + ΔH) * cos(α * e))^2)
end

# slant TEC (m⁻²) along `ray` and its pierce point pp; the time and the ray origin are checked also below the horizon
function slant_tec(m::GlobalMap, ray, t)
    I = m.ionex
    k = mapindex(I, t)
    pp = pierce(ray.origin, ray.direction, I.radius + I.height)
    up = geodetic_up(ray.origin)
    # NaN below the geodetic horizon of the ray origin
    stec = dot(ray.direction, up) > 0 ? vtec_in_interval(m, geocentric_latlon(pp), t, k) * mapping_factor(m.mapping, I, pp, ray, up) :
           nan_like(I.rms) * TECU
    (; stec, pp)
end
