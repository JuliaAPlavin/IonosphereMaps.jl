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
A map node without a value (9999 in IONEX) makes interpolation through it NaN.

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
    (; a..., epochs=[a.epochs[keep]; b.epochs], tec=catmaps(a.tec, b.tec, keep), rms=catmaps(a.rms, b.rms, keep))
end
catmaps(A, B, keep) = cat(A[:, :, keep], B; dims=3)
catmaps(::Nothing, ::Nothing, keep) = nothing

# grid cell of x in the range r: x = r[i] + p * step(r), 0 ≤ p ≤ 1
function gridcell(r, x)
    f = (x - first(r)) / step(r)
    0 <= f <= length(r) - 1 || throw(DomainError(x, "outside the IONEX grid $(first(r)) to $(last(r))"))
    i = min(floor(Int, f), length(r) - 2)
    (i + 1, f - i)
end

function bilinear(I, M, lat, lon)
    (i, p) = gridcell(I.lats, lat)
    # a full circle without a repeated end column closes on its first column;
    # either side of a repeated (e.g. ±180) column interpolates from its own copy
    closed = length(I.lons) * step(I.lons) ≈ 360
    lons = range(first(I.lons); step=step(I.lons), length=length(I.lons) + closed)
    (j, q) = gridcell(lons, first(lons) + mod(lon - first(lons), 360))
    j1 = mod1(j + 1, length(I.lons))
    (1-p)*(1-q)*M[i,j] + p*(1-q)*M[i+1,j] + (1-p)*q*M[i,j1] + p*q*M[i+1,j1]
end

# maps A[lat, lon, epoch] at t; lat, lon in degrees (geocentric)
function interpolate(I, A, t, lat, lon)
    k = searchsortedlast(I.epochs, t)
    t == I.epochs[k] && return bilinear(I, view(A, :, :, k), lat, lon)
    w = (t - I.epochs[k]) /ₜ (I.epochs[k+1] - I.epochs[k])
    e0 = bilinear(I, view(A, :, :, k), lat, lon + 360 * ((t - I.epochs[k]) /ₜ Day))
    e1 = bilinear(I, view(A, :, :, k+1), lat, lon + 360 * ((t - I.epochs[k+1]) /ₜ Day))
    (1 - w) * e0 + w * e1
end
interpolate(I, ::Nothing, t, lat, lon) = nothing

# VTEC (TECU), ± its RMS where the file has RMS maps. RMS is interpolated with the TEC weights:
# errors at neighbouring nodes of a smooth map are close to fully correlated.
vtec(I, t, lat, lon) = with_rms(interpolate(I, I.tec, t, lat, lon), interpolate(I, I.rms, t, lat, lon))
with_rms(v, σ) = v ±ᵤ σ
with_rms(v, ::Nothing) = v

# below the horizon: NaN of the type a computed RM has
nan_rm(m::GlobalMap) = with_rms(NaN, isnothing(m.ionex.rms) ? nothing : NaN)

mapping_factor(::SLM, pp, k, site) = norm(pp) / dot(pp, k)
mapping_factor(::MSLM, pp, k, site) = 1 / sqrt(1 - (6371 / (6371 + 506.7) * sin(0.9782 * acos(dot(k, up(site)))))^2)

function slant_tec(m::GlobalMap, o, site, k, t)
    I = m.ionex
    pp = pierce(o, k, I.radius + I.height)
    (; lat, lon) = geocentric_latlon(pp)
    (; stec=vtec(I, t, rad2deg(lat), rad2deg(lon)) * mapping_factor(m.mapping, pp, k, site), pp)
end
