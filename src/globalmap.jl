"""Single-layer mapping at the IONEX shell height: exact 1/cos z′ at the pierce point."""
struct SLM end
"""CODE modified single-layer mapping (H = 506.7 km, α = 0.9782, R = 6371 km) at the geodetic zenith distance."""
struct MSLM end

"""
    GlobalMap(path; mapping=SLM())

Thin-shell model from an IONEX global ionosphere map: VTEC bilinearly interpolated at the pierce point,
maps blended in time after rotating each with the Sun (IONEX recommendation).
"""
struct GlobalMap{I, M}
    ionex::I
    mapping::M
end
GlobalMap(path::AbstractString; mapping=SLM()) = GlobalMap(read_ionex(path), mapping)

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
