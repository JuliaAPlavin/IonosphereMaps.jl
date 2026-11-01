module IonosphericRM

using Dates
using LinearAlgebra
using StaticArrays
using SkyCoords: ICRSCoords
using SatelliteToolboxGeomagneticField: igrf

export ionospheric_rm, GlobalMap, SLM, MSLM

include("geometry.jl")
include("io.jl")
include("globalmap.jl")

"""
    ionospheric_rm(model, site, source::ICRSCoords, t::DateTime)

Ionospheric Faraday rotation measure (rad/m²) toward `source` (J2000/ICRS) seen from `site` at UTC `t`.
`site` is geodetic WGS84: `(; lat, lon, h)` in radians and metres.
Positive RM: field component along the line of sight points toward the observer.
"""
function ionospheric_rm(model, site, source::ICRSCoords, t::DateTime)
    o = ecef(site)
    k = source_direction(source, t)
    dot(k, up(site)) > 0 || return NaN   # below the horizon
    (; stec, pp) = slant_tec(model, o, site, k, t)
    K_RM * stec * 1e16 * bpar(pp, k, t)
end

end
