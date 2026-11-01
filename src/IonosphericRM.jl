module IonosphericRM

using Dates
using AxisKeys: KeyedArray, axiskeys
using DateFormats: /ₜ, julian_day, yeardecimal
using LinearAlgebra
using StaticArrays
using SkyCoords: AbstractSkyCoords, ICRSCoords, cartesian
using SatelliteToolboxGeomagneticField: igrf
using SatelliteToolboxTransformations: r_eci_to_ecef, J2000, PEF
using Uncertain: ±ᵤ
using Unitful: @u_str, uconvert, ustrip

export ionospheric_rm, GlobalMap, SLM, MSLM

include("geometry.jl")
include("io.jl")
include("globalmap.jl")

"""
    ionospheric_rm(model, site, source::AbstractSkyCoords, t::DateTime)

Ionospheric Faraday rotation measure (rad/m²) toward `source` (any SkyCoords frame) seen from `site` at UTC `t`.
`site` is geodetic WGS84 `(; lat, lon, h)`: angles (Unitful, or plain radians) and a Unitful height.
Returns a Unitful RM; positive when the field along the line of sight points toward the observer, NaN below the horizon.

For maps with RMS maps, returns `rm ±ᵤ σ` (`Uncertain.Value`); σ covers the map RMS only,
not the thin-shell and mapping approximations or the field model.
"""
function ionospheric_rm(model, site, source::AbstractSkyCoords, t::DateTime)
    k = source_direction(convert(ICRSCoords, source), t)
    dot(k, up(site)) > 0 || return nan_rm(model)
    (; stec, pp) = slant_tec(model, site, k, t)
    uconvert(u"rad/m^2", K_RM * stec * TECU * bpar(pp, k, t))
end

end
