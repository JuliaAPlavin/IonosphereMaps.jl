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
using Unitful: @u_str, NoUnits, uconvert, ustrip

export ionospheric_rm, ionospheric_tec, LineOfSight, Vertical, GlobalMap, SLM, MSLM

"""
    LineOfSight(site, source)

From `site`, geodetic WGS84 `(; lat, lon, h)` with angles (Unitful, or plain radians) and a Unitful height,
toward `source` in any SkyCoords frame.
"""
struct LineOfSight{S, C<:AbstractSkyCoords}
    site::S
    source::C
end

"""
    Vertical(lat, lon)

Vertical column through geocentric `lat`, `lon` on the ionospheric shell; angles Unitful, or plain radians.
"""
struct Vertical{A, B}
    lat::A
    lon::B
end

include("geometry.jl")
include("io.jl")
include("globalmap.jl")

"""
    ionospheric_tec(model, geometry, t::DateTime)

Ionospheric total electron content (m⁻²) at UTC `t`: vertical for `Vertical`, slant for `LineOfSight`
(NaN below the horizon). For maps with RMS maps, returns `tec ±ᵤ σ` (`Uncertain.Value`), σ from the map RMS.
"""
function ionospheric_tec(model, los::LineOfSight, t::DateTime)
    k = source_direction(convert(ICRSCoords, los.source), t)
    dot(k, up(los.site)) > 0 || return nan_like(model, TECU)
    slant_tec(model, los.site, k, t).stec * TECU
end

"""
    ionospheric_rm(model, los::LineOfSight, t::DateTime)

Ionospheric Faraday rotation measure (rad/m²) along `los` at UTC `t`; positive when the field along the line of sight
points toward the observer, NaN below the horizon. For maps with RMS maps, returns `rm ±ᵤ σ` (`Uncertain.Value`);
σ covers the map RMS only, not the thin-shell and mapping approximations or the field model.
"""
function ionospheric_rm(model, los::LineOfSight, t::DateTime)
    k = source_direction(convert(ICRSCoords, los.source), t)
    dot(k, up(los.site)) > 0 || return nan_like(model, u"rad/m^2")
    (; stec, pp) = slant_tec(model, los.site, k, t)
    uconvert(u"rad/m^2", K_RM * stec * TECU * bpar(pp, k, t))
end

end
