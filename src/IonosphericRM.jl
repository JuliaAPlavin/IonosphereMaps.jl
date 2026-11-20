module IonosphericRM

using Dates
using AxisKeys: KeyedArray, axiskeys
using DateFormats: /ₜ, julian_day, yeardecimal
using LinearAlgebra
using StaticArrays
using SkyCoords: AbstractSkyCoords, ICRSCoords, cartesian
using SatelliteToolboxGeomagneticField: igrf
using SatelliteToolboxTransformations: r_eci_to_ecef, J2000, PEF, ecef_to_geodetic, geodetic_to_ecef
using Uncertain: ±ᵤ
using Unitful: @u_str, Frequency, Length, NoUnits, uconvert, ustrip

export ionospheric_rm, ionospheric_tec, ionospheric_group_delay, Ray, LineOfSight, Vertical, GlobalMap, SLM, MSLM, Petrov2023

"""
    Ray(origin, direction)

From `origin`, an Earth-fixed (ECEF, ITRF/WGS84) position as an `SVector{3}` of Unitful lengths,
along `direction`, an Earth-fixed unit `SVector{3}`.

    Ray(los::LineOfSight, t::DateTime)

The ray of `los` at UTC `t`: from the site's ECEF position toward the source's Earth-fixed direction at `t`.
"""
struct Ray{O, D}
    origin::O
    direction::D

    function Ray(origin::SVector{3, <:Length}, direction::SVector{3, <:Real})
        @assert norm(direction) ≈ 1
        new{typeof(origin), typeof(direction)}(origin, direction)
    end
end

"""
    LineOfSight(site, source)

From `site`, geodetic WGS84 `(; lat, lon, h)` with angles (Unitful, or plain radians) and a Unitful height,
toward `source` in any SkyCoords frame.
"""
struct LineOfSight{S, C<:AbstractSkyCoords}
    site::S
    source::C
end

Ray(los::LineOfSight, t::DateTime) = Ray(ecef(los.site), source_direction(convert(ICRSCoords, los.source), t))

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

Ionospheric total electron content (m⁻²) at UTC `t`: vertical for `Vertical`, slant for `Ray` and `LineOfSight`
(NaN below the geodetic horizon of the ray origin). For maps with RMS maps, returns `tec ±ᵤ σ` (`Uncertain.Value`),
σ from the map RMS.
"""
function ionospheric_tec(model, ray::Ray, t::DateTime)
    above_horizon(ray) || return nan_like(model, TECU)
    slant_tec(model, ray, t).stec * TECU
end
ionospheric_tec(model, los::LineOfSight, t::DateTime) = ionospheric_tec(model, Ray(los, t), t)

"""
    ionospheric_group_delay(model, geometry, t::DateTime, ν)

Ionospheric group delay (s) at frequency `ν` (Unitful) for any geometry of [`ionospheric_tec`](@ref) at UTC `t`:
`K·TEC / (c ν²)`, `K = e² / (8π² ε₀ mₑ)`. Positive: the group is delayed, the phase advanced by the same amount
(phase delay `-K·TEC / (c ν²)`). NaN below the horizon, and `delay ±ᵤ σ` for maps with RMS maps, as for TEC.
"""
ionospheric_group_delay(model, geometry, t::DateTime, ν::Frequency) =
    uconvert(u"s", K_ION * ionospheric_tec(model, geometry, t) / (u"c" * ν^2))

"""
    ionospheric_rm(model, geometry, t::DateTime)

Ionospheric Faraday rotation measure (rad/m²) along a `Ray` or `LineOfSight` at UTC `t`; positive when the field
along the line of sight points toward the observer, NaN below the geodetic horizon of the ray origin.
For maps with RMS maps, returns `rm ±ᵤ σ` (`Uncertain.Value`); σ covers the map RMS only, not the thin-shell
and mapping approximations or the field model.
"""
function ionospheric_rm(model, ray::Ray, t::DateTime)
    above_horizon(ray) || return nan_like(model, u"rad/m^2")
    (; stec, pp) = slant_tec(model, ray, t)
    uconvert(u"rad/m^2", K_RM * stec * TECU * bpar(pp, ray.direction, t))
end
ionospheric_rm(model, los::LineOfSight, t::DateTime) = ionospheric_rm(model, Ray(los, t), t)

end
