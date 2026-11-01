module IonosphericRM

using Dates
using LinearAlgebra
using StaticArrays
using StructArrays
using DataManipulation
using CodecZlib
using SkyCoords: ICRSCoords
using SatelliteToolboxGeomagneticField: igrf

export tec_model, ionospheric_rm, GlobalMap, SingleStation, SLM, MSLM, MeanVTEC, PlanarVTEC

include("geometry.jl")
include("io.jl")
include("globalmap.jl")
include("singlestation.jl")

"""
    tec_model(path; kw...)

Load a TEC source, dispatching on the file content:
- IONEX file → [`GlobalMap`](@ref);
- RINEX 3 observation file, or a directory containing exactly one → [`SingleStation`](@ref)
  (SP3 orbits and an IONEX file with the DCB block are taken from the same directory).

Files may be gzip-compressed. Keywords are passed to the model constructor.
"""
function tec_model(path; kw...)
    isdir(path) && return SingleStation(only(filter(isrinexobs, readdir(path; join=true))); kw...)
    isionex(path) ? GlobalMap(path; kw...) :
    isrinexobs(path) ? SingleStation(path; kw...) :
    error("$path is neither an IONEX nor a RINEX observation file")
end

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
