# RM = K_RM ∫ nₑ B∥ dl
const K_RM = u"q"^3 / (8π^2 * u"ε0" * u"me"^2 * u"c"^3)
# group delay = K_ION ∫ nₑ dl / (c ν²)
const K_ION = u"q"^2 / (8π^2 * u"ε0" * u"me")
const TECU = 1e16u"m^-2"

# geodetic WGS84 site → ECEF; the library takes plain metres
ecef(site) = geodetic_to_ecef(site.lat, site.lon, ustrip(u"m", site.h)) * u"m"

# WGS84 ellipsoid normal through the ECEF position x
function geodetic_up(x)
    (lat, lon) = ecef_to_geodetic(ustrip.(u"m", x))
    SVector(cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat))
end

geocentric_latlon(x) = (; lat=asin(x[3] / norm(x)), lon=atan(x[2], x[1]))

# point where the ray o + s·k (s > 0) from inside the sphere |x| = r crosses it
function pierce(o, k, r)
    norm(o) < r || throw(DomainError(norm(o), "ray origin distance from the geocentre, at or above the ionospheric shell of radius $r"))
    b = dot(o, k)
    o + (-b + sqrt(b^2 - (dot(o, o) - r^2))) * k
end

# ICRS direction → Earth-fixed unit vector: precession, nutation, Earth rotation with UT1 ≈ UTC;
# no aberration (≤20″, far below map errors in RM), polar motion or refraction
source_direction(src::ICRSCoords, t) = r_eci_to_ecef(J2000(), PEF(), julian_day(t)) * cartesian(src).vec

# IGRF field component along -k, i.e. toward the observer, at ECEF position x
function bpar(x, k, t)
    (; lat, lon) = geocentric_latlon(x)
    b = igrf(yeardecimal(t), ustrip(u"m", norm(x)), lat, lon) * u"nT"   # geocentric (north, east, down)
    north = SVector(-sin(lat) * cos(lon), -sin(lat) * sin(lon), cos(lat))
    east = SVector(-sin(lon), cos(lon), zero(lon))
    down = -SVector(cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat))
    dot(b[1] * north + b[2] * east + b[3] * down, -k)
end
