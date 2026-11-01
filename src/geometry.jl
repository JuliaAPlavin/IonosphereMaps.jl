# RM = K_RM ∫ nₑ B∥ dl  (SI), K_RM = e³ / (8π² ε₀ mₑ² c³)
const K_RM = 1.602176634e-19^3 / (8π^2 * 8.8541878128e-12 * 9.1093837015e-31^2 * 299792458.0^3)

const WGS84 = (; a=6378137.0, f=1 / 298.257223563)

function ecef(site)
    (; lat, lon, h) = site
    e2 = WGS84.f * (2 - WGS84.f)
    N = WGS84.a / sqrt(1 - e2 * sin(lat)^2)
    SVector((N + h) * cos(lat) * cos(lon), (N + h) * cos(lat) * sin(lon), (N * (1 - e2) + h) * sin(lat))
end

up(site) = SVector(cos(site.lat) * cos(site.lon), cos(site.lat) * sin(site.lon), sin(site.lat))

geocentric_latlon(x) = (; lat=asin(x[3] / norm(x)), lon=atan(x[2], x[1]))

# point where the ray o + s·k (s > 0) crosses the sphere |x| = r
function pierce(o, k, r)
    b = dot(o, k)
    o + (-b + sqrt(b^2 - (dot(o, o) - r^2))) * k
end

# ICRS direction → Earth-fixed unit vector: precession, nutation, Earth rotation with UT1 ≈ UTC; no polar motion
source_direction(src::ICRSCoords, t) = r_eci_to_ecef(J2000(), PEF(), julian_day(t)) * cartesian(src).vec

# IGRF field component (T) along -k, i.e. toward the observer, at ECEF position x
function bpar(x, k, t)
    (; lat, lon) = geocentric_latlon(x)
    b = igrf(yeardecimal(t), norm(x), lat, lon)   # geocentric (north, east, down), nT
    north = SVector(-sin(lat) * cos(lon), -sin(lat) * sin(lon), cos(lat))
    east = SVector(-sin(lon), cos(lon), 0.0)
    down = -SVector(cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat))
    1e-9 * dot(b[1] * north + b[2] * east + b[3] * down, -k)
end
