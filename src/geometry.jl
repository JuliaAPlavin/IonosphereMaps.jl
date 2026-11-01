# RM = K_RM ∫ nₑ B∥ dl
const K_RM = u"q"^3 / (8π^2 * u"ε0" * u"me"^2 * u"c"^3)
const TECU = 1e16u"m^-2"

const WGS84 = (; a=6378137.0u"m", f=1 / 298.257223563)

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

# ICRS direction → Earth-fixed unit vector: precession, nutation, Earth rotation with UT1 ≈ UTC;
# no aberration (≤20″, far below map errors in RM), polar motion or refraction
source_direction(src::ICRSCoords, t) = r_eci_to_ecef(J2000(), PEF(), julian_day(t)) * cartesian(src).vec

# IGRF field component along -k, i.e. toward the observer, at ECEF position x
function bpar(x, k, t)
    (; lat, lon) = geocentric_latlon(x)
    b = igrf(yeardecimal(t), ustrip(u"m", norm(x)), lat, lon) * u"nT"   # geocentric (north, east, down)
    north = SVector(-sin(lat) * cos(lon), -sin(lat) * sin(lon), cos(lat))
    east = SVector(-sin(lon), cos(lon), 0.0)
    down = -SVector(cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat))
    dot(b[1] * north + b[2] * east + b[3] * down, -k)
end
