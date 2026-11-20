# IonosphereMaps.jl

Ionospheric slant TEC, Faraday rotation measure and group delay along a line of sight, and the maps' own vertical TEC,
from IONEX global ionosphere maps (GIMs) in the thin-shell model.

```julia
using IonosphereMaps, Dates, SkyCoords, Unitful

m = GlobalMap(["COD0OPSFIN_20250260000_01D_01H_GIM.INX.gz", "COD0OPSFIN_20250270000_01D_01H_GIM.INX.gz"]; mapping=MSLM())
site = (; lat=-30.71292u"°", lon=21.44380u"°", h=1050u"m")                 # MeerKAT
los = LineOfSight(site, ICRSCoords(deg2rad(80.29119), deg2rad(16.63946)))  # 3C138
t = DateTime(2025, 1, 26, 23, 30)

ionospheric_tec(m, los, t)                      # 5.39e17 m^-2 ±ᵤ 5.81e16 m^-2
ionospheric_rm(m, los, t)                       # -2.257 rad m^-2 ±ᵤ 0.243 rad m^-2
ionospheric_group_delay(m, los, t, 8.4u"GHz")   # 1.026e-9 s ±ᵤ 1.11e-10 s
vtec(m, (; lat=-30u"°", lon=21u"°"), t)         # 1.639e17 m^-2 ±ᵤ 2.21e16 m^-2
```

## Model

`GlobalMap(paths; mapping)` reads IONEX files (plain, `.gz`, or `.Z`; decompressed by Gzip_jll, no system `gzip` needed):
2-D TEC maps with ascending longitudes, and their RMS maps if present.
VTEC at the pierce point of the shell is interpolated bilinearly in each map, and in time between maps rotated with the Sun (IONEX method 3).
Slant TEC is VTEC times the mapping factor; RM is `K·TEC·B∥` with the IGRF field at the pierce point;
group delay is `K·TEC / (c ν²)` at the cyclic frequency `ν`.

Geometries:
- `LineOfSight(site, source)`: a geodetic WGS84 site `(; lat, lon, h)` toward a source in any SkyCoords frame.
- `Ray(origin, direction)`: an ECEF origin (Unitful lengths) inside the shell and an Earth-fixed unit direction. `Ray(los, t)` is the ray of `los` at time `t`:
  precession, nutation and Earth rotation with UT1 ≈ UTC; no aberration, polar motion or refraction.

Mappings:
- `SLM()`: the IONEX single layer, exact `1/cos z′` at the pierce point on the map's shell.
- `MSLM()`: CODE's modified single-layer mapping (506.7 km, α = 0.9782), the one CODE produces its maps with.
- `Petrov2023()`: Petrov (2023, AJ 165, 183) mapping fitted to VLBI astrometry, its TEC scale factor k = 0.85 included.

`vtec(m, p, t)` is the maps' own VTEC at a point `p = (; lat, lon)` on the shell (geocentric); the mapping plays no role.

## Uncertainties

For files with RMS maps, every result is `value ±ᵤ σ` (`Uncertain.Value`); without RMS maps, a plain Unitful quantity.
`Uncertain.value` returns the central value of either.
The RMS is interpolated with the same weights as the TEC (errors of neighbouring nodes taken as fully correlated) and scaled by the same factors.
σ covers the map RMS only: not the thin-shell and mapping approximations or the field model.

## Behaviour

- RMS map i belongs to TEC map i, by the map number in its `START OF RMS MAP` record; its epoch record is optional and, if present, must be the TEC map's.
- Height maps, 3-D maps and descending longitudes (negative DLON) are not supported: an error (`ArgumentError`).
- Several files of consecutive periods (e.g. the daily files around midnight) join into one model: same grid, shell and mapping function,
  RMS maps in all or none, gaps of at most one map interval.
  Where they overlap, the later file's maps are used, so TEC is continuous.
- A time outside the maps, a point outside the grid, and a ray origin at or above the shell are errors (`DomainError`).
- Below the geodetic horizon of the ray origin, results are NaN.
- In a map covering all longitudes whose outermost row is one grid step from the pole, the pole is a node valued the mean of that row: pierce points in the polar caps are interpolated, not errors.
- A missing map value (9999 or `nan`) gives NaN wherever it enters the interpolation.
