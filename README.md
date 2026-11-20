# IonosphereMaps.jl

Ionospheric slant TEC, Faraday rotation measure and group delay along a line of sight, and the maps' own vertical TEC,
from IONEX global ionosphere maps (GIMs) in the thin-shell model.

```julia
using IonosphereMaps, Dates, SkyCoords, Unitful

m = GlobalMap(["COD0OPSFIN_20250260000_01D_01H_GIM.INX.gz", "COD0OPSFIN_20250270000_01D_01H_GIM.INX.gz"]; mapping=MSLM())
site = (; lat=-30.71292u"°", lon=21.44380u"°", h=1050u"m")                 # MeerKAT
los = LineOfSight(site, ICRSCoords(80.29119u"°", 16.63946u"°"))  # 3C138
t = DateTime(2025, 1, 26, 23, 30)

ionospheric_tec(m, los, t)                      # 53.9 TECU ±ᵤ 5.81 TECU
ionospheric_rm(m, los, t)                       # -2.257 rad m^-2 ±ᵤ 0.243 rad m^-2
ionospheric_group_delay(m, los, t, 8.4u"GHz")   # 1.026e-9 s ±ᵤ 1.11e-10 s
vtec(m, (; lat=-30u"°", lon=21u"°"), t)         # 16.39 TECU ±ᵤ 2.21 TECU
```

## Model

`GlobalMap(paths; mapping)` reads IONEX files (plain, `.gz`, or `.Z`):
2-D TEC maps, and their RMS maps when present.
VTEC at the pierce point of the shell is interpolated bilinearly in each map, and in time between maps rotated with the Sun (so-called IONEX method 3).

Geometries:
- `LineOfSight(site, source)`: a geodetic WGS84 `site` `(; lat, lon, h)`, looking toward a source represented as `SkyCoords`.
- `Ray(origin, direction)`: an ECEF origin (Unitful lengths) inside the shell and an Earth-fixed unit direction; `Ray(site, direction)` from a geodetic site. `Ray(los::LineOfSight, t)` is the ray of `los` at time `t`.

Slant TEC is VTEC times the mapping factor. Available mappings:
- `SLM()`: the IONEX single layer, exact `1/cos z′` at the pierce point on the map's shell.
- `MSLM()`: CODE's modified single-layer mapping (506.7 km, α = 0.9782), the one CODE produces its maps with.
- `Petrov2023()`: Petrov (2023, AJ 165, 183) mapping fitted to VLBI astrometry, its TEC scale factor k = 0.85 included.

`vtec(m, p, t)` is the maps' own VTEC at a point `p = (; lat, lon)` on the shell (geocentric); the mapping plays no role.

Convenience computation functions:
- RM is `K·TEC·B∥` with the IGRF field at the pierce point,
- group delay is `K·TEC / (c ν²)` at the cyclic frequency `ν`.

## Uncertainties

For files with RMS maps, every result is `value ± σ` (from `Uncertain.jl`); without RMS maps, a plain Unitful quantity.
`Uncertain.value` returns the nominal value of either.

The RMS is interpolated with the same weights as the TEC and scaled by the same factors.
The reported uncertainty covers the map RMS only: not the thin-shell and mapping approximations or the field model.
