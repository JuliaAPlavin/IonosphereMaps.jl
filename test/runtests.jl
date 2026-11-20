using TestItems
using TestItemRunner
@run_package_tests

@testmodule Track begin
    using Dates, SkyCoords, Unitful, IonosphericRM
    using IonosphericRM: ecef, source_direction
    # MeerKAT array centre and 3C138 (J2000); MeerKAT 2025-07-25 scan mid-times (UTC)
    const DATA = joinpath(@__DIR__, "data")
    const SITE = (; lat=-30.71292u"°", lon=21.44380u"°", h=1050.0u"m")
    const SRC = ICRSCoords(deg2rad(15 * (5 + 21/60 + 9.886/3600)), deg2rad(16 + 38/60 + 22.05/3600))
    const SCANS = DateTime(2025, 7, 25) .+ Minute.([366, 386, 403, 424, 441, 462, 479, 500, 517, 537, 555, 574, 593, 612])
    const LOS = LineOfSight(SITE, SRC)
    # the same line of sight as an Earth-fixed ray, at time t
    ray(t) = Ray(ecef(SITE), source_direction(SRC, t))
    const GEOMETRIES = [Returns(LOS), ray]
    const TECU = 1e16u"m^-2"
end

@testitem "constants and geometry" setup=[Track] begin
    using IonosphericRM: K_RM, ecef, geodetic_up, pierce, source_direction
    using LinearAlgebra, Dates, SkyCoords, Unitful
    @test K_RM * 1u"m^-2" * 1u"T" ≈ 2.6312e-13u"rad/m^2" rtol=1e-4
    o = ecef(Track.SITE)
    @test pierce(o, normalize(o), 6821u"km") ≈ 6821u"km" * normalize(o)
    # 3C138 elevation, astropy (via spinifex) reference: 37.08263° at 06:06, 30.659° (ALBUS) at 10:12
    el(t) = rad2deg(asin(dot(source_direction(Track.SRC, t), geodetic_up(ecef(Track.SITE)))))
    @test el(Track.SCANS[1]) ≈ 37.0826 atol=0.01
    @test el(Track.SCANS[end]) ≈ 30.659 atol=0.02
    # Earth-fixed direction, SOFA c2t06a (IAU 2006/2000A, no polar motion, UT1 = UTC); also near the celestial poles
    angle(a, b) = rad2deg(acos(clamp(dot(a, b), -1, 1))) * 3600
    @testset for (src, t, ref) in [
            (ICRSCoords(deg2rad(30.0), deg2rad(89.9)), DateTime(2025, 7, 25, 6, 6), [-0.0003361819, 0.0012448534, 0.9999991687]),
            (Track.SRC, DateTime(2025, 7, 25, 6, 6), [0.6660455275, 0.6885773449, 0.2867831853]),
            (ICRSCoords(deg2rad(200.0), deg2rad(-89.5)), DateTime(2031, 1, 3, 17, 45), [-0.0055652884, -0.0021720281, -0.9999821548])]
        @test angle(source_direction(src, t), normalize(ref)) < 0.5   # arcsec
    end
end

@testitem "global map: verified reference values" setup=[Track] begin
    using Dates, SkyCoords, Uncertain, Unitful
    path = joinpath(Track.DATA, "COD0OPSFIN_20252060000_01D_01H_GIM_crop.INX.gz")
    # independently verified calculator (vs spinifex 2.0 and AIPS TECOR), CODE final GIM 2025-206
    slm = [-1.0456, -1.1890, -1.3138, -1.4600, -1.5775, -1.7246, -1.8365, -2.0016, -2.1344, -2.2655, -2.3375, -2.4113, -2.4826, -2.5611]u"rad/m^2"
    mslm = [-1.0184, -1.1611, -1.2851, -1.4303, -1.5465, -1.6913, -1.8009, -1.9613, -2.0891, -2.2135, -2.2788, -2.3440, -2.4047, -2.4701]u"rad/m^2"
    @testset for geometry in Track.GEOMETRIES
        @test U.value.(ionospheric_rm.(Ref(GlobalMap(path; mapping=SLM())), geometry.(Track.SCANS), Track.SCANS)) ≈ slm atol=3e-4u"rad/m^2"   # norm over the 14 scans; references to 1e-4, sky direction 7″ off ours
        @test U.value.(ionospheric_rm.(Ref(GlobalMap(path; mapping=MSLM())), geometry.(Track.SCANS), Track.SCANS)) ≈ mslm atol=3e-4u"rad/m^2"
    end
    # circumpolar source, time outside the map epochs
    @test_throws DomainError ionospheric_rm(GlobalMap(path; mapping=MSLM()), LineOfSight(Track.SITE, ICRSCoords(0.0, -1.2)), DateTime(2025, 7, 25, 20))
    # any sky frame: 3C138 in galactic coordinates
    m = GlobalMap(path; mapping=MSLM())
    gal, icrs = ionospheric_rm.(Ref(m), LineOfSight.(Ref(Track.SITE), [convert(GalCoords, Track.SRC), Track.SRC]), Track.SCANS[1])
    @test U.value(gal) ≈ U.value(icrs) rtol=1e-12
    @test U.uncertainty(gal) ≈ U.uncertainty(icrs) rtol=1e-12
    # site latitude and longitude as plain radians or as Unitful angles; height needs its unit
    plain = (; lat=deg2rad(-30.71292), lon=deg2rad(21.44380), h=1050.0u"m")
    @test U.value(ionospheric_rm(m, LineOfSight(plain, Track.SRC), Track.SCANS[1])) ≈ U.value(icrs) rtol=1e-12
    @test_throws Unitful.DimensionError ionospheric_rm(m, LineOfSight((; Track.SITE..., h=1050.0), Track.SRC), Track.SCANS[1])
end

@testitem "horizon and non-IONEX input" setup=[Track] begin
    using Dates, SkyCoords, Uncertain, Unitful
    m = GlobalMap(joinpath(Track.DATA, "COD0OPSFIN_20252060000_01D_01H_GIM_crop.INX.gz"); mapping=MSLM())
    # below the horizon: 3C138 is up throughout these maps (03:00–13:00), a source at +80° never rises at MeerKAT
    north = LineOfSight(Track.SITE, ICRSCoords(0.0, deg2rad(80)))
    r = @inferred ionospheric_rm(m, north, DateTime(2025, 7, 25, 8))
    s = @inferred ionospheric_tec(m, north, DateTime(2025, 7, 25, 8))
    @test isnan(U.value(r)) && isnan(U.uncertainty(r))
    @test isnan(U.value(s)) && isnan(U.uncertainty(s))
    m0 = GlobalMap(joinpath(Track.DATA, "codg0600.21i_crop.gz"); mapping=MSLM())
    @test isnan(@inferred ionospheric_rm(m0, Track.LOS, DateTime(2021, 3, 1, 3)))
    @test isnan(@inferred ionospheric_tec(m0, Track.LOS, DateTime(2021, 3, 1, 3)))
    # a time outside the maps is an error also below the horizon: 3C138 at 01:00
    @testset for geometry in Track.GEOMETRIES
        t = DateTime(2025, 7, 25, 1)
        @test_throws DomainError ionospheric_rm(m, geometry(t), t)
        @test_throws DomainError ionospheric_tec(m, geometry(t), t)
        @test_throws DomainError ionospheric_group_delay(m, geometry(t), t, 1u"GHz")
    end
    @test_throws ArgumentError("$(@__FILE__) has no MAP DIMENSION record from line 1 on") GlobalMap(@__FILE__; mapping=MSLM())
    @test_throws UndefKeywordError GlobalMap(joinpath(Track.DATA, "COD0OPSFIN_20252060000_01D_01H_GIM_crop.INX.gz"))
end

@testitem "RMS maps: VTEC and RM uncertainties" setup=[Track] begin
    using Dates, Uncertain, Unitful
    using .Track: TECU
    points = [(-23.98, 24.26), (-39.17, 26.65), (-30.85, 21.34), (-25.0, 40.0)]
    # spinifex 2.0 interpolate_ionex: (TEC, RMS) in TECU, Earth rotation on, UQRG RMS correction off
    cod = [
        DateTime(2025, 7, 25, 6, 0)     => [(14.669162, 1.314800), (11.080644, 2.633000), (9.258624, 1.382312), (18.900000, 1.600000)],
        DateTime(2025, 7, 25, 6, 7, 30) => [(15.348065, 1.385412), (11.882880, 2.718849), (10.110481, 1.489642), (19.693750, 1.600000)],
        DateTime(2025, 7, 25, 10, 20)   => [(30.792941, 1.446933), (32.429112, 2.681652), (28.865619, 1.552504), (30.133333, 1.500000)]]
    upc = [
        DateTime(2025, 7, 25, 6, 0)     => [(15.877715, 3.206038), (11.038096, 3.311288), (12.772272, 3.248312), (20.500000, 3.200000)],
        DateTime(2025, 7, 25, 6, 7, 30) => [(16.496600, 3.214119), (11.221681, 3.327341), (13.349103, 3.217131), (21.718750, 3.231250)],
        DateTime(2025, 7, 25, 10, 20)   => [(32.118028, 3.268679), (35.637708, 4.017779), (31.646952, 3.220445), (31.883333, 3.250000)],
        DateTime(2025, 7, 25, 23, 55)   => [(7.560292, 3.340718), (4.810465, 3.828603), (7.626979, 3.306528), (7.200000, 3.233333)]]
    @testset for (file, ref) in ["COD0OPSFIN_20252060000_01D_01H_GIM_crop.INX.gz" => cod, "UPC0OPSRAP_20252060000_01D_15M_GIM_crop.INX.Z" => upc]
        m = GlobalMap(joinpath(Track.DATA, file); mapping=SLM())
        v = [vtec(m, lat * u"°", lon * u"°", t) for (t, _) in ref for (lat, lon) in points]
        expected = reduce(vcat, last.(ref))
        @test U.value.(v) ≈ first.(expected) .* TECU atol=3e-6TECU   # norm over values given to 1e-6
        @test U.uncertainty.(v) ≈ last.(expected) .* TECU atol=3e-6TECU
    end
    # no RMS maps: a plain value; at a node at a map epoch it is the file's, 33 × 0.1 TECU
    @test (@inferred vtec(GlobalMap(joinpath(Track.DATA, "codg0600.21i_crop.gz"); mapping=MSLM()), -30.0u"°", 20.0u"°", DateTime(2021, 3, 1, 3)))::Unitful.Quantity ≈
          3.3e16u"m^-2" rtol=1e-14

    # spinifex 2.0 get_rm_from_skycoord, single layer at 450 km (= SLM): RM and the TEC part of its error, |K·RMS·M·B∥|,
    # rescaled from its rounded constant 2.62e-6 to K_RM
    m = GlobalMap(joinpath(Track.DATA, "COD0OPSFIN_20252060000_01D_01H_GIM_crop.INX.gz"); mapping=SLM())
    @testset for geometry in Track.GEOMETRIES
        rm = [ionospheric_rm(m, geometry(t), t) for t in DateTime(2025, 7, 25) .+ Minute.([366, 480, 612])]
        @test U.value.(rm) ≈ [-1.045650, -1.843167, -2.561181]u"rad/m^2" atol=3e-4u"rad/m^2"
        @test U.uncertainty.(rm) ≈ [0.105383, 0.106252, 0.136799]u"rad/m^2" atol=3e-5u"rad/m^2"
        t = DateTime(2025, 7, 25, 8)
        @test (@inferred ionospheric_rm(m, geometry(t), t)) isa Uncertain.Value
    end

    # VTEC at a shell point (from the spinifex values above), angles Unitful or plain radians; the mapping plays no role
    t = DateTime(2025, 7, 25, 6)
    v = @inferred vtec(m, -23.98u"°", 24.26u"°", t)
    @test U.value(v) ≈ 14.669162e16u"m^-2" rtol=1e-6
    @test U.uncertainty(v) ≈ 1.314800e16u"m^-2" rtol=1e-6
    @test vtec(m, deg2rad(-23.98), deg2rad(24.26), t) == v
    @test vtec(GlobalMap(joinpath(Track.DATA, "COD0OPSFIN_20252060000_01D_01H_GIM_crop.INX.gz"); mapping=Petrov2023()), -23.98u"°", 24.26u"°", t) == v
    # slant TEC along the line of sight, in m⁻²: spinifex VTEC and RMS × its airmass at the same three times;
    # its line of sight includes aberration, 15″ here, ours does not
    @testset for geometry in Track.GEOMETRIES
        s = [ionospheric_tec(m, geometry(t), t) for t in DateTime(2025, 7, 25) .+ Minute.([366, 480, 612])]
        @test U.value.(s) ≈ [20.35420, 31.72016, 49.29405]u"m^-2" .* 1e16 rtol=1e-4
        @test U.uncertainty.(s) ≈ [2.05135, 1.82856, 2.63291]u"m^-2" .* 1e16 rtol=1e-4
        t = DateTime(2025, 7, 25, 8)
        @test (@inferred ionospheric_tec(m, geometry(t), t)) isa Uncertain.Value
    end
end

@testitem "Ray: vertical at a grid node, geodetic horizon, construction" setup=[Track] begin
    using Dates, StaticArrays, Uncertain, Unitful
    # from the WGS84 equator at 90°E straight up: pierces the shell at the grid node 0°N 90°E;
    # its VTEC and RMS in the 06:00 maps of the file are 719 and 25 (× 0.1 TECU)
    m = GlobalMap(joinpath(Track.DATA, "COD0OPSFIN_20252060000_01D_01H_GIM_seam.INX.gz"); mapping=SLM())
    v = @inferred ionospheric_tec(m, Ray(SVector(0.0, 6378137.0, 0.0)u"m", SVector(0.0, 1.0, 0.0)), DateTime(2025, 7, 25, 6))
    @test U.value(v) == 71.9e16u"m^-2"
    @test U.uncertainty(v) == 2.5e16u"m^-2"

    # at MeerKAT, toward the north at -0.05° geodetic elevation: below the geodetic horizon,
    # above the geocentric one (+0.12°)
    lat, lon = -30.71292u"°", 21.44380u"°"
    up = SVector(cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat))
    north = SVector(-sin(lat) * cos(lon), -sin(lat) * sin(lon), cos(lat))
    ray = Ray(SVector(5109235.683, 2006794.555, -3239107.175)u"m", cosd(-0.05) * north + sind(-0.05) * up)
    t = DateTime(2025, 7, 25, 8)
    @testset for mapping in [SLM(), MSLM()]
        m = GlobalMap(joinpath(Track.DATA, "COD0OPSFIN_20252060000_01D_01H_GIM_crop.INX.gz"); mapping)
        @test isnan(U.value(@inferred ionospheric_tec(m, ray, t)))
        @test isnan(U.uncertainty(ionospheric_tec(m, ray, t)))
        @test isnan(U.value(@inferred ionospheric_rm(m, ray, t)))
        @test isnan(U.uncertainty(ionospheric_rm(m, ray, t)))
        m0 = GlobalMap(joinpath(Track.DATA, "codg0600.21i_crop.gz"); mapping)
        @test isnan(@inferred ionospheric_tec(m0, ray, DateTime(2021, 3, 1, 3)))
        @test isnan(@inferred ionospheric_rm(m0, ray, DateTime(2021, 3, 1, 3)))
    end

    # an origin at or above the maps' shell (6371 + 450 km) is an error, also below the horizon (pointing down)
    m = GlobalMap(joinpath(Track.DATA, "COD0OPSFIN_20252060000_01D_01H_GIM_seam.INX.gz"); mapping=SLM())
    @testset for (r, d) in [(6821.0, 1.0), (7000.0, 1.0), (7000.0, -1.0)]
        ray = Ray(SVector(0.0, r, 0.0)u"km", SVector(0.0, d, 0.0))
        err = DomainError(r * u"km", "ray origin distance from the geocentre, at or above the ionospheric shell of radius 6821.0 km")
        @test_throws err ionospheric_tec(m, ray, DateTime(2025, 7, 25, 6))
        @test_throws err ionospheric_rm(m, ray, DateTime(2025, 7, 25, 6))
    end

    # an ECEF origin with length units, a unit direction
    @test_throws MethodError Ray(SVector(0.0, 6378137.0, 0.0), SVector(0.0, 1.0, 0.0))
    @test_throws AssertionError Ray(SVector(0.0, 6378137.0, 0.0)u"m", SVector(0.0, 1.0, 1.0))
end

@testitem "Petrov 2023 mapping" begin
    using Dates, AxisKeys, StaticArrays, Unitful
    # a map of 1 TECU everywhere: slant TEC in TECU is the mapping factor
    t = DateTime(2025, 1, 1)
    uniform(radius, height) = GlobalMap((; tec=KeyedArray(ones(71, 73, 1); lat=(87.5:-2.5:-87.5)u"°", lon=(-180.0:5.0:180.0)u"°", epoch=[t]),
                                           rms=nothing, radius, height), Petrov2023())
    # from 45°N 0°E on the reference sphere, northward at geocentric elevation e (geodetic elevation 0.19° higher)
    ray(R, e) = Ray(R * SVector(cosd(45), 0, sind(45)), sind(e) * SVector(cosd(45), 0, sind(45)) + cosd(e) * SVector(-sind(45), 0, cosd(45)))
    factor(m, R, e) = (@inferred ionospheric_tec(m, ray(R, e), t)) / 1e16u"m^-2"
    # eq. 12 of the paper in 256-bit BigFloat, outside the package
    m = uniform(6371.0u"km", 450.0u"km")
    M = factor.(Ref(m), 6371.0u"km", [90, 45, 30, 20, 15, 10, 5])
    @test M ≈ [0.850427787832, 1.13966199380, 1.44085829434, 1.74183394892, 1.91520009405, 2.08193646078, 2.20830802078] rtol=1e-11
    # Fringy's independent implementation, M/k to 3 decimals
    @test round.(M ./ 0.85; digits=3) == [1.001, 1.341, 1.695, 2.049, 2.253, 2.449, 2.598]
    # R and H are the map's own
    @test factor(uniform(6378.137u"km", 350.0u"km"), 6378.137u"km", 10) ≈ 2.25724822057 rtol=1e-11
end

@testitem "group delay" setup=[Track] begin
    using Dates, AxisKeys, SkyCoords, StaticArrays, Uncertain, Unitful
    t = DateTime(2025, 1, 1)
    tecu = GlobalMap((; tec=KeyedArray(ones(3, 3, 1); lat=(5.0:-5.0:-5.0)u"°", lon=(-5.0:5.0:5.0)u"°", epoch=[t]), rms=nothing,
                        radius=6371.0u"km", height=450.0u"km"), SLM())
    # 1 TECU along a vertical ray from the equator at 0°E: e²/(8π²ε₀mₑ) · 10¹⁶ m⁻² / (c ν²) in 256-bit BigFloat
    # from CODATA 2018 values (Unitful's ε₀ = 1/(μ₀c²) with μ₀ = 4π·10⁻⁷ H/m is 5·10⁻¹⁰ off)
    delay(ν) = @inferred ionospheric_group_delay(tecu, Ray(SVector(6378137.0, 0.0, 0.0)u"m", SVector(1.0, 0.0, 0.0)), t, ν)
    @test delay(1u"GHz") ≈ 1.344536593446e-9u"s" rtol=1e-9
    @test delay(15.368u"GHz") ≈ 5.692957254713e-12u"s" rtol=1e-9
    @test delay(1u"GHz") ≈ 4 * delay(2e9u"Hz") rtol=1e-14

    # end to end at 1.4 GHz: spinifex slant TEC and its RMS at 06:06 (in "RMS maps" above), 20.35420 ± 2.05135 TECU,
    # times the 1 TECU delay at 1.4 GHz from the same BigFloat evaluation
    m = GlobalMap(joinpath(Track.DATA, "COD0OPSFIN_20252060000_01D_01H_GIM_crop.INX.gz"); mapping=SLM())
    @testset for geometry in Track.GEOMETRIES
        t = Track.SCANS[1]
        d = @inferred ionospheric_group_delay(m, geometry(t), t, 1.4u"GHz")
        @test U.value(d) ≈ 13.962738128e-9u"s" rtol=1e-4
        @test U.uncertainty(d) ≈ 1.4072016025e-9u"s" rtol=1e-4
    end

    # below the horizon (a source at +80° never rises at MeerKAT), with and without RMS maps
    d = @inferred ionospheric_group_delay(m, LineOfSight(Track.SITE, ICRSCoords(0.0, deg2rad(80))), DateTime(2025, 7, 25, 8), 1u"GHz")
    @test isnan(U.value(d))
    @test isnan(U.uncertainty(d))
    m0 = GlobalMap(joinpath(Track.DATA, "codg0600.21i_crop.gz"); mapping=MSLM())
    @test isnan(@inferred ionospheric_group_delay(m0, Track.LOS, DateTime(2021, 3, 1, 3), 1u"GHz"))
    @test_throws MethodError ionospheric_group_delay(m, Track.LOS, DateTime(2025, 7, 25, 8), 1e9)
end

@testitem "grid interpolation: longitude wrap, seam, grid edges, map epochs" setup=[Track] begin
    using Dates, AxisKeys
    using Unitful: °
    using .Track: TECU
    t0 = DateTime(2025, 1, 1)
    f(lat, lon) = 10 + 3cosd(lon) + 2sind(2lon) + lat / 10
    lats = 10.0:-2.5:-10.0
    # two maps an hour apart, the second rotated 15° west as if fixed to the Sun; VTEC lookup does not use the mapping
    ionex(lons; fixlon=identity) = GlobalMap((; rms=nothing,
        tec=KeyedArray(stack([fixlon(f.(lats, lons')), fixlon(f.(lats, lons' .+ 15))]); lat=lats * °, lon=lons * °, epoch=[t0, t0 + Hour(1)])), SLM())

    # a full circle with or without the repeated ±180 column, starting anywhere: the same field
    grids = ionex.([-180.0:5.0:180.0, -180.0:5.0:175.0, 0.0:5.0:355.0, 0.0:5.0:360.0])
    pts = [(t, lat, lon) for t in t0 .+ Minute.([0, 20, 60]) for lat in (7.0, -10.0) for lon in (-180, -179, -1, 0, 1, 176, 177.5, 179, 180, 181, 359)]
    v = [[vtec(m, lat * °, lon * °, t) for (t, lat, lon) in pts] for m in grids]
    @test all(≈(v[1]; rtol=1e-12), v)
    @test vtec(grids[1], 5.0°, 177.5°, t0) ≈ (f(5, 175) + f(5, 180)) / 2 * TECU
    @test vtec(grids[2], 5.0°, 177.5°, t0) ≈ (f(5, 175) + f(5, -180)) / 2 * TECU
    @test vtec(grids[1], 5.0°, 2.5°, t0 + Hour(1)) ≈ (f(5, 15) + f(5, 20)) / 2 * TECU

    # ±180 columns that differ: each side of the seam interpolates from its own column
    seam = ionex(-180.0:5.0:180.0; fixlon=M -> [M[:, 1:end-1] M[:, end] .+ 10])
    @test vtec(seam, 5.0°, 179.0°, t0) ≈ (0.2f(5, 175) + 0.8(f(5, 180) + 10)) * TECU
    @test vtec(seam, 5.0°, -179.0°, t0) ≈ (0.8f(5, -180) + 0.2f(5, -175)) * TECU

    # regional grid: grid edges included, beyond them an error, not an extrapolation
    reg = ionex(-10.0:5.0:55.0)
    @test vtec(reg, -10.0°, 55.0°, t0) ≈ f(-10, 55) * TECU
    @test vtec(reg, 10.0°, -10.0°, t0) ≈ f(10, -10) * TECU
    @test_throws DomainError vtec(reg, 0.0°, 56.0°, t0)
    @test_throws DomainError vtec(reg, 0.0°, -11.0°, t0)
    @test_throws DomainError vtec(reg, 0.0°, 50.0°, t0 + Minute(30))   # rotated out of the grid in the second map
    @test_throws DomainError vtec(grids[1], 10.5°, 0.0°, t0)
    @test_throws DomainError vtec(grids[1], -12.0°, 0.0°, t0)

    # map epochs: both ends included, a file with a single map
    @test vtec(grids[1], 0.0°, 20.0°, t0 + Hour(1)) ≈ f(0, 35) * TECU
    single = GlobalMap((; rms=nothing, tec=grids[1].ionex.tec[:, :, 1:1]), SLM())
    @test vtec(single, 0.0°, 20.0°, t0) ≈ f(0, 20) * TECU
    @test_throws DomainError vtec(single, 0.0°, 20.0°, t0 + Minute(1))
    @test_throws DomainError vtec(single, 0.0°, 20.0°, t0 - Minute(1))
    @test_throws DomainError vtec(grids[1], 0.0°, 20.0°, t0 - Minute(1))
    @test_throws DomainError vtec(grids[1], 0.0°, 20.0°, t0 + Minute(61))
end

@testitem "grid interpolation: polar caps" setup=[Track] begin
    using Dates, AxisKeys, StaticArrays, Uncertain
    using Unitful: @u_str, °
    using IonosphericRM: ecef
    using .Track: TECU
    t = DateTime(2025, 1, 1)
    # rows at 60, 30, 0, -30, -60°, columns at -180, -90, 0, 90°; ring means: TEC 3 (north), 25 (south), RMS 2 (north)
    tec = [1 2 3 6; 7 7 7 7; 8 8 8 8; 9 9 9 9; 10 20 30 40]
    rms = [1 1 1 5; 1 1 1 1; 1 1 1 1; 1 1 1 1; 1 1 1 1]
    keyed(M, lats, lons) = KeyedArray(reshape(M, size(M)..., 1); lat=lats * °, lon=lons * °, epoch=[t])
    gmap(T, R, lats, lons) = GlobalMap((; tec=keyed(T, lats, lons), rms=keyed(R, lats, lons)), SLM())
    grids = [gmap(tec, rms, 60.0:-30.0:-60.0, -180.0:90.0:90.0),
             gmap(reverse(tec; dims=1), reverse(rms; dims=1), -60.0:30.0:60.0, -180.0:90.0:90.0),   # ascending latitudes
             gmap([tec tec[:, 1]], [rms rms[:, 1]], 60.0:-30.0:-60.0, -180.0:90.0:180.0),          # repeated ±180 column
             gmap(circshift(tec, (0, -2)), circshift(rms, (0, -2)), 60.0:-30.0:-60.0, 0.0:90.0:270.0)]
    # (lat, lon) => (TEC, RMS): poles at any longitude; the 60° row from the cap side; inside the caps
    # T = w·ring(lon) + (1 − w)·P, w = (90° − |lat|) / 30°
    cases = [(90.0, 0.0) => (3, 2), (90.0, -135.0) => (3, 2), (90.0, 171.0) => (3, 2), (-90.0, 0.0) => (25, 1), (-90.0, 100.0) => (25, 1),
             (60.0, -135.0) => (1.5, 1), (60.0 + 1e-12, -135.0) => (1.5, 1), (-60.0 - 1e-12, 45.0) => (35, 1),
             (75.0, -135.0) => (2.25, 1.5), (80.0, 45.0) => (3.5, 7/3), (-70.0, 112.5) => (30, 1)]
    @testset for m in grids
        v = [vtec(m, lat * °, lon * °, t) for ((lat, lon), _) in cases]
        @test U.value.(v) ≈ [T for (_, (T, _)) in cases] .* TECU rtol=1e-12
        @test U.uncertainty.(v) ≈ [σ for (_, (_, σ)) in cases] .* TECU rtol=1e-12
        @test_throws DomainError vtec(m, 90.5°, 0.0°, t)
    end
    # beyond the grid otherwise: not all longitudes, or the outermost row not one step from the pole
    @test_throws DomainError vtec(gmap(tec[:, 1:3], rms[:, 1:3], 60.0:-30.0:-60.0, -180.0:90.0:0.0), 75.0°, -135.0°, t)
    @test_throws DomainError vtec(gmap(tec, rms, 50.0:-25.0:-50.0, -180.0:90.0:90.0), 60.0°, 0.0°, t)
    @test_throws DomainError vtec(gmap(tec, rms, 50.0:-25.0:-50.0, -180.0:90.0:90.0), -60.0°, 0.0°, t)

    # CODE final GIM 2025-243, its rows 87.5–80°N and maps 06:00–07:00: from Ny-Ålesund (78.93°N, 11.87°E, 87 m) north
    # at 12° geodetic elevation, the ray crosses the pole and pierces the shell in the cap at 88.846°N, 168.13°W.
    # Reference: numpy, own WGS84 site and pierce point, ring and pole mean read from the file text, 06:20 between the
    # two rotated maps, times 1/cos z′ at the pierce point
    m = GlobalMap(joinpath(Track.DATA, "COD0OPSFIN_20252430000_01D_01H_GIM_polar.INX.gz"); mapping=SLM())
    lat, lon = 78.93°, 11.87°
    up = SVector(cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat))
    north = SVector(-sin(lat) * cos(lon), -sin(lat) * sin(lon), cos(lat))
    s = @inferred ionospheric_tec(m, Ray(ecef((; lat, lon, h=87.0u"m")), cosd(12) * north + sind(12) * up), DateTime(2025, 8, 31, 6, 20))
    @test U.value(s) ≈ 40.244515712906505TECU rtol=1e-12
    @test U.uncertainty(s) ≈ 7.990758323319431TECU rtol=1e-12
end

@testitem "grid interpolation: real maps across the ±180 seam" setup=[Track] begin
    using Dates, Uncertain, Unitful
    using .Track: TECU
    points = [(0.0, 179.0), (0.0, -179.0), (10.0, 357.5)]
    # spinifex 2.0 interpolate_ionex (TEC, RMS), with its lookup of the two nodes around a point near ±180 fixed;
    # UPC's -180 and +180 columns differ, each side of the seam uses its own
    upc = [
        DateTime(2025, 7, 25, 6, 0)     => [(51.160000, 3.380000), (51.540000, 3.380000), (7.600000, 3.400000)],
        DateTime(2025, 7, 25, 6, 7, 30) => [(51.113750, 3.355000), (48.973750, 3.395000), (9.050000, 3.487500)],
        DateTime(2025, 7, 25, 10, 20)   => [(34.843333, 3.300000), (35.923333, 3.300000), (47.233333, 3.300000)],
        DateTime(2025, 7, 25, 23, 55)   => [(53.923333, 3.550000), (54.376667, 3.616667), (15.600000, 3.716667)]]
    cod = [
        DateTime(2025, 7, 25, 6, 0)     => [(61.040000, 2.020000), (60.660000, 2.180000), (11.400000, 2.000000)],
        DateTime(2025, 7, 25, 6, 7, 30) => [(60.136562, 2.098750), (59.391562, 2.238750), (12.407813, 2.081250)],
        DateTime(2025, 7, 25, 23, 55)   => [(60.788750, 1.861250), (60.865417, 1.984583), (11.675000, 2.097917)]]
    @testset for (file, ref) in ["UPC0OPSRAP_20252060000_01D_15M_GIM_seam.INX.Z" => upc, "COD0OPSFIN_20252060000_01D_01H_GIM_seam.INX.gz" => cod]
        m = GlobalMap(joinpath(Track.DATA, file); mapping=SLM())
        v = [vtec(m, lat * u"°", lon * u"°", t) for (t, _) in ref for (lat, lon) in points]
        expected = reduce(vcat, last.(ref))
        @test U.value.(v) ≈ first.(expected) .* TECU atol=3e-6TECU
        @test U.uncertainty.(v) ≈ last.(expected) .* TECU atol=3e-6TECU
    end
end

@testitem "gzip, unix compress and plain files read alike" setup=[Track] begin
    using IonosphericRM: read_ionex
    gz, Z = joinpath.(Track.DATA, ["codg0590.21i_crop.gz", "codg0590.21i_crop.Z"])
    plain = tempname()
    write(plain, read(`gzip -dc $gz`))
    @test isequal(read_ionex(Z), read_ionex(gz))
    @test isequal(read_ionex(plain), read_ionex(gz))
end

@testitem "multi-day maps: the later file's midnight map, continuous across midnight" setup=[Track] begin
    using Dates, SkyCoords, AxisKeys, Unitful
    using .Track: TECU
    # CODE final GIMs for 2021-02-28 (21–24 h) and 2021-03-01 (0–6 h), cropped; both hold a 03-01 00:00 map
    day59, day60 = joinpath.(Track.DATA, ["codg0590.21i_crop.gz", "codg0600.21i_crop.gz"])
    m59 = GlobalMap(day59; mapping=MSLM()); m60 = GlobalMap(day60; mapping=MSLM())
    m = GlobalMap([day59, day60]; mapping=MSLM())
    @test axiskeys(m.ionex.tec, :epoch) == DateTime(2021, 2, 28, 21):Hour(1):DateTime(2021, 3, 1, 6)
    @test m.ionex.tec == cat(m59.ionex.tec[:, :, 1:3], m60.ionex.tec; dims=3)
    @test m59.ionex.tec[:, :, end] != m60.ionex.tec[:, :, 1]   # two solutions of the same epoch
    @test GlobalMap((day59, day60); mapping=MSLM()).ionex == m.ionex
    @test GlobalMap([day59]; mapping=MSLM()).ionex == m59.ionex

    # VTEC (TECU) from spinifex 2.0 with remove_midnight_jumps, at (lat, lon) points over MeerKAT
    points = [(-23.98, 24.26), (-39.17, 26.65), (-30.85, 21.34), (-25.0, 40.0)]
    spinifex = [
        DateTime(2021, 2, 28, 23, 0)  => [4.936323, 4.000000, 3.931912, 4.000000],
        DateTime(2021, 2, 28, 23, 20) => [4.951028, 4.118504, 3.767187, 4.233333],
        DateTime(2021, 2, 28, 23, 40) => [5.071574, 4.196771, 3.658741, 4.166667],
        DateTime(2021, 2, 28, 23, 59) => [5.258914, 4.078122, 3.570119, 3.918083],
        DateTime(2021, 3, 1, 0, 0)    => [5.268362, 4.067600, 3.558712, 3.900000],
        DateTime(2021, 3, 1, 0, 20)   => [5.153949, 3.982052, 3.629691, 4.100000],
        DateTime(2021, 3, 1, 0, 40)   => [4.949695, 3.911511, 3.637171, 4.033333]]
    @test [vtec(m, lat * u"°", lon * u"°", t) for (t, _) in spinifex for (lat, lon) in points] ≈ reduce(vcat, last.(spinifex)) .* TECU atol=3e-6TECU   # norm over 28 values given to 1e-6

    # 3C286: a single daily file jumps at midnight, the merged map does not
    src = ICRSCoords(deg2rad(15 * (13 + 31/60 + 8.288/3600)), deg2rad(30 + 30/60 + 32.96/3600))
    rm(model, t) = ionospheric_rm(model, LineOfSight(Track.SITE, src), t)
    midnight = DateTime(2021, 3, 1)
    @test [rm(m59, midnight), rm(m60, midnight)] ≈ [-0.4833, -0.4795]u"rad/m^2" atol=1e-4u"rad/m^2"
    @test rm(m, midnight - Millisecond(1)) ≈ rm(m, midnight) atol=1e-7u"rad/m^2"
    @test rm(m, midnight - Hour(1)) == rm(m59, midnight - Hour(1))
    @test rm.(Ref(m), midnight .+ Minute.(0:30:360)) == rm.(Ref(m60), midnight .+ Minute.(0:30:360))
    @test_throws DomainError vtec(m, -30.0u"°", 21.0u"°", DateTime(2021, 3, 1, 6, 30))
    @test_throws DomainError vtec(m, -30.0u"°", 21.0u"°", DateTime(2021, 2, 28, 20, 30))
end

@testitem "multi-day maps: EXPONENT records apply per file and per map" setup=[Track] begin
    using AxisKeys
    # the 03-01 file rewritten at exponent -2, stated in a record inside each map (header still -1)
    day59, day60, day60e = joinpath.(Track.DATA, ["codg0590.21i_crop.gz", "codg0600.21i_crop.gz", "codg0600.21i_crop_exp2.gz"])
    a = GlobalMap([day59, day60]; mapping=MSLM()).ionex; b = GlobalMap([day59, day60e]; mapping=MSLM()).ionex
    @test axiskeys(b.tec) == axiskeys(a.tec)
    @test b.tec ≈ a.tec rtol=1e-14
end

@testitem "IONEX records: fixed-width fields, record order, completeness, missing values" setup=[Track] begin
    using Dates, AxisKeys, Unitful
    using IonosphericRM: read_ionex, headerlabel
    day59, day60, day60e2, day60e3 = joinpath.(Track.DATA, ["codg0590.21i_crop.gz", "codg0600.21i_crop.gz", "codg0600.21i_crop_exp2.gz", "codg0600.21i_crop_exp3.gz"])
    ref = read_ionex(day60)
    # exponent -3: values from 10000 fill all 5 columns and touch their neighbours
    @test read_ionex(day60e3).tec ≈ ref.tec rtol=1e-14

    function rewritten(f, path)
        out = tempname()
        write(out, join(f(readlines(`gzip -dc $path`)), "\n"))
        out
    end
    starts(lines) = findall(l -> headerlabel(l) == "START OF TEC MAP", lines)
    # EXPONENT record before EPOCH OF CURRENT MAP
    swapped = rewritten(day60e2) do lines
        foldl((ls, i) -> ls[[1:i; i+2; i+1; i+3:end]], starts(lines); init=lines)
    end
    @test isequal(read_ionex(swapped), read_ionex(day60e2))
    # cut between maps or inside a map (the last map starts at line 837); a header record missing
    @test_throws ArgumentError read_ionex(rewritten(ls -> ls[1:starts(ls)[end]-1], day60))
    cut = rewritten(ls -> ls[1:starts(ls)[end]+5], day60)
    @test_throws ArgumentError("$cut has no END OF TEC MAP record from line 837 on") read_ionex(cut)
    nolon = rewritten(ls -> filter(l -> headerlabel(l) != "LON1 / LON2 / DLON", ls), day60)
    @test_throws ArgumentError("$nolon has no LON1 / LON2 / DLON record from line 1 on") read_ionex(nolon)
    # 3-D maps
    dim3 = rewritten(ls -> map(l -> headerlabel(l) == "MAP DIMENSION" ? "     3" * l[7:end] : l, ls), day60)
    @test_throws ArgumentError("$dim3 has 3-D maps, only 2-D maps are supported") read_ionex(dim3)

    # map epochs: hour 24 is the next day's 00:00; not increasing is an error
    epochlines(lines) = findall(l -> headerlabel(l) == "EPOCH OF CURRENT MAP", lines)
    h24 = rewritten(day59) do ls
        e = epochlines(ls)[end]
        @assert ls[e][1:36] == "  2021     3     1     0     0     0"
        [ls[1:e-1]; "  2021     2    28    24     0     0" * ls[e][37:end]; ls[e+1:end]]
    end
    @test axiskeys(read_ionex(h24).tec, :epoch) == DateTime(2021, 2, 28, 21):Hour(1):DateTime(2021, 3, 1)
    @test isequal(read_ionex(h24), read_ionex(day59))
    # the epochs of maps 1 and 2 (hours 0 and 1) swapped, or map 1's in both
    @testset for order in [[2, 1, 3, 4, 5, 6, 7], [1, 1, 3, 4, 5, 6, 7]]
        reordered = rewritten(ls -> (e = epochlines(ls); ls[[1:e[1]-1; e[order[1]]; e[1]+1:e[2]-1; e[order[2]]; e[2]+1:end]]), day60)
        @test_throws ArgumentError("$reordered has TEC map epochs not strictly increasing: $(DateTime(2021, 3, 1) .+ Hour.(order .- 1))") read_ionex(reordered)
    end
    # RMS maps at the TEC map epochs: one missing, one at another epoch
    cod = joinpath(Track.DATA, "COD0OPSFIN_20252060000_01D_01H_GIM_crop.INX.gz")
    tecepochs = collect(DateTime(2025, 7, 25, 3):Hour(1):DateTime(2025, 7, 25, 13))
    rmsstarts(lines) = findall(l -> headerlabel(l) == "START OF RMS MAP", lines)
    fewer = rewritten(ls -> [ls[1:rmsstarts(ls)[end]-1]; ls[end]], cod)
    @test_throws ArgumentError("$fewer has RMS maps at $(tecepochs[1:end-1]), its TEC maps at $tecepochs") read_ionex(fewer)
    moved = rewritten(cod) do ls
        e = rmsstarts(ls)[2] + 1
        @assert ls[e][1:36] == "  2025     7    25     4     0     0"
        [ls[1:e-1]; "  2025     7    25     4    30     0" * ls[e][37:end]; ls[e+1:end]]
    end
    @test_throws ArgumentError("$moved has RMS maps at $([tecepochs[1]; DateTime(2025, 7, 25, 4, 30); tecepochs[3:end]]), its TEC maps at $tecepochs") read_ionex(moved)
    # row records that disagree with the header grid
    @test_throws AssertionError read_ionex(rewritten(ls -> replace(ls, "   -10.0 -10.0  55.0   5.0 450.0" * " "^28 * "LAT/LON1/LON2/DLON/H" =>
                                                                       "   -10.0  -5.0  60.0   5.0 450.0" * " "^28 * "LAT/LON1/LON2/DLON/H"), day60))
    @test_throws AssertionError read_ionex(rewritten(ls -> replace(ls, "   -12.5 -10.0  55.0   5.0 450.0" * " "^28 * "LAT/LON1/LON2/DLON/H" =>
                                                                       "   -12.0 -10.0  55.0   5.0 450.0" * " "^28 * "LAT/LON1/LON2/DLON/H"), day60))

    # a node without a value (9999): NaN wherever it takes part
    holed = rewritten(day60) do ls
        j = starts(ls)[1] + 3   # first data line of the first map: lat -10, lon -10 onwards
        [ls[1:j-1]; " 9999" * ls[j][6:end]; ls[j+1:end]]
    end
    m = GlobalMap(holed; mapping=SLM())
    @test findall(isnan, m.ionex.tec) == [CartesianIndex(1, 1, 1)]
    @test isnan(vtec(m, -11.0u"°", -8.0u"°", DateTime(2021, 3, 1)))
    @test vtec(m, -11.0u"°", 0.0u"°", DateTime(2021, 3, 1)) == vtec(GlobalMap(day60; mapping=SLM()), -11.0u"°", 0.0u"°", DateTime(2021, 3, 1))
end

@testitem "IONEX records: `nan` as a missing value" setup=[Track] begin
    using Dates, Uncertain
    using Unitful: °
    using IonosphericRM: read_ionex
    using .Track: TECU
    # UQRG 2022-02-20, maps 06:30–07:00 on 50–40°N, 25–45°E; the 06:45 RMS map has `nan` at 45°N 35°E
    path = joinpath(Track.DATA, "uqrg0510.22i_crop.gz")
    m = GlobalMap(path; mapping=SLM())
    I = m.ionex
    @test findall(isnan, I.rms) == [CartesianIndex(3, 3, 2)]
    @test !any(isnan, I.tec)
    # the 06:45 maps as in the file (0.1 TECU)
    @test I.tec[:, :, 2] ≈ [ 92  98 112 115 124
                            105 114 124 130 139
                            123 135 145 143 140
                            127 133 151 150 167
                            135 135 140 163 170] ./ 10 rtol=1e-14
    @test I.rms[[1, 2, 4, 5], :, 2] ≈ [23 23 24 23 23
                                       23 23 24 23 24
                                       23 24 23 23 24
                                       23 24 24 23 23] ./ 10 rtol=1e-14
    @test I.rms[3, [1, 2, 4, 5], 2] ≈ [2.3, 2.3, 2.3, 2.3] rtol=1e-14

    # at 06:45: VTEC throughout, RMS NaN at the node and the centres of its four cells, finite in other cells
    t = DateTime(2022, 2, 20, 6, 45)
    near = [vtec(m, lat * °, lon * °, t) for (lat, lon) in [(45.0, 35.0), (46.25, 32.5), (46.25, 37.5), (43.75, 32.5), (43.75, 37.5)]]
    @test U.value.(near) ≈ [14.5, 12.95, 13.55, 14.1, 14.725] .* TECU rtol=1e-12
    @test all(isnan, U.uncertainty.(near))
    far = [vtec(m, lat * °, lon * °, t) for (lat, lon) in [(48.75, 27.5), (41.25, 42.5), (50.0, 45.0), (40.0, 25.0)]]
    @test U.value.(far) ≈ [10.225, 16.25, 12.4, 13.5] .* TECU rtol=1e-12
    @test U.uncertainty.(far) ≈ [2.3, 2.325, 2.3, 2.3] .* TECU rtol=1e-12
    # the node is missing only from its own map: NaN within its two map intervals, not at the neighbouring map epochs
    @test isnan(U.uncertainty(vtec(m, 45.0°, 35.0°, DateTime(2022, 2, 20, 6, 50))))
    ends = [vtec(m, 45.0°, 35.0°, t) for t in DateTime(2022, 2, 20, 6) .+ Minute.([30, 60])]
    @test U.value.(ends) ≈ [13.7, 14.2] .* TECU rtol=1e-12
    @test U.uncertainty.(ends) ≈ [2.3, 2.3] .* TECU rtol=1e-12

    # `nan` in any case; any other non-number fails to parse
    function respelled(s)
        out = tempname()
        write(out, replace(read(`gzip -dc $path`, String), "  nan" => s))
        out
    end
    @testset for s in ["  NaN", "  NAN"]
        @test isequal(read_ionex(respelled(s)), I)
    end
    @testset for s in ["  ???", "  inf", "   na", "  9e9"]
        @test_throws ArgumentError read_ionex(respelled(s))
    end
end

@testitem "multi-day maps: incompatible or non-consecutive files" setup=[Track] begin
    using Dates, AxisKeys, Unitful
    using IonosphericRM: read_ionex, append_ionex
    day59, day60 = joinpath.(Track.DATA, ["codg0590.21i_crop.gz", "codg0600.21i_crop.gz"])
    @test_throws ArgumentError GlobalMap([day60, day59]; mapping=MSLM())
    @test_throws ArgumentError GlobalMap([day59, day59]; mapping=MSLM())
    @test_throws ArgumentError GlobalMap([day59, joinpath(Track.DATA, "COD0OPSFIN_20252060000_01D_01H_GIM_crop.INX.gz")]; mapping=MSLM())

    a, b = read_ionex.([day59, day60])
    rekey(A; kw...) = KeyedArray(AxisKeys.keyless_unname(A); named_axiskeys(A)..., kw...)
    @testset for b′ in [(; b..., tec=rekey(b.tec; lat=(-12.5:-2.5:-47.5)u"°")), (; b..., tec=rekey(b.tec; lon=(-5.0:5.0:60.0)u"°")),
                        (; b..., radius=6371.2u"km"), (; b..., height=400u"km"), (; b..., mapping="COSZ")]
        @test_throws ArgumentError append_ionex(a, b′)
    end
    # RMS maps join like TEC maps, and only with RMS maps
    @test append_ionex((; a..., rms=a.tec), (; b..., rms=2b.tec)).rms == cat(a.tec[:, :, 1:3], 2b.tec; dims=3)
    @test_throws MethodError append_ionex(a, (; b..., rms=b.tec))
    # a one-interval gap interpolates like any map interval; a longer one does not
    shifted(I, Δ) = (; I..., tec=rekey(I.tec; epoch=axiskeys(I.tec, :epoch) .+ Δ))
    @test axiskeys(append_ionex(a, shifted(b, Hour(1))).tec, :epoch) == DateTime(2021, 2, 28, 21):Hour(1):DateTime(2021, 3, 1, 7)
    @test_throws ArgumentError("gap between IONEX files: 2021-03-01T00:00:00 to 2021-03-01T02:00:00, longer than any map interval within them") append_ionex(a, shifted(b, Hour(2)))
    # single maps have no interval: any gap
    single(I, k) = (; I..., tec=I.tec[:, :, k:k])
    @test_throws ArgumentError("gap between IONEX files: 2021-03-01T00:00:00 to 2021-03-01T01:00:00, longer than any map interval within them") append_ionex(single(a, 4), single(b, 2))
    # a wider overlap: the later file from its first epoch on
    let b = shifted(b, -Hour(2)), c = append_ionex(a, b)
        @test axiskeys(c.tec, :epoch) == DateTime(2021, 2, 28, 21):Hour(1):DateTime(2021, 3, 1, 4)
        @test c.tec == cat(a.tec[:, :, 1:1], b.tec; dims=3)
    end
end

@testitem "package quality: Aqua, compat bounds" begin
    import Aqua
    Aqua.test_all(IonosphericRM)

    import CompatHelperLocal as CHL
    CHL.@check()
end
