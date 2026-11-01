using TestItems
using TestItemRunner
@run_package_tests

@testmodule Track begin
    using Dates, SkyCoords
    # MeerKAT array centre and 3C138 (J2000); MeerKAT 2025-07-25 scan mid-times (UTC)
    const DATA = joinpath(@__DIR__, "data")
    const SITE = (; lat=deg2rad(-30.71292), lon=deg2rad(21.44380), h=1050.0)
    const SRC = ICRSCoords(deg2rad(15 * (5 + 21/60 + 9.886/3600)), deg2rad(16 + 38/60 + 22.05/3600))
    const SCANS = DateTime(2025, 7, 25) .+ Minute.([366, 386, 403, 424, 441, 462, 479, 500, 517, 537, 555, 574, 593, 612])
end

@testitem "constants and geometry" setup=[Track] begin
    using IonosphericRM: K_RM, ecef, geodetic, up, pierce, source_direction
    using LinearAlgebra, Dates
    @test K_RM ≈ 2.6312e-13 rtol=1e-4
    g = geodetic(ecef(Track.SITE))
    @test [g.lat, g.lon, g.h] ≈ [Track.SITE.lat, Track.SITE.lon, Track.SITE.h] atol=1e-6
    o = ecef(Track.SITE)
    @test pierce(o, normalize(o), 6821e3) ≈ 6821e3 * normalize(o)
    # 3C138 elevation, astropy (via spinifex) reference: 37.08263° at 06:06, 30.659° (ALBUS) at 10:12
    el(t) = rad2deg(asin(dot(source_direction(Track.SRC, t), up(Track.SITE))))
    @test el(Track.SCANS[1]) ≈ 37.0826 atol=0.01
    @test el(Track.SCANS[end]) ≈ 30.659 atol=0.02
end

@testitem "IONEX DCB block" setup=[Track] begin
    using IonosphericRM: read_ionex_dcb
    dcb = read_ionex_dcb(joinpath(Track.DATA, "COD0OPSFIN_20252060000_01D_01H_GIM_crop.INX.gz"))
    @test dcb.ref == "C1W-C2W"
    @test [dcb.sat[1], dcb.sat[32], dcb.sta["SUTM"], dcb.sta["SUTH"]] == [3.408, -4.787, 10.368, -0.899]
end

@testitem "SP3 interpolation" setup=[Track] begin
    using IonosphericRM: read_sp3, satpos
    using Dates
    sp3 = read_sp3(joinpath(Track.DATA, "sutm", "COD0OPSFIN_20252060000_crop_ORB.SP3.gz"))
    @test sp3.times[[1, end]] == [DateTime(2025, 7, 25, 2, 30), DateTime(2025, 7, 25, 13, 30)]
    i = 40
    @test satpos(sp3, 5, sp3.times[i]) ≈ sp3.pos[5][i] atol=1e-6
    # leave-one-out: node i predicted from its neighbours spaced 10 min (every other node) within 1 m
    thinned = (; times=sp3.times[1:2:end], pos=Dict(p => v[1:2:end] for (p, v) in sp3.pos))
    @test satpos(thinned, 5, sp3.times[2i]) ≈ sp3.pos[5][2i] atol=1
end

@testitem "global map: verified reference values" setup=[Track] begin
    using Dates, SkyCoords
    path = joinpath(Track.DATA, "COD0OPSFIN_20252060000_01D_01H_GIM_crop.INX.gz")
    @test tec_model(path) isa GlobalMap
    # independently verified calculator (vs spinifex 2.0 and AIPS TECOR), CODE final GIM 2025-206
    slm = [-1.0456, -1.1890, -1.3138, -1.4600, -1.5775, -1.7246, -1.8365, -2.0016, -2.1344, -2.2655, -2.3375, -2.4113, -2.4826, -2.5611]
    mslm = [-1.0184, -1.1611, -1.2851, -1.4303, -1.5465, -1.6913, -1.8009, -1.9613, -2.0891, -2.2135, -2.2788, -2.3440, -2.4047, -2.4701]
    @test ionospheric_rm.(Ref(tec_model(path)), Ref(Track.SITE), Ref(Track.SRC), Track.SCANS) ≈ slm atol=2e-4   # norm over the 14 scans
    @test ionospheric_rm.(Ref(tec_model(path; mapping=MSLM())), Ref(Track.SITE), Ref(Track.SRC), Track.SCANS) ≈ mslm atol=2e-4
    # circumpolar source, time outside the map epochs
    @test_throws ErrorException ionospheric_rm(tec_model(path), Track.SITE, ICRSCoords(0.0, -1.2), DateTime(2025, 7, 25, 20))
end

@testitem "single station SUTM: ALBUS RI_G01 and full-day reference" setup=[Track] begin
    m = tec_model(joinpath(Track.DATA, "sutm"))
    @test m isa SingleStation{MeanVTEC}
    rm = ionospheric_rm.(Ref(m), Ref(Track.SITE), Ref(Track.SRC), Track.SCANS)
    # ALBUS (2025-11 master) RI_G01, station sutm, run on the same day of data
    albus = [-0.83, -0.986, -1.0896, -1.176, -1.3376, -1.5031, -1.6155, -1.6724, -1.8377, -1.9815, -2.1121, -2.2434, -2.369, -2.5752]
    @test maximum(abs, rm - albus) < 0.035
    # this package on the full-day RINEX file (fixture is cropped to 03:30-12:30)
    full = [-0.8451, -0.9995, -1.1062, -1.2044, -1.3598, -1.5259, -1.639, -1.6962, -1.8617, -2.0053, -2.141, -2.2634, -2.3844, -2.5898]
    @test maximum(abs, rm - full) < 0.006
end

@testitem "single station SUTM: planar VTEC" setup=[Track] begin
    m = tec_model(joinpath(Track.DATA, "sutm"); model=PlanarVTEC())
    full = [-1.0697, -1.2463, -1.3671, -1.4845, -1.5721, -1.7175, -1.8378, -2.0233, -2.1082, -2.2588, -2.4381, -2.6181, -2.774, -2.8374]
    @test maximum(abs, ionospheric_rm.(Ref(m), Ref(Track.SITE), Ref(Track.SRC), Track.SCANS) - full) < 0.04
end

@testitem "missing data and horizon" setup=[Track] begin
    using Dates
    m = tec_model(joinpath(Track.DATA, "sutm"))
    @test isnan(ionospheric_rm(m, Track.SITE, Track.SRC, DateTime(2025, 7, 25, 12, 45)))   # source up, after the data
    @test isnan(ionospheric_rm(m, Track.SITE, Track.SRC, DateTime(2025, 7, 25, 1)))    # 3C138 below the horizon
    @test_throws Exception tec_model(@__FILE__)
end

@testitem "_" begin
    import Aqua
    Aqua.test_all(IonosphericRM)

    import CompatHelperLocal as CHL
    CHL.@check()
end
