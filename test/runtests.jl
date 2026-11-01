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
    using IonosphericRM: K_RM, ecef, up, pierce, source_direction
    using LinearAlgebra, Dates
    @test K_RM ≈ 2.6312e-13 rtol=1e-4
    o = ecef(Track.SITE)
    @test pierce(o, normalize(o), 6821e3) ≈ 6821e3 * normalize(o)
    # 3C138 elevation, astropy (via spinifex) reference: 37.08263° at 06:06, 30.659° (ALBUS) at 10:12
    el(t) = rad2deg(asin(dot(source_direction(Track.SRC, t), up(Track.SITE))))
    @test el(Track.SCANS[1]) ≈ 37.0826 atol=0.01
    @test el(Track.SCANS[end]) ≈ 30.659 atol=0.02
end

@testitem "global map: verified reference values" setup=[Track] begin
    using Dates, SkyCoords
    path = joinpath(Track.DATA, "COD0OPSFIN_20252060000_01D_01H_GIM_crop.INX.gz")
    # independently verified calculator (vs spinifex 2.0 and AIPS TECOR), CODE final GIM 2025-206
    slm = [-1.0456, -1.1890, -1.3138, -1.4600, -1.5775, -1.7246, -1.8365, -2.0016, -2.1344, -2.2655, -2.3375, -2.4113, -2.4826, -2.5611]
    mslm = [-1.0184, -1.1611, -1.2851, -1.4303, -1.5465, -1.6913, -1.8009, -1.9613, -2.0891, -2.2135, -2.2788, -2.3440, -2.4047, -2.4701]
    @test ionospheric_rm.(Ref(GlobalMap(path)), Ref(Track.SITE), Ref(Track.SRC), Track.SCANS) ≈ slm atol=2e-4   # norm over the 14 scans
    @test ionospheric_rm.(Ref(GlobalMap(path; mapping=MSLM())), Ref(Track.SITE), Ref(Track.SRC), Track.SCANS) ≈ mslm atol=2e-4
    # circumpolar source, time outside the map epochs
    @test_throws ErrorException ionospheric_rm(GlobalMap(path), Track.SITE, ICRSCoords(0.0, -1.2), DateTime(2025, 7, 25, 20))
end

@testitem "horizon and non-IONEX input" setup=[Track] begin
    using Dates
    m = GlobalMap(joinpath(Track.DATA, "COD0OPSFIN_20252060000_01D_01H_GIM_crop.INX.gz"))
    @test isnan(ionospheric_rm(m, Track.SITE, Track.SRC, DateTime(2025, 7, 25, 1)))    # 3C138 below the horizon
    @test_throws Exception GlobalMap(@__FILE__)
end

@testitem "_" begin
    import Aqua
    Aqua.test_all(IonosphericRM)

    import CompatHelperLocal as CHL
    CHL.@check()
end
