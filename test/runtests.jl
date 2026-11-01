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
    @test ionospheric_rm.(Ref(GlobalMap(path; mapping=SLM())), Ref(Track.SITE), Ref(Track.SRC), Track.SCANS) ≈ slm atol=2e-4   # norm over the 14 scans
    @test ionospheric_rm.(Ref(GlobalMap(path; mapping=MSLM())), Ref(Track.SITE), Ref(Track.SRC), Track.SCANS) ≈ mslm atol=2e-4
    # circumpolar source, time outside the map epochs
    @test_throws ErrorException ionospheric_rm(GlobalMap(path; mapping=MSLM()), Track.SITE, ICRSCoords(0.0, -1.2), DateTime(2025, 7, 25, 20))
end

@testitem "horizon and non-IONEX input" setup=[Track] begin
    using Dates
    m = GlobalMap(joinpath(Track.DATA, "COD0OPSFIN_20252060000_01D_01H_GIM_crop.INX.gz"); mapping=MSLM())
    @test isnan(ionospheric_rm(m, Track.SITE, Track.SRC, DateTime(2025, 7, 25, 1)))    # 3C138 below the horizon
    @test_throws Exception GlobalMap(@__FILE__; mapping=MSLM())
    @test_throws UndefKeywordError GlobalMap(joinpath(Track.DATA, "COD0OPSFIN_20252060000_01D_01H_GIM_crop.INX.gz"))
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
    using Dates, SkyCoords
    using IonosphericRM: vtec
    # CODE final GIMs for 2021-02-28 (21–24 h) and 2021-03-01 (0–6 h), cropped; both hold a 03-01 00:00 map
    day59, day60 = joinpath.(Track.DATA, ["codg0590.21i_crop.gz", "codg0600.21i_crop.gz"])
    m59 = GlobalMap(day59; mapping=MSLM()); m60 = GlobalMap(day60; mapping=MSLM())
    m = GlobalMap([day59, day60]; mapping=MSLM())
    @test m.ionex.epochs == DateTime(2021, 2, 28, 21):Hour(1):DateTime(2021, 3, 1, 6)
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
    @test [vtec(m.ionex, t, lat, lon) for (t, _) in spinifex for (lat, lon) in points] ≈ reduce(vcat, last.(spinifex)) atol=3e-6   # norm over 28 values given to 1e-6

    # 3C286: a single daily file jumps at midnight, the merged map does not
    src = ICRSCoords(deg2rad(15 * (13 + 31/60 + 8.288/3600)), deg2rad(30 + 30/60 + 32.96/3600))
    rm(model, t) = ionospheric_rm(model, Track.SITE, src, t)
    midnight = DateTime(2021, 3, 1)
    @test [rm(m59, midnight), rm(m60, midnight)] ≈ [-0.4833, -0.4795] atol=1e-4
    @test rm(m, midnight - Millisecond(1)) ≈ rm(m, midnight) atol=1e-7
    @test rm(m, midnight - Hour(1)) == rm(m59, midnight - Hour(1))
    @test rm.(Ref(m), midnight .+ Minute.(0:30:360)) == rm.(Ref(m60), midnight .+ Minute.(0:30:360))
    @test_throws ErrorException vtec(m.ionex, DateTime(2021, 3, 1, 6, 30), -30.0, 21.0)
    @test_throws ErrorException vtec(m.ionex, DateTime(2021, 2, 28, 20, 30), -30.0, 21.0)
end

@testitem "multi-day maps: EXPONENT records apply per file and per map" setup=[Track] begin
    # the 03-01 file rewritten at exponent -2, stated in a record inside each map (header still -1)
    day59, day60, day60e = joinpath.(Track.DATA, ["codg0590.21i_crop.gz", "codg0600.21i_crop.gz", "codg0600.21i_crop_exp2.gz"])
    a = GlobalMap([day59, day60]; mapping=MSLM()).ionex; b = GlobalMap([day59, day60e]; mapping=MSLM()).ionex
    @test b.epochs == a.epochs
    @test b.tec ≈ a.tec rtol=1e-14
end

@testitem "multi-day maps: incompatible or non-consecutive files" setup=[Track] begin
    using Dates
    using IonosphericRM: read_ionex, append_ionex
    day59, day60 = joinpath.(Track.DATA, ["codg0590.21i_crop.gz", "codg0600.21i_crop.gz"])
    @test_throws ArgumentError GlobalMap([day60, day59]; mapping=MSLM())
    @test_throws ArgumentError GlobalMap([day59, day59]; mapping=MSLM())
    @test_throws ArgumentError GlobalMap([day59, joinpath(Track.DATA, "COD0OPSFIN_20252060000_01D_01H_GIM_crop.INX.gz")]; mapping=MSLM())

    a, b = read_ionex.([day59, day60])
    @testset for (k, v) in [:lats => -10.0:-2.5:-42.5, :lons => -5.0:5.0:55.0, :radius => 6371.2e3, :height => 400e3, :mapping => "COSZ"]
        @test_throws ArgumentError append_ionex(a, (; b..., k => v))
    end
    # a one-interval gap interpolates like any map interval; a longer one does not
    @test append_ionex(a, (; b..., epochs=b.epochs .+ Hour(1))).epochs == DateTime(2021, 2, 28, 21):Hour(1):DateTime(2021, 3, 1, 7)
    @test_throws ArgumentError append_ionex(a, (; b..., epochs=b.epochs .+ Hour(2)))
    # a wider overlap: the later file from its first epoch on
    let c = append_ionex(a, (; b..., epochs=b.epochs .- Hour(2)))
        @test c.epochs == DateTime(2021, 2, 28, 21):Hour(1):DateTime(2021, 3, 1, 4)
        @test c.tec == cat(a.tec[:, :, 1:1], b.tec; dims=3)
    end
end

@testitem "_" begin
    import Aqua
    Aqua.test_all(IonosphericRM)

    import CompatHelperLocal as CHL
    CHL.@check()
end
