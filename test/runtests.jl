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
    using Dates, SkyCoords, Uncertain
    path = joinpath(Track.DATA, "COD0OPSFIN_20252060000_01D_01H_GIM_crop.INX.gz")
    # independently verified calculator (vs spinifex 2.0 and AIPS TECOR), CODE final GIM 2025-206
    slm = [-1.0456, -1.1890, -1.3138, -1.4600, -1.5775, -1.7246, -1.8365, -2.0016, -2.1344, -2.2655, -2.3375, -2.4113, -2.4826, -2.5611]
    mslm = [-1.0184, -1.1611, -1.2851, -1.4303, -1.5465, -1.6913, -1.8009, -1.9613, -2.0891, -2.2135, -2.2788, -2.3440, -2.4047, -2.4701]
    @test U.value.(ionospheric_rm.(Ref(GlobalMap(path; mapping=SLM())), Ref(Track.SITE), Ref(Track.SRC), Track.SCANS)) ≈ slm atol=2e-4   # norm over the 14 scans
    @test U.value.(ionospheric_rm.(Ref(GlobalMap(path; mapping=MSLM())), Ref(Track.SITE), Ref(Track.SRC), Track.SCANS)) ≈ mslm atol=2e-4
    # circumpolar source, time outside the map epochs
    @test_throws BoundsError ionospheric_rm(GlobalMap(path; mapping=MSLM()), Track.SITE, ICRSCoords(0.0, -1.2), DateTime(2025, 7, 25, 20))
    # any sky frame: 3C138 in galactic coordinates
    m = GlobalMap(path; mapping=MSLM())
    gal, icrs = ionospheric_rm.(Ref(m), Ref(Track.SITE), [convert(GalCoords, Track.SRC), Track.SRC], Track.SCANS[1])
    @test U.value(gal) ≈ U.value(icrs) rtol=1e-12
    @test U.uncertainty(gal) ≈ U.uncertainty(icrs) rtol=1e-12
end

@testitem "horizon and non-IONEX input" setup=[Track] begin
    using Dates, Uncertain
    m = GlobalMap(joinpath(Track.DATA, "COD0OPSFIN_20252060000_01D_01H_GIM_crop.INX.gz"); mapping=MSLM())
    @test isnan(U.value(@inferred ionospheric_rm(m, Track.SITE, Track.SRC, DateTime(2025, 7, 25, 1))))    # 3C138 below the horizon
    @test isnan(@inferred ionospheric_rm(GlobalMap(joinpath(Track.DATA, "codg0600.21i_crop.gz"); mapping=MSLM()), Track.SITE, Track.SRC, DateTime(2021, 3, 1, 3)))
    @test_throws Exception GlobalMap(@__FILE__; mapping=MSLM())
    @test_throws UndefKeywordError GlobalMap(joinpath(Track.DATA, "COD0OPSFIN_20252060000_01D_01H_GIM_crop.INX.gz"))
end

@testitem "RMS maps: VTEC and RM uncertainties" setup=[Track] begin
    using Dates, Uncertain
    using IonosphericRM: read_ionex, vtec
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
        I = read_ionex(joinpath(Track.DATA, file))
        v = [vtec(I, t, lat, lon) for (t, _) in ref for (lat, lon) in points]
        expected = reduce(vcat, last.(ref))
        @test U.value.(v) ≈ first.(expected) atol=3e-6   # norm over values given to 1e-6
        @test U.uncertainty.(v) ≈ last.(expected) atol=3e-6
    end
    @test vtec(read_ionex(joinpath(Track.DATA, "codg0600.21i_crop.gz")), DateTime(2021, 3, 1, 3), -30.0, 21.0) isa Float64   # no RMS maps

    # spinifex 2.0 get_rm_from_skycoord, single layer at 450 km (= SLM): RM and the TEC part of its error, |K·RMS·M·B∥|,
    # rescaled from its rounded constant 2.62e-6 to K_RM
    m = GlobalMap(joinpath(Track.DATA, "COD0OPSFIN_20252060000_01D_01H_GIM_crop.INX.gz"); mapping=SLM())
    rm = [ionospheric_rm(m, Track.SITE, Track.SRC, t) for t in DateTime(2025, 7, 25) .+ Minute.([366, 480, 612])]
    @test U.value.(rm) ≈ [-1.045650, -1.843167, -2.561181] atol=3e-4
    @test U.uncertainty.(rm) ≈ [0.105383, 0.106252, 0.136799] atol=3e-5
    @test (@inferred ionospheric_rm(m, Track.SITE, Track.SRC, DateTime(2025, 7, 25, 8))) isa Uncertain.Value
end

@testitem "grid interpolation: longitude wrap, seam, grid edges, map epochs" begin
    using Dates
    using IonosphericRM: vtec
    t0 = DateTime(2025, 1, 1)
    f(lat, lon) = 10 + 3cosd(lon) + 2sind(2lon) + lat / 10
    lats = 10.0:-2.5:-10.0
    # two maps an hour apart, the second rotated 15° west as if fixed to the Sun
    ionex(lons; fixlon=identity) = (; epochs=[t0, t0 + Hour(1)], lats, lons, rms=nothing,
        tec=stack([fixlon(f.(lats, lons')), fixlon(f.(lats, lons' .+ 15))]))

    # a full circle with or without the repeated ±180 column, starting anywhere: the same field
    grids = ionex.([-180.0:5.0:180.0, -180.0:5.0:175.0, 0.0:5.0:355.0, 0.0:5.0:360.0])
    pts = [(t, lat, lon) for t in t0 .+ Minute.([0, 20, 60]) for lat in (7.0, -10.0) for lon in (-180, -179, -1, 0, 1, 176, 177.5, 179, 180, 181, 359)]
    v = [[vtec(I, pt...) for pt in pts] for I in grids]
    @test all(≈(v[1]; rtol=1e-12), v)
    @test vtec(grids[1], t0, 5.0, 177.5) ≈ (f(5, 175) + f(5, 180)) / 2
    @test vtec(grids[2], t0, 5.0, 177.5) ≈ (f(5, 175) + f(5, -180)) / 2
    @test vtec(grids[1], t0 + Hour(1), 5.0, 2.5) ≈ (f(5, 15) + f(5, 20)) / 2

    # ±180 columns that differ: each side of the seam interpolates from its own column
    seam = ionex(-180.0:5.0:180.0; fixlon=M -> [M[:, 1:end-1] M[:, end] .+ 10])
    @test vtec(seam, t0, 5.0, 179.0) ≈ 0.2f(5, 175) + 0.8(f(5, 180) + 10)
    @test vtec(seam, t0, 5.0, -179.0) ≈ 0.8f(5, -180) + 0.2f(5, -175)

    # regional grid: grid edges included, beyond them an error, not an extrapolation
    reg = ionex(-10.0:5.0:55.0)
    @test vtec(reg, t0, -10.0, 55.0) ≈ f(-10, 55)
    @test vtec(reg, t0, 10.0, -10.0) ≈ f(10, -10)
    @test_throws DomainError vtec(reg, t0, 0.0, 56.0)
    @test_throws DomainError vtec(reg, t0, 0.0, -11.0)
    @test_throws DomainError vtec(reg, t0 + Minute(30), 0.0, 50.0)   # rotated out of the grid in the second map
    @test_throws DomainError vtec(grids[1], t0, 10.5, 0.0)
    @test_throws DomainError vtec(grids[1], t0, -12.0, 0.0)

    # map epochs: both ends included, a file with a single map
    @test vtec(grids[1], t0 + Hour(1), 0.0, 20.0) ≈ f(0, 35)
    single = (; grids[1]..., epochs=[t0], tec=grids[1].tec[:, :, 1:1])
    @test vtec(single, t0, 0.0, 20.0) ≈ f(0, 20)
    @test_throws BoundsError vtec(single, t0 + Minute(1), 0.0, 20.0)
    @test_throws BoundsError vtec(grids[1], t0 - Minute(1), 0.0, 20.0)
    @test_throws BoundsError vtec(grids[1], t0 + Minute(61), 0.0, 20.0)
end

@testitem "grid interpolation: real maps across the ±180 seam" setup=[Track] begin
    using Dates, Uncertain
    using IonosphericRM: read_ionex, vtec
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
        I = read_ionex(joinpath(Track.DATA, file))
        v = [vtec(I, t, lat, lon) for (t, _) in ref for (lat, lon) in points]
        expected = reduce(vcat, last.(ref))
        @test U.value.(v) ≈ first.(expected) atol=3e-6
        @test U.uncertainty.(v) ≈ last.(expected) atol=3e-6
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
    @test_throws BoundsError vtec(m.ionex, DateTime(2021, 3, 1, 6, 30), -30.0, 21.0)
    @test_throws BoundsError vtec(m.ionex, DateTime(2021, 2, 28, 20, 30), -30.0, 21.0)
end

@testitem "multi-day maps: EXPONENT records apply per file and per map" setup=[Track] begin
    # the 03-01 file rewritten at exponent -2, stated in a record inside each map (header still -1)
    day59, day60, day60e = joinpath.(Track.DATA, ["codg0590.21i_crop.gz", "codg0600.21i_crop.gz", "codg0600.21i_crop_exp2.gz"])
    a = GlobalMap([day59, day60]; mapping=MSLM()).ionex; b = GlobalMap([day59, day60e]; mapping=MSLM()).ionex
    @test b.epochs == a.epochs
    @test b.tec ≈ a.tec rtol=1e-14
end

@testitem "IONEX records: fixed-width fields, record order, completeness, missing values" setup=[Track] begin
    using Dates
    using IonosphericRM: read_ionex, headerlabel, vtec
    day60, day60e2, day60e3 = joinpath.(Track.DATA, ["codg0600.21i_crop.gz", "codg0600.21i_crop_exp2.gz", "codg0600.21i_crop_exp3.gz"])
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
    # cut between maps or inside a map
    @test_throws ArgumentError read_ionex(rewritten(ls -> ls[1:starts(ls)[end]-1], day60))
    @test_throws Exception read_ionex(rewritten(ls -> ls[1:starts(ls)[end]+5], day60))
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
    I = read_ionex(holed)
    @test findall(isnan, I.tec) == [CartesianIndex(1, 1, 1)]
    @test isnan(vtec(I, DateTime(2021, 3, 1), -11.0, -8.0))
    @test vtec(I, DateTime(2021, 3, 1), -11.0, 0.0) == vtec(ref, DateTime(2021, 3, 1), -11.0, 0.0)
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
    # RMS maps join like TEC maps, and only with RMS maps
    @test append_ionex((; a..., rms=a.tec), (; b..., rms=2b.tec)).rms == cat(a.tec[:, :, 1:3], 2b.tec; dims=3)
    @test_throws MethodError append_ionex(a, (; b..., rms=b.tec))
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
