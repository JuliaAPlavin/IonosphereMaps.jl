# Interactive globe: the Earth lit by the Sun, the IONEX ionosphere shell coloured by VTEC, and lines of sight from
# site presets toward a source picked on a sky map, through their shell pierce points, coloured by slant TEC.
# Scene: Earth-fixed, unit = the IONEX base radius, +z = north pole.
# Run: julia --project=examples examples/ionosphere_viewer.jl [--mapping MSLM] ionex...   (--help for details)

using VolumeRendering, ColorTypes, ColorSchemes
using IonosphereMaps, RectiGrids, Dates, Unitful
import UnitfulAngles   # u"hourAngle"
import UnitfulAstro    # u"TECU"
using SkyCoords: ICRSCoords, CartesianCoords
using StaticArrays: SVector
using IntervalSets: (..)
using LinearAlgebra: normalize
import Uncertain
using NaturalEarth: naturalearth
using SatelliteToolboxTransformations: geodetic_to_ecef, ecef_to_geodetic, r_eci_to_ecef, r_eci_to_eci, MOD, PEF, J2000
using SatelliteToolboxCelestialBodies: sun_position_mod
import GeoInterface as GI
import CImGui, ImPlot, GLFW, ModernGL     # GLFW + ModernGL activate the CImGui GLFW/OpenGL backend
using CImGui: ImVec2, ImVec4
using ArgParse: ArgParseSettings, @add_arg_table!, parse_args

# ECEF positions, drawn as directions on the spherical globe
coastlines() = map(naturalearth("coastline", 110)) do feature
    map(p -> geodetic_to_ecef(GI.y(p) * u"°", GI.x(p) * u"°", 0), GI.getpoint(feature.geometry))
end

# ── sites ──

# ITRF position (m) → geodetic WGS84 site for `LineOfSight`
itrf(name, x, y, z) = let (lat, lon, h) = ecef_to_geodetic(SVector(x, y, z))
    (; name, lat, lon, h = h * u"m")
end

# ITRF positions from the EVN Observation Planner catalogue (vlbiplanobs 5.1.1, data/stations_catalog.inp);
# EVN: its default antennas (data/network_catalog.inp); VLA, MeerKAT: phased arrays
const SITE_PRESETS = [
    "VLBA" => [
        itrf("Sc", 2607848.6379, -5488069.5358, 1932739.7326),  # St. Croix
        itrf("Hn", 1446374.8658, -4447939.6774, 4322306.1838),  # Hancock
        itrf("Nl", -130872.4987, -4762317.0925, 4226851.0014),  # North Liberty
        itrf("Fd", -1324009.3266, -5332181.9547, 3231962.3949),  # Fort Davis
        itrf("La", -1449752.5839, -4975298.5757, 3709123.8459),  # Los Alamos
        itrf("Pt", -1640953.9383, -5014816.0237, 3575411.7916),  # Pie Town
        itrf("Kp", -1995678.8571, -5037317.7049, 3357328.0174),  # Kitt Peak
        itrf("Ov", -2409150.4018, -4478573.118, 3838617.3385),  # Owens Valley Radio Observatory
        itrf("Br", -2112065.2062, -3705356.5048, 4726813.6759),  # Brewster
        itrf("Mk", -5464075.1874, -2495248.1168, 2148297.3643),  # Mauna Kea
    ],
    "EVN" => [
        itrf("Ef", 4033947.2616, 486990.7943, 4900431.0017),  # Effelsberg
        itrf("Hh", 5085442.7637, 2668263.7913, -2768696.7523),  # Hartebeesthoek
        itrf("Jb2", 3822846.76, -153802.28, 5086285.9),  # Jodrell Bank (Mrk 2)
        itrf("Mc", 4461369.6954, 919597.124, 4449559.3812),  # Medicina
        itrf("Nt", 4934562.8407, 1321201.5428, 3806484.7336),  # Noto
        itrf("Ur", 228310.2129, 4631922.7656, 4367064.0638),  # Urumqi (Nanshan)
        itrf("O8", 3370965.8787, 711466.1978, 5349664.2006),  # Onsala 25 m
        itrf("O6", 3370605.7035, 711917.8146, 5349830.9852),  # Onsala 20 m
        itrf("Sr", 4865183.542, 791922.255, 4035136.024),  # Sardinia
        itrf("T6", -2826708.603, 4679237.077, 3274667.551),  # Tianma
        itrf("Tr", 3638558.51, 1221969.72, 5077036.76),  # Torun
        itrf("Ys", 4848761.8027, -261484.1746, 4123085.0529),  # Yebes
        itrf("Wb", 3828767.2647, 442446.1739, 5064921.57),  # Westerbork
        itrf("Ir", 3183649.314, 1276902.989, 5359264.71),  # Irbene
    ],
    "GBT" => [
        itrf("Gb", 882589.4102, -4924872.3416, 3943729.4062),  # Green Bank Telescope
    ],
    "VLA" => [
        itrf("Y27", -1601185.4286, -5041977.1754, 3554875.6231),  # VLA, phased array
    ],
    "MeerKAT" => [
        itrf("Me", 5109360.13332123, 2006852.58604291, -3238948.12747888),  # MeerKAT, phased array
    ],
]

# ── ionosphere ──

tecu(tec) = ustrip(u"TECU", Uncertain.value(tec))

# VTEC (TECU) at UTC t for a lon/lat GridField, sampled through `vtec`. The grid holds the nodes of standard
# 5° × 2.5° GIMs (and the poles), so the GPU's linear interpolation reproduces `vtec`: exactly in latitude;
# in longitude every 4 min (nodes rotate 0.25°/min with the Sun), otherwise within ~0.1 TECU.
vtec_grid(m, t) = map(p -> tecu(vtec(m, p, t)), grid(lon = (-180:1.0:179)u"°", lat = (-90:2.5:90)u"°"))

# ── Sun ──

# MOD position (Vallado) at UTC t; TDB ≈ UTC moves it by ~0.001°
sun_mod(t) = sun_position_mod(datetime2julian(t))

# Earth-fixed unit vector toward the Sun, UT1 ≈ UTC as the package takes it for sources
sun_direction(t) = normalize(r_eci_to_ecef(MOD(), PEF(), datetime2julian(t)) * sun_mod(t))

# J2000 ≈ ICRS, as the package takes it
sun_icrs(t) = convert(ICRSCoords, CartesianCoords{ICRSCoords}(r_eci_to_eci(MOD(), J2000(), datetime2julian(t)) * sun_mod(t)))

# cosine of the Sun's zenith angle on the spherical Earth: > 0 on the day side
sunlight(t) = GLSLField("""
    float value(vec3 p) { return dot(normalize(p), sun); }
    float stepSize(vec3 p) { return 1.0; }   // never marched: colours a surface only
    """; uniforms = (; sun = sun_direction(t)))

# ── lines of sight ──

const RAY_END = 1.5          # geocentric radius where lines of sight end, past the shell
const MARKER_RADIUS = 0.006  # sites and pierce points, ~40 km

scene(m, x) = NoUnits.(x ./ m.ionex.radius)

# per site toward `source` at UTC t: scene positions of the site, the shell pierce point and the line end, and the
# slant TEC (TECU); TEC is NaN below the site's horizon, or where a missing map node enters
function sightlines(m, t, sites, source)
    (; radius, height) = m.ionex
    map(sites) do site
        ray = Ray(LineOfSight(site, source), t)
        point(r) = scene(m, IonosphereMaps.pierce(ray.origin, ray.direction, r))
        (; site = scene(m, ray.origin), pp = point(radius + height), far = point(RAY_END * radius), tec = tecu(ionospheric_tec(m, ray, t)))
    end
end

marker(center, color) = Surface(Sphere(; center, radius = MARKER_RADIUS); color)

# site markers; lines of sight coloured by slant TEC through `tecmap` and their pierce points, where TEC is defined
function sight_marks(m, t, sites, source, tecmap)
    lines = sightlines(m, t, sites, source)
    seen = filter(l -> !isnan(l.tec), lines)
    [map(l -> marker(l.site, RGB(1.0, 0.75, 0.2)), lines);
     Lines(map(l -> [l.site, l.pp, l.far], seen); color = Mapped(map(l -> fill(l.tec, 3), seen), tecmap), width_px = 2.0,
           label = "lines of sight");
     map(l -> marker(l.pp, RGB(1.0, 1.0, 1.0)), seen)]
end

# ── scene ──

function globe_view(m, t, sites, source)
    (; tec, radius, height) = m.ionex
    overlay = Sphere(radius = 1.002)  # lines 13 km above the surface, so they don't tie with it in depth
    # the maps' maximum, which interpolation never exceeds: one window for all times, shared with the lines of sight
    tecmap = Colormap(ColorSchemes.turbo, 0.0..maximum(filter(!isnan, tec)))
    ionosphere = GridField(vtec_grid(m, t))
    # night to day across the Sun ±6° about the horizon (civil twilight)
    daylight = Mapped(sunlight(t), Colormap(ColorScheme([RGB(0.02, 0.025, 0.035), RGB(0.32, 0.36, 0.44)]), -sin(6u"°")..sin(6u"°")))
    base = Any[
        Surface(Sphere(radius = 1.0); color = daylight, label = "Earth"),
        Lines(graticule(overlay; parallels = 5, meridians = 12); color = RGB(0.42, 0.45, 0.52)),
        Lines(on_sphere(coastlines(), overlay); color = RGB(0.75, 0.80, 0.85)),
        Surface(Sphere(radius = (radius + height) / radius); color = Mapped(ionosphere, tecmap), opacity = 0.5, label = "ionosphere"),
    ]
    g = (; view = View(copy(base); camera = Camera(; radius = 1.1)), base, ionosphere, daylight, tecmap)
    setsight!(g, m, t, sites, source)
    g
end

function settime!(g, m, t)
    g.ionosphere.data = vtec_grid(m, t)   # a new array: textures upload per array
    g.daylight.source = sunlight(t)
end

setsight!(g, m, t, sites, source) = (g.view.marks = [g.base; sight_marks(m, t, sites, source, g.tecmap)])

# ── UI ──

# one button per preset, the selected one highlighted; whether the selection changed
function preset_buttons!(selected::Ref{Int}, names)
    changed = false
    for (i, name) in enumerate(names)
        i > 1 && CImGui.SameLine()
        active = selected[] == i
        active && CImGui.PushStyleColor(CImGui.ImGuiCol_Button, unsafe_load(CImGui.GetStyleColorVec4(CImGui.ImGuiCol_ButtonActive)))
        CImGui.Button(name) && !active && (selected[] = i; changed = true)
        active && CImGui.PopStyleColor()
    end
    changed
end

# RA (h, increasing leftward) × Dec (°), the source a drag point (white), the Sun yellow; whether the source moved
function skymap!(ra::Ref{Cdouble}, dec::Ref{Cdouble}, sun::ICRSCoords)
    ImPlot.BeginPlot("##sky", ImVec2(-1, -1), ImPlot.ImPlotFlags_NoMenus | ImPlot.ImPlotFlags_NoBoxSelect | ImPlot.ImPlotFlags_NoLegend) || return false
    ImPlot.SetupAxes("RA (h)", "Dec (°)", ImPlot.ImPlotAxisFlags_Invert, 0)
    ImPlot.SetupAxesLimits(0, 24, -90, 90, ImPlot.ImPlotCond_Always)
    yellow = ImVec4(1.0, 0.85, 0.2, 1.0)
    ImPlot.PlotScatter("Sun", [ustrip(u"hourAngle", sun.ra)], [ustrip(u"°", sun.dec)];
                       spec = ImPlot.ImPlotSpec(Marker = ImPlot.ImPlotMarker_Circle, MarkerSize = 7, MarkerFillColor = yellow, MarkerLineColor = yellow))
    moved = ImPlot.DragPoint(0, ra, dec, ImVec4(1.0, 1.0, 1.0, 1.0), 6)
    ImPlot.EndPlot()
    moved && (ra[] = clamp(ra[], 0, 24); dec[] = clamp(dec[], -90, 90))
    moved
end

const MAPPINGS = (; SLM = SLM(), MSLM = MSLM(), Petrov2023 = Petrov2023())

function cli()
    s = ArgParseSettings(prog = basename(@__FILE__), description = "The Earth with the ionosphere of IONEX global maps, and lines of sight from sites toward a source.")
    @add_arg_table! s begin
        "ionex"
            help = "IONEX global maps (plain, .gz or .Z); several files are consecutive periods in time order, " *
                   "e.g. CODE daily GIMs from https://www.aiub.unibe.ch/download/CODE/<year>/"
            nargs = '+'
            required = true
        "--mapping"
            help = "mapping from VTEC to the slant TEC of the lines of sight: $(join(keys(MAPPINGS), ", "))"
            arg_type = Symbol
            range_tester = in(keys(MAPPINGS))
            default = :MSLM
    end
    s
end

function (@main)(args)
    opts = parse_args(args, cli())   # on a usage error or --help: prints and exits
    m = GlobalMap(opts["ionex"]; mapping = MAPPINGS[opts["mapping"]])
    epochs = axiskeys(m.ionex.tec, :epoch)
    # UI state outlives frames, so widgets don't snap back mid-drag
    minutes = Ref(Cint(0))   # since the first map
    preset = Ref(1)
    ra, dec = Ref(12.4852), Ref(2.0524)   # hourAngle, °: as ImPlot shows them; 3C 273
    t() = first(epochs) + Minute(minutes[])
    sites() = last(SITE_PRESETS[preset[]])
    source() = ICRSCoords(ra[] * u"hourAngle", dec[] * u"°")
    g = globe_view(m, t(), sites(), source())
    CImGui.set_backend(:GlfwOpenGL3)
    ctx = CImGui.CreateContext()
    ImPlot.CreateContext()
    ImPlot.SetImGuiContext(ctx)
    io = CImGui.GetIO()
    io.ConfigFlags = unsafe_load(io.ConfigFlags) | CImGui.ImGuiConfigFlags_DockingEnable
    CImGui.render(ctx; window_title = "Globe", window_size = (1300, 900), opengl_version = v"4.1") do  # shaders are GLSL 4.1
        CImGui.DockSpaceOverViewport(0, CImGui.GetMainViewport(), CImGui.ImGuiDockNodeFlags_PassthruCentralNode)
        CImGui.Begin("Sky", C_NULL)
        moved = skymap!(ra, dec, sun_icrs(t()))
        CImGui.End()
        CImGui.Begin("Globe", C_NULL)
        picked = preset_buttons!(preset, first.(SITE_PRESETS))
        CImGui.SameLine()
        CImGui.SetNextItemWidth(-1)
        # the time is the format text itself: no % specifier, shown verbatim
        timed = CImGui.SliderInt("##time", minutes, 0, (last(epochs) - first(epochs)) ÷ Minute(1), Dates.format(t(), "yyyy-mm-dd HH:MM") * " UTC")
        timed && settime!(g, m, t())
        (moved | picked | timed) && setsight!(g, m, t(), sites(), source())
        VolumeRendering.ShowView(g.view)
        CImGui.End()
    end
end
