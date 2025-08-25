project "engine"
    location "."
    kind "StaticLib"
    language "C++"
    cppdialect "C++26"

    targetdir "../../bin/%{cfg.buildcfg}/engine"
    objdir "../../bin-int/%{cfg.buildcfg}/engine"

    includedirs {
        "../../src"
    }

    files {
        "../../src/engine/**.hpp", "../../src/engine/**.cpp",
        "../../src/geometry/**.hpp", "../../src/geometry/**.cpp",
        "../../src/physics/**.hpp", "../../src/physics/**.cpp",
        "../../src/utils/**.hpp", "../../src/utils/**.cpp"
    }

    filter "configurations:Debug"
        defines { "DEBUG" }
        symbols "On"
        runtime "Debug"

    filter "configurations:Release"
        defines { "NDEBUG" }
        optimize "Speed"
        runtime "Release"

    filter "configurations:Dist"
        defines { "NDEBUG", "DIST" }
        optimize "Full"
        runtime "Release"

    filter "configurations:Profile"
        defines { "NDEBUG", "PROFILE" }
        symbols "On"
        optimize "Debug"
        runtime "Release"

    filter {}
