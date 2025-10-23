---@diagnostic disable: undefined-global
require "cmake"  -- enable Premake's CMake exporter

workspace "Renderer"
    configurations { "Debug", "Profile", "Dist" }
    architecture "x64"
    location "build"
    filter "action:cmake"
        location "."
    filter {}
    filter "system:windows"
        systemversion "latest"
        defines { "_WINDOWS", "WIN32_LEAN_AND_MEAN", "NOMINMAX", "GLEW_STATIC", "_GLFW_WIN32" }
        links(winlibs)
    filter {}

-- Options
newoption {
    trigger = "scene",
    value = "SCENE_NAME",
    description = "Single .cpp file to build from src/render-projects"
}
local configFile = "config/build_config.lua"
local config = dofile(configFile)
if not config then error("Failed to read/parse config file: " .. configFile) end

local dialect = config.cpp_dialect or "C++23"
local project_dirs = config.projects_dir or "src/render-projects"

local selectedScene = _OPTIONS["scene"]
if selectedScene == "tests" then error("Scene name 'tests' is reserved.") end
if selectedScene == "core" then error("Scene name 'core' is reserved.") end


local inc = {
    ".",
    "src/core/include",
    "src/core/include/utils",
    "src/core/include/application",
    "src/core/include/engine",
    "src/core/include/engine/sdf-rendering",
    "src/core/include/geometry",
    "src/core/include/physics",
    "src/core/include/file-management",
    "src/core/include/openglAPI",

    "external/glew-2.1.0/include",
    "external/glfw-3.4/include",
    "external/glm-0.9.7.1",
    "external/spdlog-1.x/include",
}


project "engine"
    location "build/core"
    kind "StaticLib"
    language "C++"
    cppdialect(dialect)
    staticruntime "off"
    targetdir ("build/core/%{cfg.architecture}/bin")
    objdir    ("build/core/%{cfg.architecture}/obj")

    files {
        "src/core/**.hpp",
        "src/core/**.cpp",
        "external/glew-2.1.0/src/glew.c",
        "external/glfw-3.4/src/**.c",
        "external/glm-0.9.7.1/**.hpp",
    }
    removefiles { "src/core/**_dep.**" }

    includedirs(inc)

    filter "configurations:Debug"
        defines { "DEBUG" }
        runtime "Debug"
        symbols "On"
    filter "configurations:Profile"
        defines { "NDEBUG" }
        runtime "Release"
        symbols "On"
        optimize "Speed"
        linktimeoptimization "On"
    filter "configurations:Dist"
        defines { "NDEBUG" }
        runtime "Release"
        symbols "Off"
        optimize "Speed"
        linktimeoptimization "On"
    filter {}

if selectedScene then

    scenePath = path.join(project_dirs, selectedScene .. ".cpp")

    project(selectedScene)
        location ("build/" .. selectedScene)
        language "C++"
        cppdialect(dialect)
        staticruntime "off"
        kind "ConsoleApp"

        links { "engine" } -- Link the static library


        targetdir ("build/" .. selectedScene .. "/bin/build-%{cfg.architecture}")
        objdir    ("build/" .. selectedScene .. "/obj/build-%{cfg.architecture}")

        files { scenePath }
        includedirs(inc)

        filter "configurations:Debug"
            defines { "DEBUG" }
            runtime "Debug"
            symbols "On"
        filter "configurations:Profile"
            defines { "NDEBUG" }
            runtime "Release"
            symbols "On"
            optimize "Speed"
            linktimeoptimization "On"
        filter "configurations:Dist"
            defines { "NDEBUG" }
            runtime "Release"
            symbols "Off"
            optimize "Speed"
            linktimeoptimization "On"
        filter {}
end



project "unitTests"
    location ("build/tests")
    language "C++"
    cppdialect(dialect)
    staticruntime "off"
    kind "ConsoleApp"

    links { "engine" }

    targetdir ("build/tests/bin/build-%{cfg.architecture}")
    objdir    ("build/tests/obj/build-%{cfg.architecture}")

    files { "src/tests/runTests.cpp", "src/tests/**.cpp", "src/tests/**.hpp" }

    includedirs(inc)

    filter "configurations:Debug"
        defines { "DEBUG" }
        runtime "Debug"
        symbols "On"
    filter "configurations:Profile"
        defines { "NDEBUG" }
        runtime "Release"
        symbols "On"
        optimize "Speed"
        linktimeoptimization "On"
    filter "configurations:Dist"
        defines { "NDEBUG" }
        runtime "Release"
        symbols "Off"
        optimize "Speed"
        linktimeoptimization "On"
    filter {}
