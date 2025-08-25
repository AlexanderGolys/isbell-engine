-- One-at-a-time scene build
local base = _OPTIONS["scene"]
if not base or base == "" then
    error("Required option --scene=<BaseName> is missing.")
end

-- Simple Lua file existence check
local function file_exists(path)
    local f = io.open(path, "r")
    if f then f:close() return true end
    return false
end

local selectedFile = "../../src/render-projects/" .. base .. ".cpp"
if not file_exists(selectedFile) then
    error("Scene '" .. base .. "' not found at " .. selectedFile)
end

local projectName = base

project(projectName)
    -- Create per-scene project files under build/scene-build/<SceneName>
    location (projectName)
    kind "ConsoleApp"
    language "C++"
    cppdialect "C++26"

    targetdir ("../../bin/%{cfg.buildcfg}/" .. projectName)
    objdir ("../../bin-int/%{cfg.buildcfg}/" .. projectName)

    includedirs {
        "../../src"
    }

    dependson { "engine" }
    links { "engine" }

    -- Only the selected .cpp; allow IDE visibility of .hpp
    files { selectedFile, "../../src/render-projects/**.hpp" }

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
