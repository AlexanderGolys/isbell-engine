workspace "Renderer"
   configurations {"Debug", "Release", "Dist",  "Profile"}
   architecture "x64"
   location "build"

newoption {
   trigger = "scene",
   value = "SCENE_NAME",
   description = "Base name of a scene main file in src/render-projects (without extension)"
}

-- Pick a sensible start project only if the requested scene .cpp exists
local selectedScene = _OPTIONS["scene"]

-- Make --scene mandatory for generation
if not selectedScene or selectedScene == "" then
   error("Required option --scene=<BaseName> is missing. Choose a file base name from src/render-projects (without .cpp).")
end

-- Simple Lua file existence check (avoid tool warnings on os.isfile)
local function file_exists(path)
   local f = io.open(path, "r")
   if f then f:close() return true end
   return false
end

local function scene_exists(name)
   if not name then return false end
   return file_exists("src/render-projects/" .. name .. ".cpp")
end

-- Fail if the requested scene doesn't exist
local scene_ok = scene_exists(selectedScene)
if not scene_ok then
   error("Requested scene '" .. selectedScene .. "' not found at src/render-projects/" .. selectedScene .. ".cpp")
end

-- Start project: the selected scene (validated above)
startproject(selectedScene)

-- Include sub-projects
include "build/engine-build"

-- Include the validated scene project
include "build/scene-build"
