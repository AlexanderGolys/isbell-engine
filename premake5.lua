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

-- Start project: the selected scene (validated above)
startproject(selectedScene)

-- Include sub-projects
include "build/engine-build"

-- Include the validated scene project
include "build/scene-build"
