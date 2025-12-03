# isbell : rendering engine for real-time animations visualising mathematical concepts




## Gallery of renders done with Isbell
https://www.behance.net/alexandergolys

## Changes in Isbell 0.2.0

1. Rendering sessions now configure from a single hub. Tweaking window behaviour or frame pacing no longer requires rummaging through old scripts.
2. Window, keyboard, mouse, and UI input moves through one queue. Bursts of events stay in order instead of disappearing during frame spikes.
3. Rendering, updating, and shutdown each live in their own stage. Start-up quirks no longer bleed into the main loop.
4. Layer components mix and match like blocks. Complex behaviour combinations no longer demand glue code.
5. Scene assembly happens in one place. Cameras, lights, and meshes drop in together instead of being wired manually.
6. GPU data flows through one toolkit. Streaming dense geometry or dynamic data feels routine rather than risky.
7. OpenGL state changes ride through safety rails. Rendering code stays focused while helper utilities catch bad bindings.
8. Shader setups spin up with a quick call. Experimenting with new looks takes minutes rather than afternoons.
9. The interaction layer recognises focused listener roles. Features only hear the inputs they actually need.
10. Control hooks now accept any time-driven behaviour. Animations and simulations plug in without scaffolding.
11. Camera navigation is ready out of the box. Orbiting, strafing, zooming, and scroll tweaks feel smooth immediately.
12. Impulse-driven motion is standard. Bursts of energy or easing curves appear without hand-tuned maths.
13. Cameras keep projection and view data aligned. Even complex moves stay free from jitter and clipping.
14. Materials carry their textures and light responses together. Look development runs faster and invites experimentation.
15. Lighting data follows one format. Scaling from a single lamp to a multi-light scene becomes a checkbox task.
16. The maths toolkit supports richer transformations. New simulations or coordinate tricks are far less of a chore.
17. Interactive overlays now recognise hover and click states. Heads-up displays and control panels are finally within reach.
18. Window presets cover popular resolutions and sticky-input variants. Getting a new demo on screen is a one-minute job.
19. Logging respects significance levels. Deep dives stay available while everyday runs stay quiet.
20. A flagship demo showcases the full workflow. Older experiments remain archived for comparison.



## Implemented features:
- using vertex, fragment, geometry and compute shaders with OpenGL API
- shading with Phong material model, both vertex and face based triangular meshes
- hyperbolic geometry, conformal metric pullbacks
- Kleinian groups and hyperbolic tesselations
- Rendering triangulated surfaces embedded in hyperbolic plane with circular triangles mode (works fine for different conformal models as well, geodesics are not really required to be circular, convexity/concavity seemed to be a sufficient conditions but not always necessary)
- Pullbacking Riemannian metrics by conformal maps 
- mean, Gaussian and normal curvatures, principal directions, torsion, fundamental forms, shape operator 
- 1-parameter families of surfaces/curves realising animated deformation of leaves
- Rendering space curves as associated pipe surface with mesh generated on the GPU level guaranteeing distortion-free deformations
- FFT-based curve renderer (framework works after if the curve is reduced to finite sequence of points, uses only FFT-based differentiation and no differential geometry)
- 2D rigid bodies and their basic dynamics
- Reading and processing audio signals, Gabor transforms, spectral power density
- Bezier and spline interpolations, Coons patches, tensor splines for surfaces
- Algebra of real functions on mesh/manifold and vector fields or differential forms as its modules
- Runge-Kutta ODE solver
- Second order partial differential equations solver on &#9633 
- Integral transforms (Fourier, Laplace)
- Discrete transforms (DFT, generalised Fourier series, convolutions)
- FFT in dimensions 1 and 2, fast convolutions of discrete functions, FFT-based rendering
- Initial surface deformation encoded in PDE solution (or any non-homogeneous flow field)
- SDF based renderer 
- GLSL shader code generator realising SDF functions of primitives composed from C++ code level 
- Automatic shader code generating and manipulating code of function definitions 
- Modular shader code management and code macros substitutions 
- SPH fluid simulations 
- Foliated parametric surfaces
- Differential operators of sampled data with FFT based implementation
- Time dependent flows 
- basic windows filesystem operations and file buffering 




Isbell can be (and has been successfully) used for rendering different things as well, although it may lack certaing features or design patterns used and expected from game engines, as instead I've 
focused more on developping interfaces and ways of constructing and manipulating various types of mathematical objects, such as manifolds, hyperbolic spaces, integral transforms, differential equations, surface deformations etc. 
I also implemented some elements of physics simulation, for example 2D rigid bodies and SPH fluid sims, although the more comprehensive physics engine is one of the goals that I plan to implement as one of the next major updates.


## Build instructions:

For the build system works on Windows only, but all* the tools used in the project are cross-platform, I just didn't test it on other systems but it can be done if needed in the future with with no substancial changes 

- choose the main file with selected animation. It should be placed in src/render-projects and be a cpp file.
- run premake `.\premake\premake5 --scene=SCENE_NAME_FILE cmake` (for SCENE_NAME_FILE using the name without .cpp suffix)
- alternatively, you can do this by running the script file `rebuild.bat`, editting the scene name in it if needed
- after premake generates the CMakeLists.txt file, load it as usual CMake project in your IDE or with any other way you use to do that if you are the notepad type of guy (although C++ is probably the last language you want to look at in this case, so it is rather extremely unlikely you'd be reading this file, but if I'm wrong hit me up as I'm curious about your story)

<!-- TODO: update this readme to something nice looking  -->
