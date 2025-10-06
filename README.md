# isbell : rendering engine for real-time math visualisations 


Isbell can be (and has been successfully) used for rendering different things as well, although it may lack certaing features or design patterns used and expected from game engines, as instead I've 
focused more on developping interfaces and ways of constructing and manipulating various types of mathematical objects, such as manifolds, hyperbolic spaces, integral transforms, differential equations, surface deformations etc. 
I also implemented some elements of physics simulation, for example 2D rigid bodies and SPH fluid sims, although the more comprehensive physics engine is one of the goals that I plan to implement as one of the next major updates.


## Gallery of renders done with Isbell
https://www.behance.net/alexandergolys

## Build instructions:

For the build system works on Windows only, but all* the tools used in the project are cross-platform, I just didn't test it on other systems but it can be done if needed in the future with with no substancial changes 

- choose the main file with selected animation. It should be placed in src/render-projects and be a cpp file.
- run premake `.\premake\premake5 --scene=SCENE_NAME_FILE cmake` (for SCENE_NAME_FILE using the name without .cpp suffix)
- alternatively, you can do this by running the script file `rebuild.bat`, editting the scene name in it if needed
- after premake generates the CMakeLists.txt file, load it as usual CMake project in your IDE or with any other way you use to do that if you are the notepad type of guy (although C++ is probably the last language you want to look at in this case, so it is rather extremely unlikely you'd be reading this file, but if I'm wrong hit me up as I'm curious about your story)



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



*This is true besides two windows specific methods responsible for creating handle, mapping and flushing the file, although as I've included <windows.h> only in their cpp file and they are not using by different parts of the of isbell-engine (as for now), 
at at least in theory they should not cause problems with cross-platform compilation. 
Thus the interface for compute shaders using SSBO may require user to map a file to work as intended and tested.
These are also tested in unittests, so should be removed from tests or overriden on other systems if needed to run tests (although there is little chance anywhere will need to do that on the single-dev project besides the one guy that is writing this instructions, so probably is aware of these complications after the windows.h macro hell it caused after including this innocent little header that is overriding half of my namespace such as min and max).
