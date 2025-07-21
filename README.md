## Render gallery: 
https://www.behance.net/alexandergolys


# Implemented features:
- shading with Phong material model, both vertex and face based triangular meshes
- hyperbolic geometry, conformal metric pullbacks
- Kleinian groups and hyperbolic tesselations
- Rendering circular triangles used for triangulations in hyperbolic plane \
- Pullbacking non-standard Riemannian metrics by conformal maps 
- mean, Gaussian and normal curvatures, principal directions, torsion, fundamental forms, shape operator 
- 1-parameter families of surfaces/curves realising animated deformation of leaves
- 2D rigid bodies and their basic dynamics
- 2D diffusion schemes
- Bezier and spline interpolations, Coons patches, tensor splines for surfaces
- Algebra of real functions on mesh/manifold and vector fields or differential forms as its modules
- Runge-Kutta ODE solver
- Second order partial differential equations solver 
- Integral transforms (Fourier, Laplace)
- Discrete transforms (DFT, generalised Fourier series, convolutions)
- Initial surface deformation encoded in PDE solution (or any non-homogeneous flow field)
- SDF based renderer 
- GLSL shader code generator realising SDF functions of primitives composed from C++ code level 
- Automatic shader code generating and manipulating code of function definitions allowing their realisation 
as elements of the function field on C++ level, modular shader code management and code macros substitutions 
- Multivariable functions allowing functorial composition of differential for elementary functions provided by user
- Rendering space curves as associated pipe surface with mesh generated on the GPU level so that curve can be deformed freely 
with no side effects in any moment within the code
- FFT in dimensions 1 and 2, fast convolutions of discrete functions, FFT-based rendering 
- SPH fluid simulations 
- Reading and processing audio signals, Gabor transforms, spectral power density 
- Foliated parametric surfaces
- Curve rendering based only on FFT-based differential operatorts
(no numerical differentiation, no differential geo, better local frames, curve may be already sampled
- config file storing paths to directories used implicitly in the code, so that no hardcoding is needed
- Time dependent flows 

# Next objectives:
- Minimal surfaces with Neumann boundary condition 
- Eulerian incompressible fluids
- rigid body physics in 3D
- hyperbolic geometry in 3D
- different non-Euclidean metrics
- connections
- discrete exterior calculus
- Hamiltonian formalism in physical dynamical system 
- Group actions, principal bundles
- Integral transforms (inverse Laplace, Fourier for general locally compact abelian group)
- Symbolic operations on elementary functions with tablebased values for certain transforms 
- Symbolic transform graph backtracking (e.g. symbolic differential operators can handle also something like sin(exp(t + 3))), not elementary itself
- Render to texture flow information
- Apply infinitesimal mesh transform in vertex shader 
- Add shadow maps


## Fixes and improvements TODO list:
- improve features available by design during animations 
- screenshots
- pause
- change speed
- rotate view
- better visual distinctions of directions and orientation of ambient space
- distinguish light positions with meshes
- add timer, and options for text base debugging areas
- Display numerical parameters used to control the motions in the scene
- Get some predefined environments of different general scene brightness and colors
- improve the code structure
