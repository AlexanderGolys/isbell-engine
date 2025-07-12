# Implemented features:
- shading with Phong material model, both vertex and face based triangular meshes
- hyperbolic geometry, conformal metric pullbacks
- Kleinian groups and hyperbolic tesselations, rendering circular triangles
- curvatures, torsion, fundamental forms, mean curvature flow
- deformations of 1-parameter families
- 2D rigid bodies
- 2D diffusion schemes
- Bezier and spline interpolations, Coons patches
- Function algebras, vector fields, differential forms
- custom curve rendering not causing deformations
- Runge-Kutta ODE solver
- Second order partial differential equations solver 
- Integral transforms (Fourier, Laplace)
- Discrete transforms (DFT, generalised Fourier series, convolutions)
- Initial surface deformation encoded in PDE solution (or any non-homogeneous flow field)


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
