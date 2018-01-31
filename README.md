# IF-DDA

IF-DDA is a numerical techniques for solving the electromagnetic
scattering problem in three dimensions. IF-DDA is based on the DDA
(discrete dipole approximation) which is a volume-integral equation
method.  The DDA (also referred to as the coupled dipole method) was
originally proposed by by Purcell and Penny packer where the object
under study is discretized into a set of small subunits and the field
at each subunit position is computed through a self consistent
equation. Then the diffracted field can be computed easily.

This method can be used to arbitrarily shaped, inhomogeneous,
anisotropic particles. The radiation condition is automatically
satisfied, because the Green's function satisfies the radiation
condition.  The computation is confined to the volume of the
scatterer, hence this method does not need any PML (perfect matching
layer).

IF-DDA has a friendly guide user interface where many particles
(cuboid, sphere, ellipsoid, many spheres,...), beams (plane wave,
Gaussian wave, multiple plane waves,...) are accessibles with a
drop-down menu. The studies are selected with the mouse: 
-Cross section 
-Poynting vector 
-Microscopy 
-Optical force 
-Optical torque
-Near field

There is two branches for the code. One with FFTE (fast Fourier transform
in the east) and one in parallel (openmp) with FFTW (fast Fourier
transform in the west). For the version with openmp you should know
how to play with ulimit and fix the number of processor...
