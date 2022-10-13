> First attempt at implemeting a basic Metropolis Scheme for Monte-Carlo simulation of a hard sphere system  
> Antoine Castagnede | October 2022

---

# WHAT DO WE NEED ?

- [x] Globally define a 2D box of given size --- Generalized to 3D
- Generate a suitable starting configuration --- start w/ SC
- [x] Read starting configuration from a file
- [x] Function to calculate "energy" of a given particle --- atm only check for overlap
- [x] Function to correctly write positions of particles to a file in the same fashion of Frank's visualization code --- be careful with extensions
- [x] Function to attempt to move a particle --- works as intended with a step size delta in [-0.5;0.5[]
- [x] Random number generation --- implemented using mt19937ar.c
- [x] Function to randomly select a particle --- already implemented in moveAttempt()
- What kind of small displacement to use ? --- depends on box size, radii, see Frenkel & Smit 3.2.2. Boundary Conditions
- How many cycles of moving attemps before writing the new system configuration ? --- chose arbitrarily, depends on number of particles and displacement size
- Implement the same using a cell list (makes it faster)
- [x] Implement periodic boundary conditions
- Implement nearest image convention for overlap computations

# NOTES :

- We want to restrict to 2D systems --- let z = 1 for all particles
  
# TYPESETTING FOR ALL ASSOCIATED .sph FILES :

```txt
&N
x-length y-length z-length
type x y z r
```

Where :

- N is the number of particles
- .-length is the size of the box in the . direction
- type, x, y, z, and r are as defined for the Disk structure below

And the last line is repeated N times to describe all particles independantly.

**N.B.** File extension matters, we use `.sph` to draw spheres.

# BUILD

To build, simply run `make` inside the directory.
