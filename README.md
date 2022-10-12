> First attempt at implemeting a basic Metropolis Scheme for Monte-Carlo simulation of a hard sphere system  
> Antoine Castagnede | October 2022

---

# WHAT DO WE NEED ?

- [x] Globally define a 2D box of given size --- Generalized to 3D
- Generate a suitable starting configuration --- ask Frank about that
- [x] Read starting configuration from a file
- [x] Function to calculate "energy" of a given particle --- atm only check for overlap
- [x] Function to correctly write positions of particles to a file in the same fashion of Frank's visualization code
- Function to attempt to move a particle --- must call writeCoords(), changeCoords(), overlapCheck()
- [x] Random number generation --- meh... ask Frank about that it properly
- Function to randomly select a particle
- What kind of small displacement to use ?
- How many cycles of moving attemps before writing the new system configuration ? --- chose arbitrarily, depends on number of particles and displacement size

# NOTES :

- We want to restrict to 2D systems --- let z = 1 for all particles
  
# TYPESETTING FOR ALL ASSOCIATED .TXT FILES :

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

# BUILD

To build, simply run `make` inside the directory.
