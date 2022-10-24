/* First attempt at implemeting a basic Metropolis Scheme for Monte-Carlo
 * simulation of a hard sphere system
 * Antoine Castagnede | October 2022
 *
 * WHAT DO WE NEED ?
 * -----------------
 * (o)  Globally define a 2D box of given size --- Generalized to 3D
 * (o)  Generate a suitable starting configuration --- Implemented for SC,
 *      will remove at a later time 
 * (o)  Read starting configuration from a file
 * (o)  Function to calculate "energy" of a given particle --- atm only check
 *      for overlap
 * (o)  Function to correctly write positions of particles to a file in the
 *      fashion of Frank's visualization code --- be careful with extensions
 * (o)  Function to attempt to move a particle --- works as intended with
 *      delta in [-0.5 ; 0.5[
 * (o)  Random number generation --- removed previous implementation, now uses
 *      mt19937ar.c
 * (o)  Function to check for overlap with the whole system
 * (o)  Function to randomly select a particle --- already implemented in
 *      particleMove()
 *  -   What kind of small displacement to use ? --- depends on system density,
 *      see Frenkel & Smit 3.1.1. Translational moves
 *  -   How many cycles of moving attemps before writing the new system
 *      configuration ? --- chose arbitrarily, depends on number of particles
 *      and displacement size
 *  -   Implement the same using a cell list (makes it faster)
 * (o)  Implement PBC
 * (o)  Implement nearest image convention for overlap computations
 *
 *  TYPESETTING FOR ALL ASSOCIATED .SPH FILES :
 *  -------------------------------------------
 *  """
 *  &N
 *  x-length y-length z-length
 *  type x y z r
 *  """
 *  Where :
 *  - N is the number of particles
 *  - .-length is the size of the box in the . direction
 *  - type, x, y, z, and r are as defined for the Disk structure below
 *  And the last line is repeated N times to describe all particles
 *  independantly
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>       // need to compile with `-lm` flag
#include "mt19937ar.c"

#define T 300           // absolute temperature
#define kB 1.380649e-23 // Boltzmann constant
#define eps kB*T        // unit of energy


typedef struct Sphere Sphere;
struct Sphere
{
	double x;	// reduced x coordinate 
	double y;	// reduced y coordinate
        double z;       // reduced z coordinate
	double r;	// reduced radius --- redundant with sigma ... might remove later on
        char type;      // particle type --- for visualization purposes only

};

int overlapCheck(Sphere *sphere1, Sphere *sphere2, double L)
{
        /*
         * Function:    overlapCheck
         * -------------------------
         * Checks for overlap between two given particules by computing the
         * distance between their centers and comparing it to the sum of their
         * radii
         *
         * *sphere1:      pointer to the data relative to the first particle
         * *sphere2:      pointer to the data relative to the second particle
         *
         * return:      1 is there is overlap, 0 if not
         */

	int overlap = 0;
        double  dx = sphere1->x - sphere2->x,
                dy = sphere1->y - sphere2->y,
                dz = sphere1->z - sphere2->z;
	if ( ( ( dx - L * rint(dx / L) )
               * ( dx - L * rint(dx / L) )
               + ( dy - L * rint(dy / L) )
               * ( dy - L * rint(dy / L) )
               + ( dz - L * rint(dz / L) )
               * ( dz - L * rint(dz / L) )
             )
           < ( (sphere1->r + sphere2->r) * (sphere1->r + sphere2->r) )
           )
	        overlap = 1;
	return overlap;
}


int overlapCheckGlobal(Sphere **spheres, Sphere *sample, int n, int N, double L)
{
        /* Function:    overlapCheckGlobal
         * -------------------------------
         *  Checks for overlap between one given particle and the rest of the
         *  system
         *
         *  **spheres:  pointer to a pointer to the data of the system
         *  *sample:    pointer to the data to check overlap for
         *  n:          location of the particle that we want to check overlap
         *              for, 0 <= n <= N - 1
         *  N:          number of particles
         *  L:          reduced box size
         *
         *  return:     1 if there is overlap, 0 if not
         */
        int overlap = 0, i = 0;
        while ((overlap == 0) && (i < N))
        {
                overlap = (i != n) * overlapCheck(sample, spheres[i], L);
                i++;
        }
        return overlap;
}

void readInit(  char *filename,
                Sphere **spheres,
                Sphere ***pSpheres,
                int *pN,
                double *pL,
                double *pSigma
             )
{
        /*
         * Function:    readInit
         * ---------------------
         * Initializes the table of data of all N particles by reading
         * it from a user supplied .sph file; also retrieves unit of length for
         * the system from particles **unique** radius; particles coordinates
         * are written in reduced units; also retrieves **cubic** box size
         * NB: Reading is done only for files with typesetting indicated above
         * NB: It is assumed that the supplied init file is written in standard units so
         * conversion to reduced units is applied
         *
         * *filename:   pointer to the name of the initialization .txt file
         * **spheres:   pointer to the pointer to the table of particles data
         * ***pSpheres: pointer to the pointer to the pointer to the table of
         *              particles data
         * *pN:         pointer to the number of particles in the system
         * *pL:         pointer to the reduced box size
         * *pSigma:     pointer to the value of the system's unit of length
         *
         */

        FILE *initfile = NULL;
        int i = 0;
        double boxSize = 0.0f;
        double buffer[4] = {0.0f};
        initfile = fopen(filename, "r");
        if (initfile != NULL)
        {
                // Read value of N and allocate memory accordingly
                fscanf(initfile, "%*c%d%*c", pN);
                *spheres = malloc(*pN * sizeof(Sphere));
                *pSpheres = malloc(*pN * sizeof(Sphere));
                for (int i = 0; i < *pN; i++)
                        (*pSpheres)[i] = (*spheres) + i;
                // Read value of box size
                fscanf(initfile, "%lf %*f %*f%*c", &boxSize);
                // Read value of radius to obtain sigma
                fscanf(initfile, "%*c %*f %*f %*f %lf%*c", &buffer[0]);
                *pSigma = 2 * buffer[0];
                *pL = boxSize / *pSigma;
                // Reset cursor at the begining of file and skip through
                // redundant lines
                fseek(initfile, 0, SEEK_SET);
                fscanf(initfile, "%*c%*d%*c");
                fscanf(initfile, "%*f %*f %*f%*c");
                // Populate table of particles data in reduced units
                for (i = 0; i < *pN; i++)
                {
                        fscanf  ( initfile,
                                  "%c %lf %lf %lf %lf%*c",
                                  &(*pSpheres)[i]->type,
                                  &buffer[0],
                                  &buffer[1],
                                  &buffer[2],
                                  &buffer[3]
                                );
                        (*pSpheres)[i]->x = buffer[0] / *pSigma;
                        (*pSpheres)[i]->y = buffer[1] / *pSigma;
                        (*pSpheres)[i]->z = buffer[2] / *pSigma;
                        (*pSpheres)[i]->r = buffer[3] / *pSigma;

                }
                fclose(initfile);
        }
}

void writeCoords(char *filename,
                 Sphere **pSpheres,
                 int *pN,
                 double *pL,
                 double *pSigma
                )
{
        /*
         * Function:    writeCoords
         * -------------------------
         * Writes the position, radius, and type data for all N
         * particles in a .sph file
         * NB: Writing is done using the typesetting indicated above
         *
         * *filename:   pointer to the name of the initialization .sph file
         * **pSpheres:  pointer to the pointer to the table of particles data
         * *pN:         pointer to the number of particles in the system
         * *pL:         pointer to the reduced box size
         * *pSigma:     pointer to the value of the system's unit of length
         *
         */
        FILE *outfile = NULL;
        int i = 0;
        double boxSize = *pL * *pSigma;
        outfile = fopen(filename, "a");
        if (outfile != NULL)
        {
                fprintf(outfile, "&%d\n", *pN);
                fprintf(outfile, "%lf %lf %lf\n", boxSize, boxSize, boxSize);
                for (i = 0; i < *pN; i++)
                        fprintf(outfile,
                                "%c %lf %lf %lf %lf\n",
                                pSpheres[i]->type,
                                pSpheres[i]->x * *pSigma,
                                pSpheres[i]->y * *pSigma,
                                pSpheres[i]->z * *pSigma,
                                pSpheres[i]->r * *pSigma
                               );
                fclose(outfile);
        }
}

void particleMove(Sphere **pSpheres, int N, double L, double sigma, int *pSuccessP, double particleStepTune)
{
        /*
         * Function: particleMove                
         * ---------------------
         * Tries to move a randomly selected particle in space by a small amount
         * and checks for overlap
         * If there is no overlap, overwrites the previous particle position
         * with the new one
         *
         * **pSphere:   pointer to the pointer to the table of particles data
         * N:           number of particles
         * L:           reduced box size
         * sigma:       system unit length
         * *pSuccessP:  pointer to a value accounting whether particle move
         *              attempt was successful (+=1) or not (+=0)
         * particleStepTune:    tuning factor for particle step size; such that
         *                      acceptance rate is ~30%
         *
         */ 
        double delta[3] = {(genrand() - 0.5) / sigma * (1 + particleStepTune),
                           (genrand() - 0.5) / sigma * (1 + particleStepTune),
                           (genrand() - 0.5) / sigma * (1 + particleStepTune)
                          };
        int n = (int) (genrand() * N);
        Sphere  bufferSphere = {        fmod((pSpheres[n]->x + delta[0]) + 2 * L, L),
                                        fmod((pSpheres[n]->y + delta[1]) + 2 * L, L),
                                        fmod((pSpheres[n]->z + delta[2]) + 2 * L, L),
                                        pSpheres[n]->r,
                                        pSpheres[n]->type
                                },
                *pBufferSphere = NULL;
        pBufferSphere = &bufferSphere;
        if (!overlapCheckGlobal(pSpheres, pBufferSphere, n, N, L))
        {
                pSpheres[n]->x = bufferSphere.x;
                pSpheres[n]->y = bufferSphere.y;
                pSpheres[n]->z = bufferSphere.z;
                *pSuccessP += 1;
                
        }
}

void volumeMove(Sphere **pSpheres, double *pL, int N, double sigma, double P, int *pSuccessV, double volumeStepTune)
{
        /*
         * Function:    volumeMove
         * -----------------------
         * Attempts to change the volume of the simulation box by a small
         * amount according to the known acceptance rule
         * 
         * Notes: implementation of this function requires to revamp how the
         * box size is implemanted --- this is WIP and being gradually
         * implemented alongside the use of reduced units
         * A way to tune the step size of the volume change also needs to be
         * implemented --- this must be done by passing forward a parameter that
         * accounts for the success rate of the present function
         *
         * **pshperes:  pointer to the pointer to the table of particles data
         * L:           reduced box size
         * N:           number of particles
         * sigma:       system unit length
         * P:           reduced pressure
         * *pSuccess:   pointer to a value indicating wheter attempt to move
         *              volume was successful or not (1 is y, 0 is n)
         * *pSoE:       pointer to a shrinking/expanding indicator (1 is
         *              expanding, 0 is shrinking)
         * *pSuccessV:  pointer to a value accounting whether particle move
         *              attempt was successful (+=1) or not (+=0)
         * volumeStepTune:      tuning factor for volume step size; such that
         *                      acceptance rate is ~30%
         *
         */
        double  vol = *pL * *pL * *pL,
                delta[2] = {(genrand() - 0.5) / (800 * sigma) * (1 + volumeStepTune), genrand()},
                newVol = vol * (1 + delta[0]),                  // vol changes are in
                newL = pow(newVol, 1.0f/3.0f),                  // the order of +-5%
                ratio = newL / *pL;
        Sphere *newSpheres = malloc(N * sizeof(Sphere));
        Sphere **pNewSpheres = malloc(N * sizeof(Sphere));
        for (int i = 0; i < N; i++)
        {
                pNewSpheres[i] = newSpheres+i;
                newSpheres[i] = (Sphere) {      ratio * pSpheres[i]->x,
                                                ratio * pSpheres[i]->y,
                                                ratio * pSpheres[i]->z,
                                                pSpheres[i]->r,
                                                pSpheres[i]->type
                                         };
        }
        if (!overlapCheckGlobal(pNewSpheres, pNewSpheres[0], 0, N, newL))
        {
                double rule = pow(newVol / vol, N) * exp(- P * (newVol - vol));
                if (delta[1] < rule)
                {
                        for (int i = 0; i < N; i++)
                        {
                                pSpheres[i]->x = newSpheres[i].x;
                                pSpheres[i]->y = newSpheres[i].y;
                                pSpheres[i]->z = newSpheres[i].z;
                        }
                        *pL = newL;
                        *pSuccessV += 1;
                }
        }
        free(newSpheres);
        free(pNewSpheres);
}

void tuneStepSize(double *stepTuner, int *pNSuccess, double acceptanceRate, int nCycles)
{
        /*
         * Function:    tuneStepSize
         * -------------------------
         * Tunes the step size for volume or particle moves by +-2% depending
         * on whether the success rate is above or below the specified
         * acceptance rate
         *
         * *stepTuner:  pointer to the step size tuner
         * *pNSuccess:  pointer to the number of successes over the last nCycles
         *              cycles
         * acceptanceRate:      targeted acceptance rate
         * nCycles:     number of cycles to average the number of successes over
         *
         */
                if (pNSuccess - (int) (acceptanceRate * nCycles) > 0)
                                *stepTuner += 0.02;
                        else
                                *stepTuner -= 0.02;
                        *pNSuccess = 0;
       
}

double measurePF(int N, double sigma, double L)
{
        /*
         * Function: measurePF
         * -------------------
         *  Computes the reduced packing fraction for the given parameters
         *
         *  N:          number of particles in the system
         *  sigma:      particle size
         *  L:          reduced box size
         *
         *  return:     reduced packing fraction
         *
         */
        double PF = (double) N * sigma * sigma * sigma * M_PI / (6.0f * L * L * L);
        return PF;
}


int main(int argc, char *argv[])
{
        // System parameters
        int N = 0,                      // number of particles
            *pN = &N;
        double L = 0.0f,                // reduced box size
               *pL = &L;
        double P = 5.0f,                // reduced pressure
               *pP = &P;
        double sigma = 0.0f,            // particle size = length unit
               *pSigma = &sigma;
        double PF = 0.0f;               // system packing fraction
        double rho = 0.0f;              // reduced number density
                                        //
        // MC parameters
        init_genrand(486488);
        double proba = 0.0f;
        int i = 0, Ncycles = 10000000;
        Sphere *spheres = NULL;
        Sphere **pSpheres = NULL;
        int successP = 0,
            *pSuccessP = &successP;
        double particleStepTune = 0.0f,
               *pParticleStepTune = &particleStepTune;
        int successV = 0,
            *pSuccessV = &successV;
        double volumeStepTune = 0.0f,
               *pVolumeStepTune = &volumeStepTune;
        FILE *writefile = NULL;

        readInit("SC_init_6.sph", &spheres, &pSpheres, pN, pL, pSigma);
        writeCoords("coords.sph", pSpheres, pN, pL, pSigma);

        // 1st half of simuation time
        // Tuning of step sizes to obtain ~30% acceptance rate
        for (i = 0; i < Ncycles / 2; i++)
        {
                proba = genrand();
                if (proba < (1.0f / (float) (N + 1)))
                        volumeMove(pSpheres, pL, N, sigma, P, pSuccessV, volumeStepTune);
                else
                        particleMove(pSpheres, N, L, sigma, pSuccessP, particleStepTune);
                if (i % 10000 == 0)
                        tuneStepSize(pParticleStepTune, pSuccessP, 0.30f, 10000);
                if (i % ((N + 1) * 10000) == 0)
                        tuneStepSize(pVolumeStepTune, pSuccessV, 0.30f, (N + 1) * 10000);
        }

        // 2nd half of simulation time
        // Measurements of observables --- here: packing fraction to obtain number density
        for (i = 0; i < Ncycles / 2; i++)
        {
                proba = genrand();
                if (proba < (1.0f / (float) (N + 1)))
                        volumeMove(pSpheres, pL, N, sigma, P, pSuccessV, volumeStepTune);
                else
                        particleMove(pSpheres, N, L, sigma, pSuccessP, particleStepTune);
                PF += measurePF(N, sigma, L);
                if (i % 10000 == 0)
                {
                        writeCoords("coords.sph", pSpheres, pN, pL, pSigma);
                        PF = PF / 10000;
                        rho = 6.0f * PF / (M_PI * sigma * sigma * sigma);
                        writefile = fopen("density.txt", "a");
                        fprintf(writefile, "%d\t%lf\n", i, rho);
                        fclose(writefile);
                        PF = 0.0f;

                }
        } 

        free(spheres);
        free(pSpheres);
        return 0;
}
