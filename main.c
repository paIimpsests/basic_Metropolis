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
 *      moveAttempt()
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

#define MAX_SIZE 1000
#define X_LENGTH 11.0f
#define Y_LENGTH 11.0f
#define Z_LENGTH 11.0f

typedef struct Sphere Sphere;
struct Sphere
{
	double x;	// x position of the center of the sphere
	double y;	// y position of the center of the sphere
        double z;       // z position of the center of the sphere
	double r;	// radius of the sphere
        char type;      // particle type --- for visualization purposes only

};

int overlapCheck(Sphere *sphere1, Sphere *sphere2)
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
	if ( ( ( dx - X_LENGTH * rint(dx / X_LENGTH) )
               * ( dx - X_LENGTH * rint(dx / X_LENGTH) )
               + ( dy - Y_LENGTH * rint(dy / Y_LENGTH) )
               * ( dy - Y_LENGTH * rint(dy / Y_LENGTH) )
               + ( dz - Z_LENGTH * rint(dz / Z_LENGTH) )
               * ( dz - Z_LENGTH * rint(dz / Z_LENGTH) )
             )
           < ( (sphere1->r + sphere2->r) * (sphere1->r + sphere2->r) )
           )
	        overlap = 1;
	return overlap;
}

int overlapCheckGlobal(Sphere **spheres, Sphere *sample, int n)
{
        /* Function:    overlapCheckGlobal
         * -------------------------------
         *  Checks for overlap between one given particle and the rest of the
         *  system
         *
         *  **spheres:  pointer to a pointer to the data of the system
         *  *sample:    pointer to the data to check overlap for
         *  n:          location of the particle that we want to check overlap
         *              for, 0 <= n <= MAX_SIZE - 1
         *
         *  return:     1 if there is overlap, 0 if not
         */
        int overlap = 0, i = 0;
        while ((overlap == 0) && (i < MAX_SIZE))
        {
                overlap = (i != n) * overlapCheck(sample, spheres[i]);
                i++;
        }
        return overlap;
}

void readInit(char *filename, Sphere **spheres)
{
        /*
         * Function:    readInit
         * ---------------------
         * Initializes the table of data of all MAX_SIZE particles by reading
         * it from a user supplied .txt file
         * NB1: Some parameters (total number of particles and box size) from
         * the initialization .sph file must match those defined here
         * NB2: Reading is done only for files with typesetting indicated above
         *
         * *filename:   pointer to the name of the initialization .txt file
         * **spheres:   pointer to the pointer to the table of particles data
         *
         */
        FILE *initfile = NULL;
        int nParticles = 0, i = 0;
        double xlength = 0.0f, ylength = 0.0f, zlength = 0.0f;
        initfile = fopen(filename, "r");
        if (initfile != NULL)
        {
                fscanf(initfile, "%*c%d%*c", &nParticles);
                fscanf(initfile, "%lf %lf %lf%*c", &xlength, &ylength, &zlength);
                if ((nParticles == MAX_SIZE) && (xlength == X_LENGTH) && (ylength == Y_LENGTH) && (zlength == Z_LENGTH))
                {
                        for (i = 0; i < MAX_SIZE; i++)
                        {
                             fscanf(initfile, "%c %lf %lf %lf %lf%*c", &spheres[i]->type, &spheres[i]->x, &spheres[i]->y, &spheres[i]->z, &spheres[i]->r);
                        }

                }
                else
                        printf("Initialization did not happen, check the initialization source file.\n"); 
                fclose(initfile);
        }
}

void writeCoords(char *filename, Sphere **spheres)
{
        /*
         * Function:    writeCoords
         * -------------------------
         * Writes the position, radius, and type data for all MAX_SIZE
         * particles in a .sph file
         * NB: Writing is done using the typesetting indicated above
         *
         * *filename:   pointer to the name of the initialization .sph file
         * **spheres:   pointer to the pointer to the table of particles data
         */
        FILE *outfile = NULL;
        int i = 0;
        outfile = fopen(filename, "a");
        if (outfile != NULL)
        {
                fprintf(outfile, "&%d\n", MAX_SIZE);
                fprintf(outfile, "%lf %lf %lf\n", X_LENGTH, Y_LENGTH, Z_LENGTH);
                for (i = 0; i < MAX_SIZE; i++)
                        fprintf(outfile, "%c %lf %lf %lf %lf\n", spheres[i]->type, spheres[i]->x, spheres[i]->y, spheres[i]->z, spheres[i]->r);
                fclose(outfile);
        }
}

void moveAttempt(Sphere **spheres)
{
        /*
         * Function: moveAttempt                
         * ---------------------
         * Tries to move a randomly selected particle in space by a small amount
         * and checks for overlap
         * If there is no overlap, overwrites the previous particle position
         * with the new one
         *
         * **sphere:    pointer to the pointer to the table of particles data
         */ 
        double delta[3] = {(genrand() - 0.5), (genrand() - 0.5), (genrand() - 0.5)};
        int n = (int) (genrand() * MAX_SIZE);
        Sphere  bufferSphere = {        fmod((spheres[n]->x + delta[0]) + 2 * X_LENGTH, X_LENGTH),
                                        fmod((spheres[n]->y + delta[1]) + 2 * Y_LENGTH, Y_LENGTH),
                                        fmod((spheres[n]->z + delta[2]) + 2 * Z_LENGTH, Z_LENGTH),
                                        spheres[n]->r,
                                        spheres[n]->type
                                },
                *pBufferSphere = NULL;
        pBufferSphere = &bufferSphere;
        if (!overlapCheckGlobal(spheres, pBufferSphere, n))
        {
                spheres[n]->x = bufferSphere.x;
                spheres[n]->y = bufferSphere.y;
                spheres[n]->z = bufferSphere.z;
        }
}

void initSCConfig(Sphere **spheres)
{
        /*
         * Function:    initSCConfig
         * -------------------------
         * Generates a simple cubic starting configuratin based on globally
         * defined parameters such as box size, number of particles and
         * particles radius
         * NB: For now we assume the use of a single fixed radius for all
         * particles
         *
         * **sheres:    pointer to the pointer to the table of particles data
         *
         */

        int i = 0, j = 1, k = 1, l = 1;
        double r = 0.5f, x = 0.0f, y = 0.0f, z = 0.0f;
        double sizeup = ceil(pow(MAX_SIZE, 1.0f/3.0f));
        if ( !(2 * r * (sizeup + 1) < X_LENGTH)
             && !(2 * r * (sizeup + 1) < Y_LENGTH)
             && !(2 * r * (sizeup + 1) < Z_LENGTH)
           )    // Checks if all particles can fit inside the box
        {
                while (i < MAX_SIZE)
                {
                        for (j = 1; j < (int) sizeup; j++)
                        {
                                for (k = 1; k < (int) sizeup; k++)
                                {
                                        for (l = 1; l < (int) sizeup; l++)
                                        {
                                                *spheres[i] = (Sphere) {r * (1 + 2 * j),
                                                r * (1 + 2 * k),
                                                r * (1 + 2 * l),
                                                r,
                                                'a'};
                                                i++;
                                        }
                                }
                        }
                }

        }
}

void tuneStepSize(double *ptune, double *pAcceptanceRate)
{
        /*
         * Function:    tuneStepSize
          * -------------------------
         * Tunes up the step size employed in trying to displace a particle by
         * 5% if the acceptance rate is < 30% and down by 5% if the latter is
         * < 30%
         *
         * *ptune:              pointer to the value of the tuning factor
         * *acceptanceRate:     pointer to the value of the acceptance rate for
         *                      moving a particle in the system
         *
         */
        if (*pAcceptanceRate < 0.3)
                *ptune -= 0.05;
        else
                *ptune += 0.05;
}

int main(int argc, char *argv[])
{
        init_genrand(920188);
        int i = 0, n = 0;
        double  tune = 0.0f,
                *pTune = &tune,
                acceptanceRate = 0.0f,
                *pAcceptanceRate = &acceptanceRate;
        Sphere *spheres = malloc(MAX_SIZE * sizeof(Sphere));
        Sphere **pSpheres = malloc(MAX_SIZE * sizeof(Sphere));
        for (int i = 0; i < MAX_SIZE; ++i)
                pSpheres[i] = spheres+i;


        // Initialization of the spheres coordinates one by one
        /*
        spheres[0] = (Sphere) {0.0f, 0.0f, 1.0f};
        spheres[1] = (Sphere) {0.0f, 1.9f, 1.0f};
        for (int i = 0; i < MAX_SIZE; i++)
               spheres[i] = (Sphere) {0.0f, 0.0f, 1.0f};
        */

        // Checks for moveAttempt()
        readInit("SC_init.sph", pSpheres);
        writeCoords("coords.sph", pSpheres);
        for (i = 0; i < 10000000; i++)
        {
                moveAttempt(pSpheres);        
                if (i % 10000 == 0)
                        writeCoords("coords.sph", pSpheres);
        }

        free(spheres);
        free(pSpheres);
        return 0;
}
