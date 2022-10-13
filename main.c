/* First attempt at implemeting a basic Metropolis Scheme for Monte-Carlo
 * simulation of a hard sphere system
 * Antoine Castagnede | October 2022
 *
 * WHAT DO WE NEED ?
 * -----------------
 * (o)  Globally define a 2D box of given size --- Generalized to 3D
 *  -   Generate a suitable starting configuration --- start w/ SC
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
 *  NOTES :
 *  -------
 *  - We want to restrict to 2D systems --- let z = 1 for all particles 
 *
 *  TYPESETTING FOR ALL ASSOCIATED .TXT FILES :
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

#define MAX_SIZE 2
#define X_LENGTH 5
#define Y_LENGTH 5
#define Z_LENGTH 1

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
	if ( ( (sphere1->x - sphere2->x) * (sphere1->x - sphere2->x)
               + (sphere1->y - sphere2->y) * (sphere1->y - sphere2->y)
               + (sphere1->z - sphere2->z) * (sphere1->z - sphere2->z)
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
         *              for, 0 <= n <= 99
         *
         *  return:     1 if there is overlap, 0 if not
         */
        int overlap = 0, i = 0;
        while ((overlap == 0) && (i < MAX_SIZE - 1))
        {
                overlap = (i != n) * overlapCheck(sample, spheres[i]);
                i++;
        }
        return overlap;
}

void changeCoords(Sphere *sphere, double x, double y, double z)
{
        /*
         * Function:    changeCoords
         * -------------------------
         * Modifies the position data for a given particle
         * 
         * *disk:       pointer to the data relative to the particle of
         *              interest
         * x:           new x position of the particle
         * y:           new y position of the particle
         * z:           new z position of the particle
         *
         */
	sphere->x = x;
	sphere->y = y;
        sphere->z = z;
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
         * **spheres:     pointer to the pointer to the table of particles data
         */
        FILE *outfile = NULL;
        int i = 0;
        outfile = fopen(filename, "a");
        if (outfile != NULL)
        {
                fprintf(outfile, "&%d\n", MAX_SIZE);
                fprintf(outfile, "%d %d %d\n", (int)X_LENGTH, (int)Y_LENGTH, (int)Z_LENGTH);
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
         * **sphere:    pointer to a pointer to the data relevant for a given
         *              particle
         */ 
        double delta[3] = {genrand() - 0.5, genrand() - 0.5, genrand() - 0.5};
        int n = (int) (genrand() * MAX_SIZE);
        Sphere  bufferSphere = {        fmod((spheres[n]->x + delta[0]) + 2 * X_LENGTH, X_LENGTH),
                                        fmod((spheres[n]->y + delta[1]) + 2 * Y_LENGTH, Y_LENGTH),
                                        fmod((spheres[n]->z) + 2 * Z_LENGTH, Z_LENGTH),
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



int main(int argc, char *argv[])
{
        init_genrand(2467109);
        int i = 0, n = 0;
        Sphere *spheres = malloc(MAX_SIZE * sizeof(Sphere));
        Sphere **pSpheres = malloc(MAX_SIZE * sizeof(Sphere));
        for (int i = 0; i < MAX_SIZE; ++i)
                pSpheres[i] = spheres+i;
       
        /*
        for (i = 0; i < 100; i++)
        {
                n = abs((int) (genrand() * MAX_SIZE));
                if (overlapCheckGlobal(pSpheres, pSpheres[n], n))
                        printf("There is no overlap ! (as should be) (also n = %d)\n", n);
                else
                        printf("There is overlap, and it shouldn't ... (also n = %d)\n", n);

        }
        */

        // Initialization of the spheres coordinates one by one
        /*
        spheres[0] = (Sphere) {0.0f, 0.0f, 1.0f};
        spheres[1] = (Sphere) {0.0f, 1.9f, 1.0f};
        for (int i = 0; i < MAX_SIZE; i++)
               spheres[i] = (Sphere) {0.0f, 0.0f, 1.0f};
        */

        // Checks for moveAttempt()
        readInit("init.sph", pSpheres);
        writeCoords("coords.sph", pSpheres);
        for (i = 0; i < 100; i++)
        {
                moveAttempt(pSpheres);        
                writeCoords("coords.sph", pSpheres);
        }

        free(spheres);
        free(pSpheres);
        return 0;
}

