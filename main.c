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
 *  -   Function to attempt to move a particle --- WIP
 * (o)  Random number generation --- implemented using mt19937ar.c
 * (o)  Function to check for overlap with the whole system
 * (o)  Function to randomly select a particle --- already implemented in
 *      moveAttempt()
 *  -   What kind of small displacement to use ? --- depends on box size, radii,
 *      see Frenkel & Smit 3.2.2. Boundary Conditions
 *  -   How many cycles of moving attemps before writing the new system
 *      configuration ? --- chose arbitrarily, depends on number of particles
 *      and displacement size
 *  -   Implement the same using a cell list (makes it faster)
 *  -   How to deal with box borders ? 
 *  -   Follow up : implement periodic boundary conditions
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
#include <inttypes.h>
#include <time.h>
#include <math.h>
#include "mt19937ar.c"

#define MAX_SIZE 100
#define X_LENGTH 10.034544
#define Y_LENGTH 10.034544
#define Z_LENGTH 1

typedef struct Disk Disk;
struct Disk
{
	double x;	// x position of the center of the disk	
	double y;	// y position of the center of the disk
        double z;       // z position of the center of the disk
	double r;	// radius of the disk
        char type;      // particle type --- for visualization purposes only

};

int overlapCheck(Disk *disk1, Disk *disk2)
{
        /*
         * Function:    overlapCheck
         * -------------------------
         * Checks for overlap between two given particules by computing the
         * distance between their centers and comparing it to the sum of their
         * radii
         *
         * *disk1:      pointer to the data relative to the first particle
         * *disk2:      pointer to the data relative to the second particle
         *
         * return:      1 is there is overlap, 0 if not
         */
	int overlap = 1;
	if (((disk1->x - disk2->x) * (disk1->x - disk2->x) + (disk1->y - disk2->y) * (disk1->y - disk2->y) + (disk1->z - disk2->z) * (disk1->z - disk2->z)) < ((disk1->r + disk2->r) * (disk1->r + disk2->r)))
		overlap = 0;
	return overlap;
}

int overlapCheckGlobal(Disk **disks, Disk *sample, int n)
{
        /* Function:    overlapCheckGlobal
         * -------------------------------
         *  Checks for overlap between one given particle and the rest of the
         *  system
         *
         *  **disks:    pointer to a pointer to the data of the system
         *  n:          location of the particle that we want to check overlap
         *              for, 0 <= n <= 99
         *
         *  return:     1 if there is overlap, 0 if not
         */
        int overlap = 0, i = 0;
        while ((overlap == 0) && (i < MAX_SIZE))
        {
                overlap = !(i == n) * overlapCheck(sample, disks[i]);
                i++;
        }
        return overlap;
}

void changeCoords(Disk *disk, double x, double y, double z)
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
	disk->x = x;
	disk->y = y;
        disk->z = z;
}

void readInit(char *filename, Disk **disks)
{
        /*
         * Function:    readInit
         * ---------------------
         * Initializes the table of data of all MAX_SIZE particles by reading
         * it from a user supplied .txt file
         * NB1: Some parameters (total number of particles and box size) from
         * the initialization .txt file must match those defined here
         * NB2: Reading is done only for files with typesetting indicated above
         *
         * *filename:   pointer to the name of the initialization .txt file
         * **disks:     pointer to the pointer to the table of particles data
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
                             fscanf(initfile, "%c %lf %lf %lf %lf%*c", &disks[i]->type, &disks[i]->x, &disks[i]->y, &disks[i]->z, &disks[i]->r);
                        }

                }
                else
                        printf("Initialization did not happen, check the initialization source file.\n"); 
                fclose(initfile);
        }
}

void writeCoords(char *filename, Disk **disks)
{
        /*
         * Function:    writeCoords
         * -------------------------
         * Writes the position, radius, and type data for all MAX_SIZE
         * particles in a .txt file
         * NB: Writing is done using the typesetting indicated above
         *
         * *filename:   pointer to the name of the initialization .txt file
         * **disks:     pointer to the pointer to the table of particles data
         */
        FILE *outfile = NULL;
        int i = 0;
        outfile = fopen(filename, "a");
        if (outfile != NULL)
        {
                fprintf(outfile, "&%d\n", MAX_SIZE);
                fprintf(outfile, "%d %d %d\n", (int)X_LENGTH, (int)Y_LENGTH, (int)Z_LENGTH);
                for (i = 0; i < MAX_SIZE; i++)
                        fprintf(outfile, "%c %lf %lf %lf %lf\n", disks[i]->type, disks[i]->x, disks[i]->y, disks[i]->z, disks[i]->r);

        }
}

double gendouble(int *pseed)
{
        /*
         * Function:    gendouble
         * ----------------------
         * Random generator wrapper for a double in [0,1[, it also resets the
         * speed from the randomly generated number
         * It employs the mt19937ar.c generator based on the Mersenne Twister
         * algorithm
         *
         * *pseed:      pointer to the seed for the generator
         *
         * return: randomly generated double in [0,1[
         */
        double m = genrand(), n = 0.0f;       
        *pseed = (int) m * 1000000;
        n = genrand();
        *pseed = (int) n * 1000000;
        return n > 0.5 ? m : -m;
}

void moveAttempt(Disk **disks, int *pseed)
{
        /*
         * Function: moveAttempt                
         * |
         * | WARNING: implementation not finished - unsure if working as intended atm
         * |
         * ---------------------
         * Tries to move a randomly selected particle in space by a small amount
         * and checks for overlap
         * If there is no overlap, overwrites the previous particle position
         * with the new one
         *
         * *disk:       pointer to the data relevant for a given particle
         * *pseed:      pointer to the seed for the random number generator
         */ 
        double delta[3] = {gendouble(pseed) / 100, gendouble(pseed) / 100, gendouble(pseed) / 100};
        int n = abs((int) (gendouble(pseed) * 1000000) % 100);
        Disk bufferDisk = {disks[n]->x + delta[0], disks[n]->y + delta[1], disks[n]->z + delta[2], disks[n]->r, disks[n]->type}, *pBufferDisk = NULL;
        pBufferDisk = &bufferDisk;
        printf("%lf\n", delta[0]);
        printf("%lf\n", delta[1]);
        printf("%lf\n", delta[2]);
        printf("%lf\n", bufferDisk.x);
        printf("%lf\n", bufferDisk.y);
        printf("%lf\n", bufferDisk.z);
        printf("%lf\n", bufferDisk.r);
        printf("%c\n", bufferDisk.type);
        printf("n = %d\n", n);

        if (!overlapCheckGlobal(disks, pBufferDisk, n))
        {
                disks[n]->x = bufferDisk.x;
                disks[n]->y = bufferDisk.y;
                disks[n]->z = bufferDisk.z;
                printf("It moves !\n");
        }
        printf("%lf\n", disks[n]->x);
        printf("%lf\n", disks[n]->y);
        printf("%lf\n", disks[n]->z);
        printf("%lf\n", disks[n]->r);
        printf("%c\n", disks[n]->type);
}



int main(int argc, char *argv[])
{
        srand(time(NULL));
        int i = 0, n = 0, p = 0, seed = 0, *pseed = &seed;
        double m = 0.0f;
        Disk *disks = malloc(MAX_SIZE * sizeof(Disk));
        Disk **pdisks = malloc(MAX_SIZE * sizeof(Disk));
        for (int i = 0; i < MAX_SIZE; ++i)
                pdisks[i] = disks+i;
       
        // Initialization of the random number generator
        *pseed = rand();
        init_genrand(*pseed);

        /*
        for (i = 0; i < 100; i++)
        {
                n = abs((int) (gendouble(pseed) * 1000000) % 100);
                if (overlapCheckGlobal(pdisks, n))
                        printf("There is no overlap ! (as should be) (also n = %d)\n", n);
                else
                        printf("There is overlap, and it shouldn't ... (also n = %d)\n", n);

        } 
        */

        // Initialization of the disks coordinates one by one
        /*
        disks[0] = (Disk) {0.0f, 0.0f, 1.0f};
        disks[1] = (Disk) {0.0f, 1.9f, 1.0f};
        for (int i = 0; i < MAX_SIZE; i++)
               disks[i] = (Disk) {0.0f, 0.0f, 1.0f};
        */

        readInit("init.sph", pdisks);
        //writeCoords("coords.sph", pdisks);
        

        moveAttempt(pdisks, pseed);

        // Simple tests for overlap check
        /*
	if (overlapCheck(pdisks[0], pdisks[1]))
		printf("Il n'y a pas d'overlap !\n");
	else
		printf("Il y a overlap !\n");

	changeCoords(pdisks[1], 2.0f, 2.0f, 2.0f);
	if (overlapCheck(pdisks[0], pdisks[1]))
		printf("Il n'y a pas d'overlap !\n");
	else
		printf("Il y a overlap !\n");
        */ 

        free(disks);
        free(pdisks);
        return 0;
}

