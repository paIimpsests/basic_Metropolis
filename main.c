/* First attempt at implemeting a basic Metropolis Scheme for Monte-Carlo
 * simulation of a hard sphere system
 * Antoine Castagnede | October 2022
 *
 * WHAT DO WE NEED ?
 * -----------------
 * (o)  Globally define a 2D box of given size --- Generalized to 3D
 *  -   Generate a suitable starting configuration --- ask Frank about that
 * (o)  Read starting configuration from a file
 * (o)  Function to calculate "energy" of a given particle --- atm only check
 *      for overlap
 * (o)  Function to correctly write positions of particles to a file in the
 *      fashion of Frank's visualization code
 *  -   Function to attempt to move a particle --- must call writeCoords(),
 *      changeCoords(), overlapCheck()
 * (o)  Random number generation --- meh... ask Frank about that
 *      it properly
 *  -   Function to randomly select a particle
 *  -   What kind of small displacement to use ?
 *  -   How many cycles of moving attemps before writing the new system
 *      configuration ? --- chose arbitrarily, depends on number of particles
 *      and displacement size
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
#include "splitmix64.h"
#include "xoshiro256plusplus.h"

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
	if (sqrt(pow(disk1->x - disk2->x, 2) + pow(disk1->y - disk2->y, 2) + pow(disk1->z - disk2->z, 2)) < disk1->r + disk2->r)
		overlap = 0;
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

int generateRandom()
{
        /*
         * Function: generateRandom
         * ------------------------
         * Genates a random integer between 1 and 100 using pseudo-random
         * generation based on xoshiro256plusplus and splitmix64 generators
         *
         * return:      integer between 1 and 100
         */
        srand(time(NULL));
        uint64_t s[4] = {(uint64_t) rand()}, n = 0;
        s[0] = next_splitmix(s[0]);
        s[1] = next_splitmix(s[0]);
        s[2] = next_splitmix(s[1]);
        s[3] = next_splitmix(s[2]);
        return (int) (next_xoshiro(s) % 100) + 1;
}


int main(int argc, char *argv[])
{
        int i = 0;
        Disk *disks = malloc(MAX_SIZE * sizeof(Disk));
        Disk **pdisks = malloc(MAX_SIZE * sizeof(Disk));
        for (int i = 0; i < MAX_SIZE; ++i)
                pdisks[i] = disks+i;
        
        // Initialisation of the disks coordinates one by one
        /*
        disks[0] = (Disk) {0.0f, 0.0f, 1.0f};
        disks[1] = (Disk) {0.0f, 1.9f, 1.0f};
        for (int i = 0; i < MAX_SIZE; i++)
               disks[i] = (Disk) {0.0f, 0.0f, 1.0f};
        */

        //readInit("init.txt", pdisks);
        //writeCoords("coords.txt", pdisks);
        //printf("%lf\n", pdisks[0]->r);
        printf("Pseudo-randomly generated number between 1 and 100 : %d\n", generateRandom());

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

