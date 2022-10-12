all: splitmix64.c xoshiro256plusplus.c main.c
	gcc -o output splitmix64.c xoshiro256plusplus.c main.c -lm

