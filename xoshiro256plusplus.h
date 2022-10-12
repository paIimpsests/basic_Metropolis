// Header for xoshiro256plusplus.c

#ifndef xoshiro256plusplus_c
#define xoshiro256plusplus_c

#include <stdint.h>

static inline uint64_t rotl(const uint64_t x, int k);
uint64_t next_xoshiro(uint64_t s[]);
void jump(void);
void long_jump(void);

#endif
