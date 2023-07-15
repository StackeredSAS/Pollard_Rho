#ifndef PRH
#define PRH
#include <gmp.h>
#include "ecc.h"

/* Parameters must be supplied as HEX strings */
// curve modulus
#define P "92ddcf14cb9e71f4489a2e9ba350ae29454d98cb93bdbcc07d62b502ea12238ee904a8b20d017197aae0c103b32713a9"
// curve weierstrass A
#define A "01"
// curve weierstrass B is not needed
// Generator coordinates
#define GX "46E3775ECE21B0898D39BEA57050D422A0AF989E497962BAEE2CB17E0A28D5360D5476B8DC966443E37A14F1AEF37742"
#define GY "7C8E741D2C34F4478E325469CD491603D807222C9C4AC09DDB2B31B3CE3F7CC191B3580079932BC6BEF70BE27604F65E"
// Generator order
#define ORDER "DB6B4C58EFBAFD"
// Public point coordinates
#define QX "5D8DBE75198015EC41C45AAB6143542EB098F6A5CC9CE4178A1B8A1E7ABBB5BC64DF64FAF6177DC1B0988AB00BA94BF8"
#define QY "23A2909A0B4803C89F910C7191758B48746CEA4D5FF07667444ACDB9512080DBCA55E6EBF30433672B894F44ACE92BFA"

#define NPETS 100
// distinguished point has 24 zero LSB (this is arbitrary)
#define DIST (1u << 24)

mpz_t bs[NPETS];

// Pet (turtle, snail, hare)
typedef struct {
    // U = aG + bQ
    Point_t U;
	mpz_t a;
	mpz_t b;
    mpz_t r;
} Pet_t[1], *Pet;

void init_randomPet(ECC_ctx ctx, gmp_randstate_t state, Pet p) {
    /*
    a = random()
    b = random()
    U = aG + bQ
    */
    Point_t tmp;
    mpz_init(p->a);
    mpz_init(p->b);
    mpz_init(p->r);
    mpz_urandomm(p->a, state, ctx->q);
    mpz_urandomm(p->b, state, ctx->q);
    copy_point(p->U, ctx->G); // U = G
    copy_point(tmp, ctx->Q); // tmp = Q
    pointMul(ctx, p->U, p->a); // U = aG
    pointMul(ctx, tmp, p->b); // tmp = bQ
    pointAdd_slow(ctx, p->U, tmp); // U = aG + bQ
    clear_point(tmp);
}

void clear_Pet(Pet p) {
    clear_point(p->U);
    mpz_clear(p->a);
    mpz_clear(p->b);
    mpz_clear(p->r);
}

void print_Pet(Pet p) {
    // x,a,b
    mpz_out_str(stdout, 16, p->U->x);
    printf(",");
    mpz_out_str(stdout, 16, p->a);
    printf(",");
    mpz_out_str(stdout, 16, p->b);
    printf("\n");
}

#endif