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

/*
// Toy values for Q = 133337*G
// Solving for the above values takes approximately 1h30 on average (can be less if you're lucky)
#define GX "7f712cb991cec4aa1b8314c8c95e2443eb77e2803c3db09ac759b67324e66cbe1ad5cff08a008658de90142cdd8a892d"
#define GY "3e962b59ef58dea1aff4b629caffdf671b64ea6d98ad55f6cc123a17e602f10ec510349dfb98a1cbb6de2432cc37ec25"
#define ORDER "8df01"
#define QX "4906867b4742830958d1608adab80d0745b71176ee11f6bb8ab11b238adc484fd82add2d35eed1f96cb65da2dfcb4359"
#define QY "197dede35023959ee69abdc256981dbaac6259f53ffe330e67d52a40b2b7b3be3c7eddbc7c6820b2186c96b5cd228c54"
*/

// Maximum number of pets to remember
// only powers of 2 are remembered, so 50 allows to go up to 2^50 iterations
#define PETS 50

// Pet (turtle, snail, hare)
typedef struct {
    // U = aG + bQ
    Point_t U;
	mpz_t a;
	mpz_t b;
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
    mpz_urandomm(p->a, state, ctx->q);
    mpz_urandomm(p->b, state, ctx->q);
    copy_point(p->U, ctx->G); // U = G
    copy_point(tmp, ctx->Q); // tmp = Q
    pointMul(ctx, p->U, p->a); // U = aG
    pointMul(ctx, tmp, p->b); // tmp = bQ
    pointAdd(ctx, p->U, tmp); // U = aG + bQ
    clear_point(tmp);
}

void copy_Pet(Pet p, Pet q) {
    copy_point(p->U, q->U);
    mpz_init_set(p->a, q->a);
    mpz_init_set(p->b, q->b);
}

void clear_Pet(Pet p) {
    clear_point(p->U);
    mpz_clear(p->a);
    mpz_clear(p->b);
}

#endif