#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include "ecc.h"
#include "pollard_rho_worker.h"


void iter(ECC_ctx ctx, Pet p) {
    /*
        Q = xG
        U = aG + bQ
                    | G + U,    if H(U) = 0
        iter(U) =>  | 2U,       if H(U) = 1
                    | Q + U,    if H(U) = 2
    
        H(U) => U.x % 3
    */
    unsigned char H = p->U->x->_mp_d[0] % 3; // _mp_d[0] is the least significant limb

    if (H == 0) {
        // U = U + G
        pointAdd(ctx, p->U, ctx->G);
        // this means we have incremented a
        mpz_add_ui(p->a, p->a, 1);
        // we work modulo the order of G but it should never be necesseray to reduce
        // mpz_mod(p->a, p->a, ctx->q);
    }
    else if (H == 1) {
        // U = 2U
        pointDouble(ctx, p->U);
        // this means we have doubled a and b
        mpz_mul_ui(p->a, p->a, 2);
        mpz_mul_ui(p->b, p->b, 2);
        // we work modulo the order of G
        mpz_mod(p->a, p->a, ctx->q);
        mpz_mod(p->b, p->b, ctx->q);
    }
    else if (H == 2) {
        // U = U + Q
        pointAdd(ctx, p->U, ctx->Q);
        // this means we have incremented b
        mpz_add_ui(p->b, p->b, 1);
        // we work modulo the order of G but it should never be necesseray to reduce
        // mpz_mod(p->b, p->b, ctx->q);
    }
}


int pollard_rho_worker(mpz_t seed) {
    /*
        Parallelized Pollard's Rho algorithm worker process.
        Walk a pseudo-random sequence until a distinguished point is found.
    */
    ECC_ctx_t ctx;
    Point_t G, Q;
    Pet_t tortoise;
    gmp_randstate_t state;

    gmp_randinit_default(state); // mersenne twister by default, which is fast
	gmp_randseed(state, seed);

    // context initialization
    init_ctx(ctx, P, A, ORDER);
    init_point_coord(G, GX, GY);
    init_point_coord(Q, QX, QY);
    ctx->G = G;
    ctx->Q = Q;

    init_randomPet(ctx, state, tortoise);

    do {        
        iter(ctx, tortoise);
    // distinguished point has 24 zero LSB (this is arbitrary)
    } while (tortoise->U->x->_mp_d[0] % 0x1000000);

    print_Pet(tortoise);
    
    clear_point(G);
    clear_point(Q);
    clear_ctx(ctx);
    clear_Pet(tortoise);

    return 0;
}


int main(int argc, char *argv[]) {
    mpz_t seed;
    if (argc < 2) {
        // in case no seed is given, just use the current time
        mpz_init_set_ui(seed, time(NULL));
    }
    else {
        mpz_init_set_str(seed, argv[1], 10);
    }
    int ret = pollard_rho_worker(seed);
    mpz_clear(seed);
    return ret;
}