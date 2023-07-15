#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include "ecc.h"
#include "pollard_rho_worker.h"


void batch_invert(Pet_t p[], mpz_t prime, int size) {
    for (size_t i = 0; i < size; i++)
    {
        mpz_invert(p[i]->r, p[i]->r, prime);
    }
}

void iter_pre(ECC_ctx ctx, Pet p) {
    /*
        Just compute the value we will need to invert during the batch inversion process.
    */
    unsigned char H = p->U->x->_mp_d[0] % 3; // _mp_d[0] is the least significant limb

    if (H == 0) {
        // U = U + G
        // (Xq - Xp)^-1 mod p
        mpz_sub(p->r, ctx->G->x, p->U->x);
    }
    else if (H == 1) {
        // U = 2U
        // 2y^-1 mod p
        mpz_mul_ui(p->r, p->U->y, 2);
    }
    else if (H == 2) {
        // U = U + Q
        // (Xq - Xp)^-1 mod p
        mpz_sub(p->r, ctx->Q->x, p->U->x);
    }
}

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
        pointAdd(ctx, p->U, ctx->G, p->r);
        // this means we have incremented a
        mpz_add_ui(p->a, p->a, 1);
        // we work modulo the order of G but it should never be necesseray to reduce
        // mpz_mod(p->a, p->a, ctx->q);
    }
    else if (H == 1) {
        // U = 2U
        pointDouble(ctx, p->U, p->r);
        // this means we have doubled a and b
        mpz_mul_ui(p->a, p->a, 2);
        mpz_mul_ui(p->b, p->b, 2);
        // we work modulo the order of G
        mpz_mod(p->a, p->a, ctx->q);
        mpz_mod(p->b, p->b, ctx->q);
    }
    else if (H == 2) {
        // U = U + Q
        pointAdd(ctx, p->U, ctx->Q, p->r);
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
    int NPETS = 100; // 1: 5_720_567  2: 11_441_133   100: 13_603_690
    ECC_ctx_t ctx;
    Point_t G, Q;
    Pet_t pets[NPETS];
    gmp_randstate_t state;
    unsigned int count = 0;
    size_t i;
    int found;

    gmp_randinit_default(state); // mersenne twister by default, which is fast
	gmp_randseed(state, seed);

    // context initialization
    init_ctx(ctx, P, A, ORDER);
    init_point_coord(G, GX, GY);
    init_point_coord(Q, QX, QY);
    ctx->G = G;
    ctx->Q = Q;

    // initialize all pets
    for (i = 0; i < NPETS; i++) {
        init_randomPet(ctx, state, pets[i]);
    }
    
    do {
        // no dinstinguished point found in this batch
        found = 0;

        // compute the value to invert for all pets
        for (i = 0; i < NPETS && !found; i++) {
            iter_pre(ctx, pets[i]);
        }

        // compute the inverses in batch to speed things up
        batch_invert(pets, ctx->p, NPETS);

        // iter all pets and check for distinguised point
        for (i = 0; i < NPETS && !found; i++) {
            iter(ctx, pets[i]);
            // distinguished point has 24 zero LSB (this is arbitrary)
            found |= (pets[i]->U->x->_mp_d[0] % 0x1000000 == 0);
            count++;
        }

    // next batch
    } while (!found);

    print_Pet(pets[i-1]); // i-1 points to the distinguished point Pet
    printf("nb iter: %u\n", count);
    
    clear_point(G);
    clear_point(Q);
    clear_ctx(ctx);
    // clear all pets
    for (i = 0; i < NPETS && !found; i++) {
        clear_Pet(pets[i]);
    }

    return 0;
}


int main(int argc, char *argv[]) {
    mpz_t seed;
    if (argc < 2) {
        // in case no seed is given, just use the current time
        // mpz_init_set_ui(seed, time(NULL));
        mpz_init_set_ui(seed, 11);
    }
    else {
        mpz_init_set_str(seed, argv[1], 10);
    }
    int ret = pollard_rho_worker(seed);
    mpz_clear(seed);
    return ret;
}