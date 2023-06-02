#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include "ecc.h"
#include "pollard_rho.h"


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
        // we work modulo the order of G
        mpz_mod(p->a, p->a, ctx->q);
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
        // we work modulo the order of G
        mpz_mod(p->b, p->b, ctx->q);
    }
}


int found(ECC_ctx ctx, Pet a, Pet_t *table, int total) {
    /*
        Check if a cycle has been found and compute the private key.
    */
    mpz_t priv;
    // Check if the turtle has matched with any snail
    for (int i=0; i<total; i++) {
        if (pointEqual(a->U, table[i]->U)) {

                mpz_init(priv);
                mpz_sub(priv, a->a, table[i]->a); // priv = a - a_s
                mpz_sub(a->a, table[i]->b, a->b); // a = b_s - b
                mpz_invert(a->a, a->a, ctx->q); // a = (b_s - b)^-1 mod q
                mpz_mul(priv, a->a, priv); // priv = (a - a_s) / (b_s - b)
                mpz_mod(priv, priv, ctx->q); // priv = (a - a_s) / (b_s - b) mod q

                // print the key
                printf("Found private key !\nx = ");
                mpz_out_str(stdout, 10, priv);
                printf("\n");

                mpz_clear(priv);

            return 1;
        }
    }
    return 0;
}

int pollard_rho_Brent(unsigned long int seed) {
    /*
        Pollard's Rho algorithm using Brent's cycle-finding algorithm
    */
    ECC_ctx_t ctx;
    Point_t G, Q;
    Pet_t tortoise, snail;
    Pet_t table[PETS];
    // loop variables
    mpz_t nextPower, i, stop;
    int j = 0, return_code = 0;
    gmp_randstate_t state;

    gmp_randinit_default(state); // mersenne twister by default, which is fast
	gmp_randseed_ui(state, seed);

    // context initialization
    init_ctx(ctx, P, A, ORDER);
    init_point_coord(G, GX, GY);
    init_point_coord(Q, QX, QY);
    ctx->G = G;
    ctx->Q = Q;

    init_randomPet(ctx, state, tortoise);
    // snail = tortoise
    copy_Pet(snail, tortoise);
    
    // Main loop
    mpz_init_set_ui(nextPower, 1);
    mpz_init(i);
    mpz_init(stop);
    // if we have not found anything after 3*sqrt(ORDER), simply retry, it should be faster than continuing
    // this bound is purely arbitrary, but it seems to yield good results
    // The bound must be chosen so it offers a good probability of finding a cycle in less this amount of tries, 
    // while being the lowest possible.
    mpz_sqrt(stop, ctx->q);
    mpz_mul_ui(stop, stop, 3);
    
    do {

        // Fail condition
        if (mpz_cmp(i, stop) >= 0) {
            return_code = -1;
            goto end;
        }
        
        // We store the snail values of Brent's algorithm
        if (mpz_cmp(i, nextPower) == 0) {
            // power of two, store previous tortoise (snail)
            copy_Pet(table[j], tortoise);
            // compute next power
            mpz_mul_ui(nextPower, nextPower, 2);
            // keep track of the number of snails
            j++; // should never be >= PETS
        }
        
        // move tortoise one step
        iter(ctx, tortoise);
        mpz_add_ui(i, i, 1);

    } while (!found(ctx, tortoise, table, j));

    printf("Found log in ");
    mpz_out_str(stdout, 10, i);
    printf(" iterations.\n");
    
end:
    clear_point(G);
    clear_point(Q);
    clear_ctx(ctx);
    clear_Pet(tortoise);
    clear_Pet(snail);
    mpz_clear(i);
    mpz_clear(stop);
    mpz_clear(nextPower);

    return return_code;
}


int main(int argc, char *argv[]) {
    unsigned long int seed;
    if (argc < 2) {
        // in case no seed is given, just use the current time
        seed = time(NULL);
    }
    else {
        seed = atoll(argv[1]);
    }
    return pollard_rho_Brent(seed);
}