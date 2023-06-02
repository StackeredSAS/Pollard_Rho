#include <stdio.h>
#include "ecc.h"

void init_point_coord(Point P, const char *x, const char *y) {
    // Init the point with the given coordinates in HEX
    mpz_init_set_ui(P->z, 1); // always 1, when not the point at infinity
    mpz_init_set_str(P->x, x, 16);
    mpz_init_set_str(P->y, y, 16);
}

void copy_point(Point P, Point Q) {
    // Init the point with the values of Q
    mpz_init_set(P->z, Q->z);
    mpz_init_set(P->x, Q->x);
    mpz_init_set(P->y, Q->y);
}

void clear_point(Point P) {
    // Clear the point
    mpz_clear(P->z);
    mpz_clear(P->x);
    mpz_clear(P->y);
}

void init_ctx(ECC_ctx ctx, const char *p, const char *a, const char *q) {
    // Init the ECC context with the given parameteres in HEX
    mpz_init_set_str(ctx->p, p, 16);
    mpz_init_set_str(ctx->a, a, 16);
    mpz_init_set_str(ctx->q, q, 16);
}

// Clears an ECC context
void clear_ctx(ECC_ctx ctx) {
    // Clear the context
    mpz_clear(ctx->p);
    mpz_clear(ctx->a);
    mpz_clear(ctx->q);
}

int pointEqual(Point A, Point B) {
    // only car about X and Y for this use case
    return mpz_cmp(A->x, B->x) == 0 && mpz_cmp(A->y, B->y) == 0;
}

int isInfinity(Point A) {
    return mpz_cmp_ui(A->z, 0) == 0 && mpz_cmp_ui(A->y, 1) == 0;
}

void pointDouble(ECC_ctx ctx, Point P) {
    // P = P*2
    mpz_t lambda, tmp;

    if (isInfinity(P)) {
        // do nothing
        return;
    }

    // compute lambda = 3x^2 + a / 2y mod p
    mpz_init_set(lambda, ctx->a); // lambda = a
    mpz_init_set(tmp, P->x); // tmp = x

    mpz_powm_ui(tmp, tmp, 2, ctx->p); // tmp = x^2 mod p
    mpz_addmul_ui(lambda, tmp, 3); // lambda = 3x^2 + a
    mpz_mod(lambda, lambda, ctx->p); // lambda = 3x^2 + a mod p

    mpz_set(tmp, P->y); // tmp = y
    mpz_mul_ui(tmp, tmp, 2); // tmp = 2y
    mpz_invert(tmp, tmp, ctx->p); // tmp = 2y^-1 mod p

    mpz_mul(lambda, lambda, tmp); // lamda = 3x^2 + a / 2y
    mpz_mod(lambda, lambda, ctx->p); // lambda = 3x^2 + a / 2y mod p

    // New X
    mpz_powm_ui(tmp, lambda, 2, ctx->p); // tmp = l^2 mod p
    mpz_sub(tmp, tmp, P->x); // tmp = l^2 - x
    mpz_sub(tmp, tmp, P->x); // tmp = l^2 - 2x
    mpz_mod(tmp, tmp, ctx->p); // tmp = l^2 - 2x mod p

    // New Y
    mpz_sub(P->x, P->x, tmp); // x = x - tmp (we do this to avoid coping data to a new variable)
    mpz_mul(P->x, lambda, P->x); // x = l(x - tmp)
    mpz_sub(P->y, P->x, P->y); // y = l(x - tmp) - y
    mpz_mod(P->y, P->y, ctx->p); // y = l(x - tmp) - y mod p

    // store new X
    mpz_mod(P->x, tmp, ctx->p);

    mpz_clear(lambda);
    mpz_clear(tmp);
}

void pointAdd(ECC_ctx ctx, Point P, Point Q) {
    // P = P + Q
    mpz_t lambda, tmp;

    if (pointEqual(P, Q)) {
        // Double and return
        pointDouble(ctx, P);
        return;
    }

    if (isInfinity(Q)) {
        // do nothing
        return;
    }

    if (isInfinity(P)) {
        // return Q
        mpz_set(P->x, Q->x);
        mpz_set(P->y, Q->y);
        mpz_set(P->z, Q->z);
        return;
    }

    // compute lambda = (Yq - Yp) / (Xq - Xp)mod p
    mpz_init_set(lambda, Q->x); // lambda = Xq
    mpz_init_set(tmp, Q->y); // tmp = Yq

    mpz_sub(lambda, lambda, P->x); // lambda = Xq - Xp
    mpz_invert(lambda, lambda, ctx->p); // lambda = (Xq - Xp)^-1 mod p
    mpz_sub(tmp, tmp, P->y); // tmp = Yq - Yp
    
    mpz_mul(lambda, lambda, tmp); // lamda = (Xq - Xp) / (Yq - Yp)
    mpz_mod(lambda, lambda, ctx->p); // lambda =(Xq - Xp) / (Yq - Yp) mod p

    // New X
    mpz_powm_ui(tmp, lambda, 2, ctx->p); // tmp = l^2 mod p
    mpz_sub(tmp, tmp, P->x); // tmp = l^2 - Xp
    mpz_sub(tmp, tmp, Q->x); // tmp = l^2 - Xp - Xq
    mpz_mod(tmp, tmp, ctx->p); // tmp = l^2 - Xp - Xq mod p

    // New Y
    mpz_sub(P->x, P->x, tmp); // x = x - tmp (we do this to avoid coping data to a new variable)
    mpz_mul(P->x, lambda, P->x); // x = l(x - tmp)
    mpz_sub(P->y, P->x, P->y); // y = l(x - tmp) - y
    mpz_mod(P->y, P->y, ctx->p); // y = l(x - tmp) - y mod p

    // store new X
    mpz_mod(P->x, tmp, ctx->p);

    mpz_clear(lambda);
    mpz_clear(tmp);
}

void pointMul(ECC_ctx ctx, Point P, mpz_t k) {
    // P = P*k
    mp_bitcnt_t index;
	mp_size_t size, i;
    Point_t tmp;
    size = mpz_sizeinbase(k, 2);
    mpz_init(tmp->x);
    mpz_init(tmp->z);
    mpz_init_set_ui(tmp->y, 1);
   
	for (i = 0; i<size; i++) {
		index = mpz_scan1(k, i);
		if (index == i) {
			pointAdd(ctx, tmp, P);
		}
        pointDouble(ctx, P);
	}

    mpz_set(P->x, tmp->x);
    mpz_set(P->y, tmp->y);
    clear_point(tmp);
}