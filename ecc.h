#ifndef ECCH
#define ECCH
#include <gmp.h>


// EC point 
typedef struct {
	mpz_t x;
	mpz_t y;
	mpz_t z;
} Point_t[1], *Point;

// ECC context
typedef struct {
    // curve modulus
	mpz_t p;
    // curve weierstrass a
	mpz_t a;
    // Generator order
	mpz_t q;
    // Generator
    Point G;
    // Public point
    Point Q;
} ECC_ctx_t[1], *ECC_ctx;

void init_point_coord(Point P, const char *x, const char *y);
void copy_point(Point P, Point Q);
void clear_point(Point P);

void init_ctx(ECC_ctx ctx, const char *p, const char *a, const char *q);
void clear_ctx(ECC_ctx ctx);

int pointEqual(Point A, Point B);
int isInfinity(Point A);

void pointDouble(ECC_ctx ctx, Point P, mpz_t r);
void pointDouble_slow(ECC_ctx ctx, Point P);
void pointAdd(ECC_ctx ctx, Point P, Point Q, mpz_t r);
void pointAdd_slow(ECC_ctx ctx, Point P, Point Q);
void pointMul(ECC_ctx ctx, Point P, mpz_t k);

#endif