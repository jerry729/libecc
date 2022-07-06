/*
定义ECC域需要的参数
1. 有限域大小p(素数)
2. 椭圆曲线系数 a, b
3. 子群基准点G
4. 子群的阶n
5. 子群的协因子h

a, b, p
有限域上椭圆曲线的阶N schoof's algorithm
h = N / n 算术除法
G = hP 椭圆曲线标量积


私钥[1, n-1]随机选d
公钥H = dG
*/

#ifndef ECC_H
#define ECC_H

#include "libtommath/tommath.h"
#include <stdio.h>
#include <string.h>

#define PRIME_FIELD_LEN 200 //bit
#define PRIVATE_KEY_LEN 128 //bit
#define BYTE_ARRAY_MAX_LEN 800

typedef struct
{
    mp_int _x;
    mp_int _y;
}Point;

/*
y^2 ≡ x^3 + ax + b (mod p)
*/
typedef struct
{
    mp_int _p;
    mp_int _a;
    mp_int _b;
}EllipticCurve;


/*
    n 为素数
    p 为素数 越大越安全，但平衡计算速度，选择200bit左右
    p != h * n 已知可以用smart方法攻击
    4a^3 + 27b^2 != 0 (mod p)
    h <= 4 使得n足够大 (h=N/n)
    pt != 1 (mod n), 1<=t<20
*/


int ecc(EllipticCurve *ec, Point *G, mp_int *n, mp_int *h);

unsigned long random_gen(unsigned long m,unsigned long n);

/*
get finite field prime
*/
int get_prime(mp_int *p, int len);

/*
get the order N of elliptic curve
*/
int schoofs(mp_int *N, EllipticCurve *ec);

/**/
mp_int init_all();

/**/
int all_primes(mp_int l[], mp_int* q);

int is_prime(const mp_int *x);

int MR_is_prime(const mp_int *x);

/*
param: A and p
ret: B and G(x, y)
*/
int get_B_n_G(mp_int *A, mp_int *p, mp_int *B, Point *G);

int scalar_mul(Point *pub_K, Point *G, const EllipticCurve *ec, mp_int *d);

/*R <- P + Q
if zero is true, always make q as the zero elment
*/
int points_add(Point *p, Point *q, Point *r, const EllipticCurve *ec, bool *zero);

/**
 * @brief deep copy a -> b
 * 
 * @param a 
 * @param b 
 * @return int 
 */
int point_copy(Point *a, Point *b);

/**
 * @brief initialize a(x, y) by mp_init
 * 
 * @param a 
 * @return int 
 */
int point_init(Point *a);

/**
 * @brief init a(x,y) by mp_init_set, 
 * x,y is digit
 * 
 * @param a 
 * @param x 
 * @param y 
 * @return int 
 */
int point_init_set(Point *a, mp_digit x, mp_digit y);

/**
 * @brief init a and copy a <- b
 * 
 * @param a 
 * @param b 
 * @return int 
 */
int point_init_copy(Point *a, Point *b);

#endif //ECC_H