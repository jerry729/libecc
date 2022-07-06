#include "ecc.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>


#define IS_PRIME_TRUE 2
#define IS_PRIME_FALSE 3
#define NO_PRIMES 2

unsigned long random_gen(unsigned long m,unsigned long n)  
{  
    srand((unsigned long)time(NULL));  
    return(unsigned long)(m+rand()%n);  
}  


int get_prime(mp_int *p, int len){

    mp_err result;

    if((result = mp_prime_rand(p, 10, len, (rand()&1) ? : MP_PRIME_2MSB_ON)) != MP_OKAY){
        printf("Error getting radom prime. %s", mp_error_to_string(result));\
        return EXIT_FAILURE;
    }

    return MP_OKAY;
}

int schoofs(mp_int *N, EllipticCurve *ec){
    mp_int M, t;
    mp_err result;

    if((result = mp_init_multi(&M, &t, NULL)) != MP_OKAY){
        printf("Error in initializing the number. %s", mp_error_to_string(result));
        return EXIT_FAILURE;
    }

    mp_set(&M, 1);
    mp_set(&t, 0);
    mp_int sqrtp;
    if((result = mp_init(&sqrtp)) != MP_OKAY){
        printf("Error in initializing the number. %s", mp_error_to_string(result));
        return EXIT_FAILURE;
    }

    if((result = mp_root_n(&(ec->_p), 2, &sqrtp)) != MP_OKAY){
        printf("Error in getting sqrt of p. %s", mp_error_to_string(result));
        return EXIT_FAILURE;
    } //p^(1/2)

    if((result = mp_mul_2d(&sqrtp, 2, &sqrtp)) != MP_OKAY){
        printf("Error in getting 4 times sqrtp. %s", mp_error_to_string(result));
        return EXIT_FAILURE;
    } // 4 * p^(1/2)

    while(mp_cmp(&M, &sqrtp)){
        //TODO: calculate trace of frobenius
        //      then using Chinese remainder theorem
    }

    return MP_OKAY;
}


int all_primes(mp_int l[], mp_int* q){
    if(mp_cmp_d(q, 2) == MP_LT){
        return NO_PRIMES;
    }

    mp_int i, cnt;
    mp_err result;

    if((result = mp_init_set(&i, 2)) != MP_OKAY){
        printf("Error in initializing the iterator. %s", mp_error_to_string(result));
        return EXIT_FAILURE;
    }

    if((result = mp_init_set(&cnt, 0)) != MP_OKAY){
        printf("Error in initializing the counter. %s", mp_error_to_string(result));
        return EXIT_FAILURE;
    }

    while(mp_cmp(&i, q) == MP_LT || mp_cmp(&i, q) == MP_EQ){
        if(MR_is_prime(&i) == IS_PRIME_TRUE){
            u_int64_t count = mp_get_mag_u64(&cnt);
            l[count] = i;
            if((result = mp_add_d(&cnt, 1, &cnt)) != MP_OKAY){
                printf("Error when adding one to the counter. %s", mp_error_to_string(result));
                return EXIT_FAILURE;
            }
        }

        if((result = mp_add_d(&i, 1, &i)) != MP_OKAY){
            printf("Error when adding 1 to the iterator. %s", mp_error_to_string(result));
            return EXIT_FAILURE;
        }
    }
    mp_clear_multi(&i, &cnt, NULL);
    return MP_OKAY;
}

int MR_is_prime(const mp_int *x){
    size_t size;
    size = mp_ubin_size(x) * 8;
    int times;
    bool res;
    mp_err err;
    times = mp_prime_rabin_miller_trials(size);
    unsigned long ul_x = mp_get_ul(x), base;
    base = random_gen(2, ul_x);

    mp_int b;
    if((err = mp_init_set(&b, base)) != MP_OKAY){
         printf("Error in initializing the base of Miller-Rabin. %s", mp_error_to_string(err));
        return EXIT_FAILURE;
    }

    for(int i = 0; i < times; i++){
        if((err = mp_prime_miller_rabin(x, &b, &res)) != MP_OKAY){
            printf("Error in a Miller-Rabin test. %s", mp_error_to_string(err));
        return EXIT_FAILURE;
        }
        if(res == false){
            return IS_PRIME_FALSE;
        }
    }

    mp_clear(&b);
    return IS_PRIME_TRUE;
}

int is_prime(const mp_int *x){
    mp_int i, i2, mod_red;
    mp_err result;
    if((result = mp_init_set(&i, 2)) != MP_OKAY){
        printf("Error in initializing the iterator. %s", mp_error_to_string(result));
        return EXIT_FAILURE;
    }

    if((result = mp_init_multi(&i2, mod_red, NULL)) != MP_OKAY){
        printf("Error in initializing the i^2. %s", mp_error_to_string(result));
        return EXIT_FAILURE;
    }

    if((result = mp_sqr(&i, &i2) != MP_OKAY)){
        printf("Error when calculating the i^2. %s", mp_error_to_string(result));
        return EXIT_FAILURE;
    }
    
    while(mp_cmp(&i2, x) == MP_LT || mp_cmp(&i2, x) == MP_EQ){
        
        if((result = mp_mod(x, &i, &mod_red)) != MP_OKAY){
            printf("Error when x mod i. %s", mp_error_to_string(result));
            return EXIT_FAILURE;
        }
        
        if(mp_cmp_d(&mod_red, 0) == MP_EQ){
            return IS_PRIME_FALSE;
        }

        if((result = mp_add_d(&i, 1, &i)) != MP_OKAY){
            printf("Error when adding 1 to the iterator. %s", mp_error_to_string(result));
            return EXIT_FAILURE;
        }

        if((result = mp_sqr(&i, &i2) != MP_OKAY)){
            printf("Error when calculating the i^2. %s", mp_error_to_string(result));
            return EXIT_FAILURE;
        }
    }

    mp_clear_multi(&i, &i2, &mod_red, NULL);

    return IS_PRIME_TRUE;
}


int get_B_n_G(mp_int *A, mp_int *p, mp_int *B, Point *G){

    mp_int temp1, temp2, temp3, temp4, temp5, temp;
    mp_err err;
    
    if((err = mp_init_multi(&temp1, &temp2, &temp3, &temp4, &temp5, &temp, NULL)) != MP_OKAY){
        printf("Error in initializing at get_B_n_G1. %s", mp_error_to_string(err));
        return EXIT_FAILURE;
    }

    // 4A^3 + 27B^2 ≠ 0 (mod p)
    while(true){
        get_prime(B, 40);
        if((err = mp_expt_n(A, 3, &temp1)) != MP_OKAY){
            printf("Error in A^3. %s", mp_error_to_string(err));
            return EXIT_FAILURE;
        }

        if((err = mp_sqr(B, &temp2)) != MP_OKAY){
            printf("Error in B^2. %s", mp_error_to_string(err));
            return EXIT_FAILURE;
        }

        if((err = mp_mul_d(&temp1, 4, &temp3)) != MP_OKAY){
            printf("Error in 4A. %s", mp_error_to_string(err));
            return EXIT_FAILURE;
        }

        if((err = mp_mul_d(&temp2, 27, &temp4)) != MP_OKAY){
            printf("Error in 27B. %s", mp_error_to_string(err));
            return EXIT_FAILURE;
        }
        
        if((err = mp_add(&temp3, &temp4, &temp5)) != MP_OKAY){
            printf("Error in 4A+27B. %s", mp_error_to_string(err));
            return EXIT_FAILURE;
        }

        if((err = mp_mod(&temp5, p, &temp)) != MP_OKAY){
            printf("Error in 4A+27B (mod p). %s", mp_error_to_string(err));
            return EXIT_FAILURE;
        }

        if(mp_cmp_d(&temp, 0) != MP_EQ){
            break;
        }
    }

    mp_int temp6, temp7, temp8, temp9;

    if((err = mp_init_multi(&temp6, &temp7, &temp8, &temp9, NULL)) != MP_OKAY){
        printf("Error in initializing at get_B_n_G2. %s", mp_error_to_string(err));
        return EXIT_FAILURE;
    }

    // y^2 = x^3 + Ax + B, 
    // randomly generate x of G, 
    // then calculate y by x
    get_prime(&(G->_x), 30);
    if((err = mp_expt_n(&(G->_x), 3, &temp6)) != MP_OKAY){
        printf("Error in x^3. %s", mp_error_to_string(err));
        return EXIT_FAILURE;
    }

    if((err = mp_mul(A, &(G->_x), &temp7)) != MP_OKAY){
        printf("Error in Ax. %s", mp_error_to_string(err));
        return EXIT_FAILURE;
    }

    if((err = mp_add(&temp6, &temp7, &temp8)) != MP_OKAY){
        printf("Error in x^3 + Ax. %s", mp_error_to_string(err));
        return EXIT_FAILURE;
    }

    if((err = mp_add(&temp8, B, &temp9)) != MP_OKAY){
        printf("Error in (...) + B. %s", mp_error_to_string(err));
        return EXIT_FAILURE;
    }

    if((err = mp_sqrt(&temp9, &(G->_y))) != MP_OKAY){
        printf("Error in sqrt(x^3 + Ax + B). %s", mp_error_to_string(err));
        return EXIT_FAILURE;
    }

    mp_clear_multi(&temp1, &temp2, &temp3, &temp4, &temp5, &temp6, &temp7, &temp8, &temp9, &temp, NULL);

    return MP_OKAY;
}

int points_add(Point *p, Point *q, Point *r, const EllipticCurve *ec, bool *zero){

    mp_err err;
    mp_int tmp_zero, tmp1, tmp2, xdiff, ydiff, xdifinv, lambda, tmpy, tmpyinv, sqrx1, tmpsqrx1, tmp3, sqrlmd, tmp4, tmp5, tmp6;
    if((err = mp_init_multi(&tmp_zero, &tmp1, &tmp2, &xdiff, &ydiff, &xdifinv, &lambda, &tmpy, &tmpyinv, &sqrx1, &tmpsqrx1, &tmp3, &sqrlmd, &tmp4, &tmp5, &tmp6, NULL)) != MP_OKAY){
        printf("Error in initializing at points addition. %s", mp_error_to_string(err));
        return EXIT_FAILURE;
    }

    mp_int t;
    if((err = mp_init_set(&t, 2)) != MP_OKAY){
        printf("Error in initializing constant 2. %s", mp_error_to_string(err));
        return EXIT_FAILURE;
    }

    if(*zero){
        // Q = 0,
        // P -> R
        if((err = mp_copy(&(p->_x), &(r->_x))) != MP_OKAY){
            printf("Error in copying px to rx. %s", mp_error_to_string(err));
            return EXIT_FAILURE;
        }
        if((err = mp_copy(&(p->_y), &(r->_y))) != MP_OKAY){
            printf("Error in copying py to ry. %s", mp_error_to_string(err));
            return EXIT_FAILURE;
        }
    }else{
        mp_zero(&tmp_zero);

        if((err = mp_sub(&(q->_x), &(p->_x), &xdiff)) != MP_OKAY){
            printf("Error in x2 - x1. %s", mp_error_to_string(err));
            return EXIT_FAILURE;
        }

        if(mp_cmp(&xdiff, &tmp_zero) == MP_LT){
            if((err = mp_add(&xdiff, &(ec->_p), &tmp1)) != MP_OKAY){
                printf("Error in xdiff + p. %s", mp_error_to_string(err));
                return EXIT_FAILURE;
            }
            mp_zero(&xdiff);
            if((err = mp_copy(&tmp1, &xdiff)) != MP_OKAY){
                printf("Error in copy xdiff+p -> xdiff. %s", mp_error_to_string(err));
                return EXIT_FAILURE;
            }
        }

        if((err = mp_sub(&(q->_y), &(p->_y), &ydiff)) != MP_OKAY){
            printf("Error in y2 - y1. %s", mp_error_to_string(err));
            return EXIT_FAILURE;
        }

        if(mp_cmp(&ydiff, &tmp_zero) == MP_LT){
            if((err = mp_add(&ydiff, &(ec->_p), &tmp2)) != MP_OKAY){
                printf("Error in ydiff + p. %s", mp_error_to_string(err));
                return EXIT_FAILURE;
            }
            mp_zero(&ydiff);
            if((err = mp_copy(&tmp2, &ydiff)) != MP_OKAY){
                printf("Error in copy ydiff+p -> ydiff. %s", mp_error_to_string(err));
                return EXIT_FAILURE;
            }
        }

        if(mp_cmp(&xdiff, &tmp_zero) != MP_EQ){
            // Px != Qx
            // λ = (Yq - Yp)/(Xq - Xp) mod p
            if((err = mp_invmod(&xdiff, &(ec->_p), &xdifinv)) != MP_OKAY){
                printf("Error in xdiff^-1 mod p. %s", mp_error_to_string(err));
                return EXIT_FAILURE;
            }

            if((err = mp_mulmod(&ydiff, &xdifinv, &(ec->_p), &lambda)) != MP_OKAY){
                printf("Error in ydiff/xdiff mod p. %s", mp_error_to_string(err));
                return EXIT_FAILURE;
            }
        }else{
            // Px=Qx

            if(mp_cmp(&ydiff, &tmp_zero) == MP_EQ){
                // Py = Qy, 
                // => P = Q
                // λ = (3Xp² + a) / 2Yp mod p

                if((err = mp_mulmod(&t, &(p->_y), &(ec->_p), &tmpy)) != MP_OKAY){
                    printf("Error in 2*y1 mod p. %s", mp_error_to_string(err));
                    return EXIT_FAILURE;
                }

                if((err = mp_invmod(&tmpy, &(ec->_p), &tmpyinv)) != MP_OKAY){
                    printf("Error in tmpy1^-1 mod p. %s", mp_error_to_string(err));
                    return EXIT_FAILURE;
                }

                if((err = mp_sqr(&(p->_x), &sqrx1)) != MP_OKAY){
                    printf("Error in sqr x1. %s", mp_error_to_string(err));
                    return EXIT_FAILURE;
                }

                if((err = mp_mul_d(&sqrx1, 3, &tmpsqrx1)) != MP_OKAY){
                    printf("Error in 3 * sqrx1. %s", mp_error_to_string(err));
                    return EXIT_FAILURE;
                }

                if((err = mp_add(&tmpsqrx1, &(ec->_a), &tmp3)) != MP_OKAY){
                    printf("Error in (3 * sqrx1) + A. %s", mp_error_to_string(err));
                    return EXIT_FAILURE;
                }

                if((err = mp_mulmod(&tmp3, &tmpyinv, &(ec->_p), &lambda)) != MP_OKAY){
                    printf("Error in ((3*sqrx1)+A)/2y1 mod p. %s", mp_error_to_string(err));
                    return EXIT_FAILURE;
                }
            }else{
                // Px = Qx, Py != Qy
                // => vertical, null->R
                *zero = true;
            }
        }

        if(!(*zero)){
            if((err = mp_sqr(&lambda, &sqrlmd)) != MP_OKAY){
                printf("Error in sqr lambda. %s", mp_error_to_string(err));
                return EXIT_FAILURE;
            }

            if((err = mp_sub(&sqrlmd, &(p->_x), &tmp4)) != MP_OKAY){
                printf("Error in sqrlmd - x1. %s", mp_error_to_string(err));
                return EXIT_FAILURE;
            }

            if((err = mp_submod(&tmp4, &(q->_x), &(ec->_p), &(r->_x))) != MP_OKAY){
                printf("Error in sqrlmd-x1-x2 mod p. %s", mp_error_to_string(err));
                return EXIT_FAILURE;
            }

            if((err = mp_sub(&(p->_x), &(r->_x), &tmp5)) != MP_OKAY){
                printf("Error in x1 - x3. %s", mp_error_to_string(err));
                return EXIT_FAILURE;
            }

            if((err = mp_mul(&tmp5, &lambda, &tmp6)) != MP_OKAY){
                printf("Error in lmd*(x1 - x3). %s", mp_error_to_string(err));
                return EXIT_FAILURE;
            }

            if((err = mp_submod(&tmp6, &(p->_y), &(ec->_p), &(r->_y))) != MP_OKAY){
                printf("Error in (lmd*(x1-x3)-y1) mod p. %s", mp_error_to_string(err));
                return EXIT_FAILURE;
            }
        }

    }

    mp_clear_multi(&tmp_zero, &tmp1, &tmp2, &xdiff, &ydiff, &xdifinv, &lambda, &tmpy, &tmpyinv, &sqrx1, &tmpsqrx1, &tmp3, &sqrlmd, &tmp4, &tmp5, &tmp6, &t, NULL);

    return MP_OKAY;

}

int point_copy(Point *a, Point *b){
    mp_err err;
    if((err = mp_copy(&(a->_x), &(b->_x))) != MP_OKAY){
        printf("Error in copy ax to bx. %s", mp_error_to_string(err));
        return EXIT_FAILURE;
    }
    if((err = mp_copy(&(a->_y), &(b->_y))) != MP_OKAY){
        printf("Error in copy ay to by. %s", mp_error_to_string(err));
        return EXIT_FAILURE;
    }

    return MP_OKAY;
}

int point_init(Point *a){
    mp_err err;
    if((err = mp_init(&(a->_x))) != MP_OKAY){
        printf("Error in init x of point. %s", mp_error_to_string(err));
        return EXIT_FAILURE;
    }
    if((err = mp_init(&(a->_y))) != MP_OKAY){
        printf("Error in init y of point. %s", mp_error_to_string(err));
        return EXIT_FAILURE;
    }

    return MP_OKAY;
}

int point_init_set(Point *a, mp_digit x, mp_digit y){
    mp_err err;
    if((err = mp_init_set(&(a->_x), x)) != MP_OKAY){
        printf("Error in init x3. %s", mp_error_to_string(err));
        return EXIT_FAILURE;
    }
    if((err = mp_init_set(&(a->_y), y)) != MP_OKAY){
        printf("Error in init y3. %s", mp_error_to_string(err));
        return EXIT_FAILURE;
    }
    return MP_OKAY;
}

int point_init_copy(Point *a, Point *b){
    point_init(a);
    point_copy(b, a);
}

int scalar_mul(Point *pub_K, Point *G, const EllipticCurve *ec, mp_int *d){
    mp_int x1, x2, x3, y1, y2, y3;
    mp_int A, p;

    Point p1, p2, p3;
    point_init_copy(&p1, G);
    point_init_copy(&p2, G);
    point_init_set(&p3, 0, 0);


    char byte_arr[BYTE_ARRAY_MAX_LEN] = {0};

    mp_err err;

    if((err = mp_to_radix(d, byte_arr, BYTE_ARRAY_MAX_LEN, NULL, 2)) != MP_OKAY){
        printf("Error in d to byte array. %s", mp_error_to_string(err));
        return EXIT_FAILURE;
    }

    if((err = mp_init_set(&x3, 0)) != MP_OKAY){
        printf("Error in init x3. %s", mp_error_to_string(err));
        return EXIT_FAILURE;
    }
    if((err = mp_init_set(&y3, 0)) != MP_OKAY){
        printf("Error in init y3. %s", mp_error_to_string(err));
        return EXIT_FAILURE;
    }

    if((err = mp_init_copy(&x1, &(G->_x))) != MP_OKAY){
        printf("Error in init x1. %s", mp_error_to_string(err));
        return EXIT_FAILURE;
    }
    if((err = mp_init_copy(&x2, &(G->_x))) != MP_OKAY){
        printf("Error in init x2. %s", mp_error_to_string(err));
        return EXIT_FAILURE;
    }

    if((err = mp_init_copy(&y1, &(G->_y))) != MP_OKAY){
        printf("Error in init y1. %s", mp_error_to_string(err));
        return EXIT_FAILURE;
    }
    if((err = mp_init_copy(&y2, &(G->_y))) != MP_OKAY){
        printf("Error in init y2. %s", mp_error_to_string(err));
        return EXIT_FAILURE;
    }

    /**
     * @brief 倍乘法
     * ex. 141P = 2^7P + 2^3P + 2^2P + 2^0P
     * 2P=P+P, 4P=2P+2P, 8P=4P+4P, ..., 128P=64P+64P
     * 2^7P + 2^3P + 2^2P + 2^0P
     */
    char c = '1';
    for(int i = 0; i < PRIVATE_KEY_LEN; i++){
        if(byte_arr[i] == c){
            points_add()
        }
    }
}