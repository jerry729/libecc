# 齐次坐标



# Schoof's algorithm



# 椭圆曲线原理

## 椭圆曲线
适合加密的椭圆曲线 y^2 = x^3 + ax + b, 4a^3 + 27b^2 != 0, 0为无穷远点
椭圆曲线上的群

## 有限域
加法，二倍运算，逆运算
有限域 y^2 同余 x^3 + ax + b (mod p), 离散化成GF(p)的椭圆曲线

## 子群
标量积 xG 2G二倍运算, 3G = 2G + G加法运算
阶
基准点的向量倍数 子群

## 子群的阶
P是基准点
子群的阶就是子群中的元素个数，等价于 nP=0，其n是正整数，且n是其中最小的一个，这个n就是子群的阶
根据拉格朗日定理，子群P的阶，是父群阶的一个因子，比如一个椭圆曲线有N个元素，它的一个子群有n个元素，n能整除N

接着我们开始算子群P的阶：

用Schoof's算法计算椭圆曲线的阶N
找出N的所有因子
对N的每一个因子n，计算nP
找出使nP=0的最小的n，这个n就是子群的阶

## 阶N找较大因子n，根据n确定基准点
也就是说，我们先找子群的阶，再根据子群的阶找其中的基准点。而不是先找基准点，再算子群的阶。

再加一个信息：根据拉格朗日定理，h=N/n肯定是个整数（因为n是N的一个因子）。这个h叫做子群的协因子(cofactor).

对一个椭圆曲线上的**任意一点**来说， NP=0,因为N总是任何一个n的倍数；
也就是 n(hP)=0
取n是素数，点G=hP生成一个阶为n的子群。
这样，我们就根据选择的阶n，算出了一个基准点



# libtommath

## Data type

```c
typedef struct { 
  int used, alloc; 
  mp_sign sign; 
  mp_digit *dp; //28-bits unsigned long, 32-28 reamining bits for carrier/borrower
} mp_int;
```



## Error

| code    | meaning                                     |
| ------- | ------------------------------------------- |
| MP_OKAY | The function succeeded.                     |
| MP_VAL  | The function input was invalid.             |
| MP_MEM  | Heap memory exhausted.                      |
| MP_ITER | Maximum iterations reached.                 |
| MP_BUF  | Buffer overflow, supplied buffer too small. |

`char *mp_error_to_string(mp_err code);`



## Arithmetic functions

```c
mp_add(&a, &b, &c); /* c = a + b */
mp_mul(&a, &a, &c); /* c = a * a */
mp_div(&a, &b, &c, &d); /* c = [a/b], d = a mod b */
```

```c
mp_add(&a, &b, &b);       /* b = a + b */
mp_div(&a, &b, &a, &c);   /* a = [a/b], c = a mod b */
```



## Initialization

- single 
    ```c
    mp_err mp_init (mp_int *a);
    void mp_clear (mp_int *a);
    ```
- multiple
    ```c
    mp_err mp_init_multi(mp_int *mp, ...);
    void mp_clear_multi(mp_int *mp, ...);
    ```
- other
    ```c
    mp_err mp_init_copy (mp_int *a, mp_int *b);
    mp_err mp_init_size (mp_int *a, int size);
    ```



## Maintainance Functions

- clear leading zeros

  `void mp_clamp(mp_int *a);`

- zero out

  `void mp_zero(mp_int *a);`

- reducing memory usage

  `mp_err mp_shrink (mp_int *a);`

- adding additional digits

  `mp_err mp_grow (mp_int *a, int size);`



## Basic operations

- copying

  `mp_err mp_copy (const mp_int *a, mp_int *b);`

  `void mp_exch (mp_int *a, mp_int *b); // swap pointer with a temporary pointer`

- bit counting

  `int mp_cnt_lsb(const mp_int *a); // LSB, the Lowest Significant Bit `

  `int mp_count_bits(const mp_int *a); // MSB, the Most Significant Bit `

- small constants

  - single digit

    `void mp_set (mp_int *a, mp_digit b);`

  - Int32 and Int64 Constants

    ```c
    void mp_set_i32 (mp_int *a, int32_t b);
    void mp_set_u32 (mp_int *a, uint32_t b);
    void mp_set_i64 (mp_int *a, int64_t b);
    void mp_set_u64 (mp_int *a, uint64_t b);
    ```

    ```c
    //obtain these values again
    int32_t mp_get_i32 (const mp_int *a);
    uint32_t mp_get_u32 (const mp_int *a);
    uint32_t mp_get_mag_u32 (const mp_int *a);
    int64_t mp_get_i64 (const mp_int *a);
    uint64_t mp_get_u64 (const mp_int *a);
    uint64_t mp_get_mag_u64 (const mp_int *a); // mag* absolute value
    ```

  - long constants \- platform dependant

    ```c
    void mp_set_l (mp_int *a, long b);
    void mp_set_ul (mp_int *a, unsigned long b);
    ```

    ```c
    long mp_get_l (const mp_int *a);
    unsigned long mp_get_ul (const mp_int *a);
    unsigned long mp_get_mag_ul (const mp_int *a);
    ```

  - Floating Point Constants - platform dependant

    `mp_err mp_set_double(mp_int *a, double b); // assign interger part of b to a``

    ``double mp_get_double(const mp_int *a);`

  - Initialize and Setting Constants

    ```c
    mp_err mp_init_set (mp_int *a, mp_digit b);
    mp_err mp_init_i32 (mp_int *a, int32_t b);
    mp_err mp_init_u32 (mp_int *a, uint32_t b);
    mp_err mp_init_i64 (mp_int *a, int64_t b);
    mp_err mp_init_u64 (mp_int *a, uint64_t b);
    mp_err mp_init_l   (mp_int *a, long b);
    mp_err mp_init_ul  (mp_int *a, unsigned long b);
    ```

- comparisons

  | Result code | Meaning |
  | ----------- | ------- |
  | MP_GT       | a > b   |
  | MP_EQ       | a = b   |
  | MP_LT       | a < b   |

  - unsigned comparison

    `mp_ord mp_cmp_mag(mp_int *a, mp_int *b);`

  - Signed comparison

    `mp_ord mp_cmp(mp_int *a, mp_int *b);`

  - Single Digit

    `mp_ord mp_cmp_d(mp_int *a, mp_digit b);`

- Logical Operations

  - Multiplication by two

    ```c
    mp_err mp_mul_2(const mp_int *a, mp_int *b);
    mp_err mp_div_2(const mp_int *a, mp_int *b);
    ```

    `mp_err mp_mul_2d(const mp_int *a, int b, mp_int *c); //  a * 2^b -> c`

    `mp_err mp_div_2d (const mp_int *a, int b, mp_int *c, mp_int *d); // a/2^b -> c...d`

    `mp_err mp_2expt(mp_int *a, int b); // the power of two 2^b`

  - Polynomial Basis Operations

    `mp_err mp_lshd (mp_int *a, int b);`

    `void mp_rshd (mp_int *a, int b)`

  - AND, OR, XOR and COMPLEMENT Operations

    ```c
    mp_err mp_or  (const mp_int *a, mp_int *b, mp_int *c);
    mp_err mp_and (const mp_int *a, mp_int *b, mp_int *c);
    mp_err mp_xor (const mp_int *a, mp_int *b, mp_int *c);
    mp_err mp_complement(const mp_int *a, mp_int *b);
    mp_err mp_signed_rsh(const mp_int *a, int b, mp_int *c, mp_int *d);
    ```

- Addition and Subtraction

  ```c
  mp_err mp_add (const mp_int *a, const mp_int *b, mp_int *c);
  mp_err mp_sub (const mp_int *a, const mp_int *b, mp_int *c)
  ```

- Sign Manipulation

  - Negation

    `mp_err mp_neg (const mp_int *a, mp_int *b);`

  - Absolute

    `mp_err mp_abs (const mp_int *a, mp_int *b);`

- Integer Division and Remainder

  `mp_err mp_div (const mp_int *a, const mp_int *b, mp_int *c, mp_int *d);`

- Hashing

  `mp_err mp_hash (mp_int *a, mp_hval *hash);`



## Multiplication and Squaring

- Multiplication

  `mp_err mp_mul (const mp_int *a, const mp_int *b, mp_int *c);`



