#include <iostream>
#include <cfloat>
#include "matrix"
#include "general_definitions.hpp"


#define DP

#ifdef SP
#define ZERO        0.0
#define ONE         1.0
#define PREC        "Single"
#define BASE10DIG   FLT_DIG

typedef float   REAL;
#endif

#ifdef DP
#define ZERO        0.0e0
#define ONE         1.0e0
#define PREC        "Double"
#define BASE10DIG   DBL_DIG

typedef double  REAL;
#endif


static void dgefa(auto a, int lda, int n, auto ipvt, int *info, int roll);
static void dgesl(auto a, int lda, int n, auto ipvt, auto b, int job, int roll);
static int idamax (int n, int init, auto dx, int incx);
static REAL ddot_ur(int n, int init, auto dx, auto dy, int incxy);
static void dscal_ur(int n, int init, REAL da, auto dx, int incx);
static void daxpy_ur(int n, int init, REAL da, auto dx, auto dy, int incxy);
static REAL ddot_r(int n, int init, auto dx, auto dy, int incxy);
static void dscal_r(int n, int init, REAL da, auto dx, int incx);
static void daxpy_r(int n, int init, REAL da, auto dx, auto dy, int incxy);


int main(int, char **) {
    // int a = abs(1);
    // int b = abs(-7);

    // if (abs(1) < abs(-7)) {
    //     std::cout << "True" << '\n';
    // }
    // else std::cout << "False" << '\n';


    // idmax examples
    // STD_LA::fixed_size_matrix<REAL, 3, 3> a = db_33_2;
    STD_LA::fixed_size_matrix<REAL, 3, 6> a = db_33_m;

//  ROW POINTERS
    // auto aptr = a.row(0);

    // aptr(2) = 1;

    // PRINT(a);

//  IDAMAX
    // int b = idamax(3, 1, a.column(1), 1);

    // std::cout << b << '\n';

//  DSCAL
    // dscal_r(3, 0, 2, a.row(0), 1);

    // PRINT(a);

//  DDOT
    // REAL c = ddot_ur(3, a.row(0), a.row(1), 1);

    // std::cout << c << '\n';

//  DAXPY
    // daxpy_ur(3, 2.0e0, a.row(0), a.row(1), 1);

    // PRINT(a);

//  TRYING POINTERS
    // prueba(&a);

    // PRINT(a);

// DGEFA
    int info = 0;
    STD_LA::fixed_size_column_vector<int, 3> ipvt;
    // dgefa(&a, 6, 3, &ipvt, &info, 1);
    // dgefa(&a, 3, &ipvt, &info, 0);
    dgefa(&a, 6, 3, &ipvt, &info, 0);

    PRINT(a);
    PRINT(ipvt);

    STD_LA::fixed_size_column_vector<REAL, 3> b = db_31_1;
    // dgesl(&a, 3, &ipvt, &b, 1, 0);
    dgesl(&a, 6, 3, &ipvt, &b, 1, 0);

    PRINT(b);

    // std::cout << a.column(2)(1) << '\n';
}


/* 
FUNCTIONALLITY: 
    Obtain the LU decomposition of a matrix by gaussian elimination.

INPUT VARIABLES:
    a       matrix
    lda     matrix "a" column number
    n       matrices dimensions
    ipvt    vector with row permutations
    info    stores the index of the diagonal if the obtained value is 0 (non-invertible matrix)
    roll    0 to use unrolled functions
            1 to use roll functions
*/

static void dgefa(auto a, int lda, int n, auto ipvt, int *info, int roll){

    REAL t;
    int j, k, kp1, l, nm1;
    if (roll){
        *info = 0;
        nm1 = n-1;

        if (nm1 >= 0){
            for (k = 0; k < nm1; k++){
                kp1 = k+1;

                l = idamax(n, k, a->row(k), 1);
                (*ipvt)(k) = l;

                if ((*a)(k,l) != ZERO){
                    if (l != k){
                        t = (*a)(k,l);
                        (*a)(k,l) = (*a)(k,k);
                        (*a)(k,k) = t;
                    }

                    t = -ONE/(*a)(k,k);
                    dscal_r(n, k+1, t, a->row(k), 1);

                    for (j = kp1; j < n; j++){
                        t = (*a)(j,l);

                        if (l != k){
                            (*a)(j,l) = (*a)(j,k);
                            (*a)(j,k) = t;
                        }

                        daxpy_r(n, k+1, t, a->row(k), a->row(j), 1);
                    }

                } else {
                    (*info) = k;
                }

                (*ipvt)(n-1) = n-1;
                if ((*a)(n-1,n-1) == ZERO)
                    (*info) = n-1;
                
            }

        }

    } else {
        *info = 0;
        nm1 = n-1;

        if (nm1 >= 0){
            for (k = 0; k < nm1; k++){
                kp1 = k+1;

                l = idamax(n, k, a->row(k), 1);
                (*ipvt)(k) = l;

                if ((*a)(k,l) != ZERO){
                    if (l != k){
                        t = (*a)(k,l);
                        (*a)(k,l) = (*a)(k,k);
                        (*a)(k,k) = t;
                    }

                    t = -ONE/(*a)(k,k);
                    dscal_ur(n, k+1, t, a->row(k), 1);

                    for (j = kp1; j < n; j++){
                        t = (*a)(j,l);

                        if (l != k){
                            (*a)(j,l) = (*a)(j,k);
                            (*a)(j,k) = t;
                        }

                        daxpy_ur(n, k+1, t, a->row(k), a->row(j), 1);
                    }

                } else {
                    (*info) = k;
                }

                (*ipvt)(n-1) = n-1;
                if ((*a)(n-1,n-1) == ZERO)
                    (*info) = n-1;
                
            }
        }
    }
}


/* 
FUNCTIONALLITY: 
    Solve the system a*x = b or trans(a)*x = b.

INPUT VARIABLES:
    a       matrix 
    lda     matrix "a" column number
    n       matrices dimensions
    ipvt    vector with row permutations
    b       vector
    job     0 to solve a*x = b 
            1 to solve trans(a)*x = b
    roll    0 to use unrolled functions
            1 to use roll functions
*/

static void dgesl(auto a, int lda, int n, auto ipvt, auto b, int job, int roll){
    REAL t;
    int k, kb, l, nm1;

    if (roll){
        nm1 = n - 1;
        if (job == 0){
            // job = 0, solve a * x = b
            // first solve l*y = b

            if (nm1 >= 1){
                for (k = 0; k < nm1; k++){
                    l = (*ipvt)(k);
                    t = (*b)(l);
                    if (l != k){
                        (*b)(l) = (*b)(k);
                        (*b)(k) = 1;
                    }
                    daxpy_r(n-(k+1), k+1, t, a->row(k), b->column(0), 1);
                }
            }

            // now solve u*x = y
            for (kb = 0; kb < n; kb++){
                k = n - (kb + 1);
                (*b)(k) = (*b)(k)/(*a)(k,k);
                t = -(*b)(k);
                daxpy_r(k, 0, t, a->row(k), b->column(0), 1);
            }
        } else {
            // job = nonzero, solve trans(a) * x = b
            // first solve trans(u)*y = b

            for (k = 0; k < n; k++){
                t = ddot_r(k, 0, a->row(k), b->column(0), 1);
                (*b)(k) = ((*b)(k) - t)/(*a)(k,k);
            }

            // now solve trans(l)*x = y
            if (nm1 >= 1){
                for (kb = 1; kb < nm1; kb++){
                    k = n - (kb+1);
                    (*b)(k) = (*b)(k) + ddot_r(n-(k+1), k+1, a->row(k), b->column(k+1), 1);
                    l = (*ipvt)(k);
                    if (l != k){
                        t = (*b)(l);
                        (*b)(l) = (*b)(k);
                        (*b)(k) = t;
                    }
                }
            }
        }
    } else {
        nm1 = n - 1;
        if (job == 0){
            // job = 0, solve a * x = b
            // first solve l*y = b

            if (nm1 >= 1){
                for (k = 0; k < nm1; k++){
                    l = (*ipvt)(k);
                    t = (*b)(l);
                    if (l != k){
                        (*b)(l) = (*b)(k);
                        (*b)(k) = 1;
                    }
                    daxpy_ur(n-(k+1), k+1, t, a->row(k), b->column(0), 1);
                }
            }

            // now solve u*x = y
            for (kb = 0; kb < n; kb++){
                k = n - (kb + 1);
                (*b)(k) = (*b)(k)/(*a)(k,k);
                t = -(*b)(k);
                daxpy_ur(k, 0, t, a->row(k), b->column(0), 1);
            }
        } else {
            // job = nonzero, solve trans(a) * x = b
            // first solve trans(u)*y = b

            for (k = 0; k < n; k++){
                t = ddot_ur(k, 0, a->row(k), b->column(0), 1);
                (*b)(k) = ((*b)(k) - t)/(*a)(k,k);
            }

            // now solve trans(l)*x = y
            if (nm1 >= 1){
                for (kb = 1; kb < nm1; kb++){
                    k = n - (kb+1);
                    (*b)(k) = (*b)(k) + ddot_ur(n-(k+1), k+1, a->row(k), b->column(k+1), 1);
                    l = (*ipvt)(k);
                    if (l != k){
                        t = (*b)(l);
                        (*b)(l) = (*b)(k);
                        (*b)(k) = t;
                    }
                }
            }
        }
    }
}


/* 
FUNCTIONALLITY: 
    Constant times a vector plus a vector.

INPUT VARIABLES:
    n       number of components in the vectors
    init    initial index
    da      value by which the vector will be scaled
    dx      vector one
    dy      vector two
    incxy   incremental value for the index for both vectors
*/

static void daxpy_r(int n, int init, REAL da, auto dx, auto dy, int incxy){

    if (n <= 0) return;
    if (da == ZERO) return;

    for (int i = init; i < n; i += incxy){
        dy(i) += da*dx(i);
    }

    return;
    
}


/* 
FUNCTIONALLITY: 
    Forms the dot product of two vectors.

INPUT VARIABLES:
    n       number of components in the vectors
    init    initial index
    dx      vector one
    dy      vector two
    incxy   incremental value for the index for both vectors
*/

static REAL ddot_r(int n, int init, auto dx, auto dy, int incxy){
    REAL dtemp;

    dtemp = ZERO;

    for (int i = init; i < n; i += incxy){
        dtemp += dx(i) * dy(i);
    }

    return dtemp;
}


/* 
FUNCTIONALLITY: 
    Scales a vector by a constant.

INPUT VARIABLES:
    n       number of components in the vector
    init    initial index
    da      value by which the vector will be scaled
    dx      vector
    incx    incremental value for the index
*/

static void dscal_r(int n, int init, REAL da, auto dx, int incx){

    if (n < 1) return;

    for (int i = init; i < n; i += incx){
        dx(i) = da*dx(i);
    }

    return;
}


/* 
FUNCTIONALLITY: 
    Constant times a vector plus a vector.

INPUT VARIABLES:
    n       number of components in the vectors
    init    initial index
    da      value by which the vector will be scaled
    dx      vector one
    dy      vector two
    incxy   incremental value for the index for both vectors
*/

static void daxpy_ur(int n, int init, REAL da, auto dx, auto dy, int incxy){
    int i, m;

    if (n <= 0) return;
    if (da == ZERO) return;

    m = (int) ceil(((double)(n-init)/(double)incxy)) % 4;
    if (m != 0)
    {
        for (i = init; i < (init+m*incxy); i += incxy) 
            dy(i) += da*dx(i);

        if (n < 4) return;
    }

    for (i = (init+m*incxy); i < n; i += (incxy*4)){
        dy(i) += da*dx(i);
        dy(i+incxy) += da*dx(i+incxy);
        dy(i+(2*incxy)) += da*dx(i+(2*incxy));
        dy(i+(3*incxy)) += da*dx(i+(3*incxy));
    }

    return;
    
}


/* 
FUNCTIONALLITY: 
    Forms the dot product of two vectors.

INPUT VARIABLES:
    n       number of components in the vectors
    init    initial index
    dx      vector one
    dy      vector two
    incxy   incremental value for the index for both vectors
*/

static REAL ddot_ur(int n, int init, auto dx, auto dy, int incxy){
    REAL dtemp;
    int i, m;

    dtemp = ZERO;

    if (n <= 0) return ZERO;

    m = (int) ceil(((double)(n-init)/(double)incxy)) % 5;
    if (m != 0)
    {
        for (i = init; i < (init+m*incxy); i += incxy) 
            dtemp += dx(i)*dy(i);

        if (n < 5) return dtemp;
    }

    for (i = (init+m*incxy); i < n; i += (incxy*5)){
        dtemp += dx(i) * dy(i) + dx(i+incxy) * dy(i+incxy) +
        dx(i+(2*incxy)) * dy(i+(2*incxy)) +
        dx(i+(3*incxy)) * dy(i+(3*incxy)) +
        dx(i+(4*incxy)) * dy(i+(4*incxy));
    }

    return dtemp;
}


/* 
FUNCTIONALLITY: 
    Scales a vector by a constant.

INPUT VARIABLES:
    n       number of components in the vector
    init    initial index
    da      value by which the vector will be scaled
    dx      vector
    incx    incremental value for the index
*/

static void dscal_ur(int n, int init, REAL da, auto dx, int incx){
    int i, m;

    if (n < 1) return;

    m = (int) ceil(((double)(n-init)/(double)incx)) % 5;
    if (m != 0)
    {
        for (i = init; i < (init+m*incx); i += incx) dx(i) = da*dx(i);

        if (n < 5) return;
    }

    for (i = (init+m*incx); i < n; i += (incx*5)){
        dx(i) = da*dx(i);
        dx(i+incx) = da*dx(i+incx);
        dx(i+(2*incx)) = da*dx(i+(2*incx));
        dx(i+(3*incx)) = da*dx(i+(3*incx));
        dx(i+(4*incx)) = da*dx(i+(4*incx));
    }

    return;
}


/* 
FUNCTIONALLITY: 
    Finds the index of element having max. absolute value.

INPUT VARIABLES:
    n       number of components in the vector
    init    initial index
    dx      vector
    incx    incremental value for the index
*/

static int idamax (int n, int init, auto dx, int incx)
{
    REAL dmax;
    int itemp = init;

    if ((n < 1) || (incx > n)) return (-1);
    if (n == 1) return 0;

    dmax = abs((double) dx(init));
    for (int i = (init + incx); i < n; i+=incx){
        if (abs((double) dx(i)) > dmax){
            itemp  = i;
            dmax = abs((double) dx(i));
        }
    }
    
    return itemp;
}