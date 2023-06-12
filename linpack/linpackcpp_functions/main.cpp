#include <iostream>
#include <chrono>
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

#define nn 10


static REAL linpack (long nreps, int job, int roll, int n, char *matOpt);
static void matgen(auto a, auto ipvt, auto b, int n);
static void dgefa(auto a, int n, auto ipvt, auto info, int roll);
static void dgesl(auto a, int n, auto ipvt, auto b, int job, int roll);
static int idamax (int n, int init, auto dx, int incx);
static REAL ddot_ur(int n, int init, auto dx, auto dy, int incxy);
static void dscal_ur(int n, int init, REAL da, auto dx, int incx);
static void daxpy_ur(int n, int init, REAL da, auto dx, auto dy, int incxy);
static REAL ddot_r(int n, int init, auto dx, auto dy, int incxy);
static void dscal_r(int n, int init, REAL da, auto dx, int incx);
static void daxpy_r(int n, int init, REAL da, auto dx, auto dy, int incxy);

using namespace std::chrono_literals;

int main (void){
    long nreps;
    char buf[5], opt[5];
    int n;

    while(1){
        // MENU
        std::cout << "Introduce one of the following options: \n\n"         \
                  << "a - Dinamic matrix\n"                                 \
                  << "b - Fixed size matrix (" << nn << "x" << nn << ")\n"    \
                  << "c - General matrix (" << nn << "x" << nn << ")\n"       \
                  << "q - QUIT\n\n";
        std::cin >> opt;

        if (opt[0] == 'q' || opt[0] == 'Q') break;

        if (opt[0] == 'a' || opt[0] == 'A'){
            std::cout << "\n\nIntroduce matrix size: ";
            std:: cin >> buf;

            if (buf[0]=='\0' || buf[0]=='\n') n = nn;
            else n = atoi(buf);
        } else {
            n = nn;
        }

        std::cout << "\n\nLINPACK benchmark, " << PREC <<" precision. \n";
        std::cout << "Machine precision: " << BASE10DIG << " digits. \n";
        std::cout << "Array size " << n << " X " << n << ". \n";
        std::cout << "Average rolled and unrolled performance:\n\n";
        std::cout << "    Reps Time(s) DGEFA   DGESL  OVERHEAD    KFLOPS\n";
        std::cout << "----------------------------------------------------\n";

        nreps = 1;
        while (linpack(nreps,0, 0, n, opt)<10.){
            nreps *= 2;
        }
        
        std::cout << '\n';
    }
}


/* 
FUNCTIONALLITY: 
    Measure execution times.

INPUT VARIABLES:
    nreps     number of times to repeat the computations
    job       0 to solve a*x = b 
              1 to solve trans(a)*x = b
    roll      0 to use unrolled functions
              1 to use roll functions
    n         matrices dimensions
    matOpt    matrix type to use
*/

static REAL linpack (long nreps, int job, int roll, int n, char *matOpt){

    REAL kflops, toverhead, ops, tdgefa = 0.0, tdgesl = 0.0, totalt = 0.0;
    int info = 0;

    STD_LA::general_matrix<REAL, nn, 2*nn> a;
    STD_LA::general_matrix<int, nn, nn> ipvt;
    STD_LA::general_column_vector<REAL, nn> b;

    // if (matOpt[0] == 'a' || matOpt[0] == 'A'){
    //     STD_LA::dynamic_matrix<REAL> a;
    //     STD_LA::dynamic_matrix<int> ipvt;
    //     STD_LA::dynamic_column_vector<REAL> b;

    //     a.resize(n,2*n);
    //     ipvt.resize(n,n);
    //     b.resize_rows(n);
    // }

    // if (matOpt[0] == 'b' || matOpt[0] == 'B'){
    //     STD_LA::fixed_size_matrix<REAL, nn, 2*nn> a;
    //     STD_LA::fixed_size_matrix<int, nn, nn> ipvt;
    //     STD_LA::fixed_size_column_vector<REAL, nn> b;
    // }


    ops = ((2.0*n*n*n)/3.0+2.0*n*n);

    using clk = std::chrono::high_resolution_clock;
    using secs = std::chrono::seconds;
    auto t0 = clk::now();
    for (int i = 0; i < nreps; ++i){
        matgen(&a, &ipvt, &b, n);
        auto t1 = clk::now();
        dgefa(&a, n, &ipvt, &info, roll);
        auto t2 = clk::now();
        tdgefa += (t2-t1)/1.s;
        auto t3 = clk::now();
        dgesl(&a, n, &ipvt, &b, job, roll);
        auto t4 = clk::now();
        tdgesl += (t4-t3)/1.s;
    }
    auto t5 = clk::now();
    totalt += (t5-t0)/1.s;

    if (totalt<0.5 || tdgefa+tdgesl<0.2) return (0.);
    kflops = 2.* nreps*ops/(1000. * (tdgefa+tdgesl));
    toverhead = totalt-tdgesl;
    if (tdgefa<0.) tdgefa=0.;
    if (tdgesl<0.) tdgesl=0.;
    if (toverhead<0.) toverhead=0.;

    std::cout << nreps << " " << totalt << " " << 100.*tdgefa/totalt << " "         \
    << 100.*tdgesl/totalt << " " << 100.*toverhead/totalt << " " << kflops << '\n';
    
    return totalt;
}


/* 
FUNCTIONALLITY: 
    Insert values in the matrix.

INPUT VARIABLES:
    a       matrix (leftside upper, rightside lower)
    ipvt    matrix with row permutations
    b       vector
    n       matrices dimensions
*/

static void matgen(auto a, auto ipvt, auto b, int n){
    int init, i, j;

    init = 1325;
    for (i = 0; i < n; ++i){
        for (j = 0; j < n; ++j){
            init = (int) ((long) 3125 * (long)init % 65536L);
            (*a)(i,j) = (init - 32768.0)/16384.0;
            (*a)(i,j+n) = ZERO;
            (*ipvt)(i,j) = 0;
        }
    }

    for (i = 0; i < n; ++i){
        (*b)(i) = 0.0;
        (*ipvt)(i,i) = 1;
        for (j = 0; j < n; ++j){
            (*b)(i) = (*b)(i) + (*a)(i,j);
        }
    }
}


/* 
FUNCTIONALLITY: 
    Obtain the LU decomposition of a matrix by gaussian elimination.

INPUT VARIABLES:
    a       matrix (leftside upper, rightside lower)
    n       matrices dimensions
    ipvt    matrix with row permutations
    info    stores the index of the diagonal if the obtained value is 0 (non-invertible matrix)
    roll    0 to use unrolled functions
            1 to use roll functions
*/

static void dgefa(auto a, int n, auto ipvt, auto info, int roll){
    REAL t, s;
    int l, t2;

    // gaussian elimination with partial pivoting

    if (roll)
    {
        *info = 0;
        for (int k = 0; k < n; ++k){
            // find l = maximum value pivot index

            l = idamax(n, k, a->column(k), 1);

            s = (*a)(l,k);

            if (s != ZERO){

                // Interchange rows if necessary
                if (l != k){
                    for (int j = k; j < (n + k); ++j){
                        t = (*a)(l,j);
                        (*a)(l,j) = (*a)(k,j);
                        (*a)(k,j) = t;

                        if (j < n){
                            t2 = (*ipvt)(l,j);
                            (*ipvt)(l,j) = (*ipvt)(k,j);
                            (*ipvt)(k,j) = t2;
                        }
                    }
                }

                // Compute multipliers
                if (s != ONE){
                    t = ONE/s;
                    (*a)(k,n+k) = s;

                    dscal_r(n, k, t, a->row(k), 1);
                }

                // Row elimination
                for (int j = (k+1); j < n; ++j){
                    t = - (*a)(j,k);
                    (*a)(j,n+k) = -t;
                    daxpy_r(n, k, t, a->row(k), a->row(j), 1);
                }

            } else {
                (*info) = k;
            }

        } 
    } else {

        *info = 0;
        for (int k = 0; k < n; ++k){
            // find l = maximum value pivot index

            l = idamax(n, k, a->column(k), 1);

            s = (*a)(l,k);

            if (s != ZERO){

                // Interchange rows if necessary
                if (l != k){
                    for (int j = k; j < (n + k); ++j){
                        t = (*a)(l,j);
                        (*a)(l,j) = (*a)(k,j);
                        (*a)(k,j) = t;

                        if (j < n){
                            t2 = (*ipvt)(l,j);
                            (*ipvt)(l,j) = (*ipvt)(k,j);
                            (*ipvt)(k,j) = t2;
                        }
                    }
                }

                // Compute multipliers
                if (s != 1){
                    t = ONE/s;
                    (*a)(k,n+k) = s;

                    dscal_ur(n, k, t, a->row(k), 1);
                }

                // Row elimination
                for (int j = (k+1); j < n; ++j){
                    t = - (*a)(j,k);
                    (*a)(j,n+k) = -t;
                    daxpy_ur(n, k, t, a->row(k), a->row(j), 1);
                }

            } else {
                (*info) = k;
            }

        } 
    }

}


/* 
FUNCTIONALLITY: 
    Solve the system a*x = b or trans(a)*x = b.

INPUT VARIABLES:
    a       matrix (leftside upper, rightside lower)
    n       matrices dimensions
    ipvt    matrix with row permutations
    b       vector
    job     0 to solve a*x = b 
            1 to solve trans(a)*x = b
    roll    0 to use unrolled functions
            1 to use roll functions
*/

static void dgesl(auto a, int n, auto ipvt, auto b, int job, int roll){
    int nm1;

    if (roll){
        if (job == 0){
            // job = 0, solve a * x = b

            // Reorder vector b
            (*b) = (*ipvt) * (*b);

            // first solve l*y = b
            (*b)(0) = (*b)(0)/(*a)(0,n);
            for (int k = 1; k < n; ++k){
                daxpy_r(n, k, -(*b)(k-1), a->column(n+k-1), b->column(0), 1);
                (*b)(k) = (*b)(k)/(*a)(k,n+k);
            }

            // now solve u*x = y
            nm1 = n - 1;
            // Not necessary in this case to divide by the diagonal, since it is 1
            for (int k = nm1; k > 0; --k){
                daxpy_r(k, 0, -(*b)(k), a->column(k), b->column(0), 1);
            }
        } else {

            // job = 0, solve a * x = b
            auto a_t = a->t();
            REAL result = 0;

            // first solve trans(u)*y = b
            for (int k = 1; k < n; ++k){
                result = ddot_r(k, 0, a_t.row(k), b->column(0), 1);
                (*b)(k) = (*b)(k) - result;
            }

            // now solve trans(l)*x = y
            nm1 = n - 1;
            // Not necessary in this case to divide by the diagonal, since it is 1
            (*b)(nm1) = (*b)(nm1)/a_t(n+nm1,nm1);
            for (int k = nm1; k > 0; --k){
                result = ddot_r(n, k, a_t.row(n+k-1), b->column(0), 1);
                (*b)(k-1) = ((*b)(k-1) - result)/a_t(n+k-1,k-1);
            }

            // Reorder vector x
            (*b) = ipvt->t() * (*b);

        }
    } else {

        if (job == 0){
            // job = 0, solve a * x = b

            // Reorder vector b
            (*b) = (*ipvt) * (*b);

            // first solve l*y = b
            (*b)(0) = (*b)(0)/(*a)(0,n);
            for (int k = 1; k < n; ++k){
                daxpy_ur(n, k, -(*b)(k-1), a->column(n+k-1), b->column(0), 1);
                (*b)(k) = (*b)(k)/(*a)(k,n+k);
            }

            // now solve u*x = y
            nm1 = n - 1;
            // Not necessary in this case to divide by the diagonal, since it is 1
            for (int k = nm1; k > 0; --k){
                daxpy_ur(k, 0, -(*b)(k), a->column(k), b->column(0), 1);
            }
        } else {

            // job = 0, solve a * x = b
            auto a_t = a->t();
            REAL result = 0;

            // first solve trans(u)*y = b
            for (int k = 1; k < n; ++k){
                result = ddot_ur(k, 0, a_t.row(k), b->column(0), 1);
                (*b)(k) = (*b)(k) - result;
            }

            // now solve trans(l)*x = y
            nm1 = n - 1;
            // Not necessary in this case to divide by the diagonal, since it is 1
            (*b)(nm1) = (*b)(nm1)/a_t(n+nm1,nm1);
            for (int k = nm1; k > 0; --k){
                result = ddot_ur(n, k, a_t.row(n+k-1), b->column(0), 1);
                (*b)(k-1) = ((*b)(k-1) - result)/a_t(n+k-1,k-1);
            }

            // Reorder vector x
            (*b) = ipvt->t() * (*b);

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