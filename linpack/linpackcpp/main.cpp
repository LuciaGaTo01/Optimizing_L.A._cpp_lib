#include <iostream>
#include <chrono>
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

#define m 10


static REAL linpack (long nreps, int job, int n, char *matOpt);
static void matgen(auto a, auto l, auto ipvt, auto b, int n);
static void dgefa(auto u, auto l, int n, auto ipvt, int *info);
static void dgesl(auto u, auto l, int n, auto ipvt, auto b, int job);
static int idamax (int n, int init, auto dx, int incx);

using namespace std::chrono_literals;

int main (void){
    long nreps;
    char buf[5], opt[5];
    int n;

    while(1){
        // MENU
        std::cout << "Introduce one of the following options: \n\n"         \
                  << "a - Dinamic matrix\n"                                 \
                  << "b - Fixed size matrix (" << m << "x" << m << ")\n"    \
                  << "c - General matrix (" << m << "x" << m << ")\n"       \
                  << "q - QUIT\n\n";
        std::cin >> opt;

        if (opt[0] == 'q' || opt[0] == 'Q') break;

        if (opt[0] == 'a' || opt[0] == 'A'){
            std::cout << "\n\nIntroduce matrix size: ";
            std:: cin >> buf;

            if (buf[0]=='\0' || buf[0]=='\n') n = m;
            else n = atoi(buf);
        } else {
            n = m;
        }

        std::cout << "\n\nLINPACK benchmark, " << PREC <<" precision. \n";
        std::cout << "Machine precision: " << BASE10DIG << " digits. \n";
        std::cout << "Array size " << n << " X " << n << ". \n";
        std::cout << "Average rolled and unrolled performance:\n\n";
        std::cout << "    Reps Time(s) DGEFA   DGESL  OVERHEAD    KFLOPS\n";
        std::cout << "----------------------------------------------------\n";

        nreps = 1;
        while (linpack(nreps,0, n, opt)<10.){
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
    n         matrices dimensions
    matOpt    matrix type to use
*/

static REAL linpack (long nreps, int job, int n, char *matOpt){

    REAL kflops, toverhead, ops, tdgefa = 0.0, tdgesl = 0.0, totalt = 0.0;
    int info = 0;

    STD_LA::general_matrix<REAL, m, m> l, u;
    STD_LA::general_matrix<int, m, m> ipvt;
    STD_LA::general_column_vector<REAL, m> b;

    if (matOpt[0] == 'a' || matOpt[0] == 'A'){
        STD_LA::dynamic_matrix<REAL> l, u;
        STD_LA::dynamic_matrix<int> ipvt;
        STD_LA::dynamic_column_vector<REAL> b;

        l.resize(n,n);
        u.resize(n,n);
        ipvt.resize(n,n);
        b.resize_rows(n);
    }

    if (matOpt[0] == 'b' || matOpt[0] == 'B'){
        STD_LA::fixed_size_matrix<REAL, m, m> l, u;
        STD_LA::fixed_size_matrix<int, m, m> ipvt;
        STD_LA::fixed_size_column_vector<REAL, m> b;
    }


    ops = ((2.0*n*n*n)/3.0+2.0*n*n);

    using clk = std::chrono::high_resolution_clock;
    using secs = std::chrono::seconds;
    auto t0 = clk::now();
    for (int i = 0; i < nreps; ++i){
        matgen(&u, &l, &ipvt, &b, n);
        auto t1 = clk::now();
        dgefa(&u, &l, n, &ipvt, &info);
        auto t2 = clk::now();
        tdgefa += (t2-t1)/1.s;
        auto t3 = clk::now();
        dgesl(&u, &l, n, &ipvt, &b, job);
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
    a       matrix
    l       lower matrix
    ipvt    matrix with row permutations
    b       vector
    n       matrices dimensions
*/

static void matgen(auto a, auto l, auto ipvt, auto b, int n){
    int init, i, j;

    init = 1325;
    for (i = 0; i < n; ++i){
        for (j = 0; j < n; ++j){
            init = (int) ((long) 3125 * (long)init % 65536L);
            (*a)(i,j) = (init - 32768.0)/16384.0;
            (*l)(i,j) = ZERO;
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
    u       upper triangular matrix
    l       lower triangular matrix
    n       matrices dimensions
    ipvt    matrix with row permutations
    info    stores the index of the diagonal if the obtained value is 0 (non-invertible matrix)
*/

static void dgefa(auto u, auto l, int n, auto ipvt, int *info){
    REAL t, s;
    int i;

    // gaussian elimination with partial pivoting

    *info = 0;
    for (int k = 0; k < n; ++k){

        // find l = maximum value pivot index
        i = idamax(n, k, u->column(k), 1);

        s = (*u)(i,k);

        if (s != ZERO){

            // Interchange rows if necessary
            if (i != k){
                u->swap_rows(i,k);
                l->swap_rows(i,k);
                ipvt->swap_rows(i,k);
            }

            // Compute multipliers
            if (s != ONE){
                t = ONE/s;
                (*l)(k,k) = s;

                u->row(k) = t * u->row(k);
            }

            // Row elimination
            for (int j = (k+1); j < n; ++j){
                t = - (*u)(j,k);
                (*l)(j,k) = -t;
                u->row(j) = u->row(j) + t*u->row(k);
            }

        } else {
            (*info) = k;
        }

    } 

}


/* 
FUNCTIONALLITY: 
    Solve the system a*x = b or trans(a)*x = b.

INPUT VARIABLES:
    u       upper triangular matrix
    l       lower triangular matrix
    n       matrices dimensions
    ipvt    matrix with row permutations
    b       vector
    job     0 to solve a*x = b 
            1 to solve trans(a)*x = b
*/

static void dgesl(auto u, auto l, int n, auto ipvt, auto b, int job){
    REAL t;
    // int i;
    int nm1;

    if (job == 0){

        // job = 0, solve a * x = b
        // Reorder vector b
        (*b) = (*ipvt) * (*b);

        // first solve l*y = b
        (*b)(0) = (*b)(0)/(*l)(0,0);
        for (int k = 1; k < n; ++k){
            t = (*l)(k,k);
            (*l)(k,k) = ZERO;
            auto mult = l->row(k) * b->column(0);
            (*b)(k) = ((*b)(k) - mult(0,0))/t;
            (*l)(k,k) = t;
        }

        // now solve u*x = y
        // Not necessary in this case to divide by the diagonal, since it is 1
        for (int k = (n-2); k >= 0; --k){
            t = (*u)(k,k);
            (*u)(k,k) = ZERO;
            auto mult = u->row(k) * b->column(0);
            (*b)(k) = (*b)(k) - mult(0,0);
            (*u)(k,k) = t;
        }
    } else {

        // job = 0, solve trans(a) * x = b

        auto u_t = u->t();
        auto l_t = l->t();

        // first solve trans(u)*y = b
        // Not necessary in this case to divide by the diagonal, since it is 1
        for (int k = 1; k < n; ++k){
            t = u_t(k,k);
            u_t(k,k) = ZERO;
            auto mult = u_t.row(k) * b->column(0);
            (*b)(k) = (*b)(k) - mult(0,0);
            u_t(k,k) = t;
        }

        // now solve trans(l)*x = y
        nm1 = n - 1;
        (*b)(nm1) = (*b)(nm1)/l_t(nm1,nm1);
        for (int k = (nm1-1); k >= 0; --k){
            t = l_t(k,k);
            l_t(k,k) = ZERO;
            auto mult = l_t.row(k) * b->column(0);
            (*b)(k) = ((*b)(k) - mult(0,0))/t;
            l_t(k,k) = t;
        }

        // Reorder vector x
        (*b) = ipvt->t() * (*b);

    }
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