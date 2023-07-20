
/*
**
** LINPACK.CPP      Linpack benchmark, calculates FLOPS.
**                  (FLoating Point Operations Per Second)
**
**
*/



#include <iostream>
#include <cfloat>
#include <cmath>
#include <string_view>
#include <limits>
#include <iomanip>
#include <sys/stat.h>
#include <vector>


template <typename T>
constexpr T ZERO = 0.0;
template <typename T>
constexpr T ONE = 1.0;
template <typename T>
constexpr int BASE10DIG = std::numeric_limits<T>::digits10;


/*
** Constant times a vector plus a vector.
** Jack Dongarra, linpack, 3/11/78.
** ROLLED version
*/
template <typename T>
void daxpy_r(int n,T da,T *dx,T *dy)

    {
    if (n <= 0)
        return;
    if (da == ZERO<T>)
        return;

    for (int i = 0;i < n; i++)
        dy[i] = dy[i] + da*dx[i];
    }

/*
** constant times a vector plus a vector.
** Jack Dongarra, linpack, 3/11/78.
** UNROLLED version
*/
template <typename T>
void daxpy_ur(int n,T da,T *dx,T *dy)

    {
    if (n <= 0)
        return;
    if (da == ZERO<T>)
        return;

    int m = n % 4;
    if ( m != 0)
        {
        for (int i = 0; i < m; i++)
            dy[i] = dy[i] + da*dx[i];
        if (n < 4)
            return;
        }
    for (int i = m; i < n; i = i + 4)
        {
        dy[i] = dy[i] + da*dx[i];
        dy[i+1] = dy[i+1] + da*dx[i+1];
        dy[i+2] = dy[i+2] + da*dx[i+2];
        dy[i+3] = dy[i+3] + da*dx[i+3];
        }
    }

template <typename T, bool is_rolled>
void daxpy (int n, T da, T *dx, T *dy){
    if constexpr (is_rolled){
        daxpy_r(n, da, dx, dy);
    } else {
        daxpy_ur(n, da, dx, dy);
    }
}

/*
** Forms the dot product of two vectors.
** Jack Dongarra, linpack, 3/11/78.
** ROLLED version
*/
template <typename T>
T ddot_r(int n,T *dx,T *dy)

    {
    T dtemp = ZERO<T>;

    if (n <= 0)
        return(ZERO<T>);

    for (int i=0;i < n; i++)
        dtemp = dtemp + dx[i]*dy[i];
    return(dtemp);
    }

/*
** Forms the dot product of two vectors.
** Jack Dongarra, linpack, 3/11/78.
** UNROLLED version
*/
template <typename T>
T ddot_ur(int n,T *dx,T *dy)

    {
    T dtemp = ZERO<T>;

    if (n <= 0)
        return(ZERO<T>);

    int m = n % 5;
    if (m != 0)
        {
        for (int i = 0; i < m; i++)
            dtemp = dtemp + dx[i]*dy[i];
        if (n < 5)
            return(dtemp);
        }
    for (int i = m; i < n; i = i + 5)
        {
        dtemp = dtemp + dx[i]*dy[i] +
        dx[i+1]*dy[i+1] + dx[i+2]*dy[i+2] +
        dx[i+3]*dy[i+3] + dx[i+4]*dy[i+4];
        }
    return(dtemp);
    }

template <typename T, bool is_rolled>
T ddot (int n, T *dx, T *dy){
    if constexpr (is_rolled){
        return (ddot_r(n, dx, dy));
    } else {
        return (ddot_ur(n, dx, dy));
    }
}


/*
** Scales a vector by a constant.
** Jack Dongarra, linpack, 3/11/78.
** ROLLED version
*/
template <typename T>
void dscal_r(int n,T da,T *dx)

    {
    if (n <= 0)
        return;

    for (int i = 0; i < n; i++)
        dx[i] = da*dx[i];
    }

/*
** Scales a vector by a constant.
** Jack Dongarra, linpack, 3/11/78.
** UNROLLED version
*/
template <typename T>
void dscal_ur(int n,T da,T *dx)

    {
    if (n <= 0)
        return;

    int m = n % 5;
    if (m != 0)
        {
        for (int i = 0; i < m; i++)
            dx[i] = da*dx[i];
        if (n < 5)
            return;
        }
    for (int i = m; i < n; i = i + 5)
        {
        dx[i] = da*dx[i];
        dx[i+1] = da*dx[i+1];
        dx[i+2] = da*dx[i+2];
        dx[i+3] = da*dx[i+3];
        dx[i+4] = da*dx[i+4];
        }
    }

template <typename T, bool is_rolled>
void dscal (int n,T da,T *dx){
    if constexpr (is_rolled){
        dscal_r(n, da, dx);
    } else {
        dscal_ur(n, da, dx);
    }
}


/*
** Finds the index of element having max. absolute value.
** Jack Dongarra, linpack, 3/11/78.
*/
template <typename T>
int idamax(int n,T *dx)

    {
    if (n < 1)
        return(-1);
    if (n ==1 )
        return(0);
    
    int itemp = 0;
    T dmax = std::fabs(dx[0]);
    for (int i = 1; i < n; i++)
        if(std::fabs(dx[i]) > dmax)
            {
            itemp = i;
            dmax = std::fabs(dx[i]);
            }
    return (itemp);
    }


/*
**
** DGEFA benchmark
**
** We would like to declare a[][lda], but c does not allow it.  In this
** function, references to a[i][j] are written a[lda*i+j].
**
**   dgefa factors a double precision matrix by gaussian elimination.
**
**   dgefa is usually called by dgeco, but it can be called
**   directly with a saving in time if  rcond  is not needed.
**   (time for dgeco) = (1 + 9/n)*(time for dgefa) .
**
**   on entry
**
**      a       REAL precision[n][lda]
**              the matrix to be factored.
**
**      lda     integer
**              the leading dimension of the array  a .
**
**      n       integer
**              the order of the matrix  a .
**
**   on return
**
**      a       an upper triangular matrix and the multipliers
**              which were used to obtain it.
**              the factorization can be written  a = l*u  where
**              l  is a product of permutation and unit lower
**              triangular matrices and  u  is upper triangular.
**
**      ipvt    integer[n]
**              an integer vector of pivot indices.
**
**      info    integer
**              = 0  normal value.
**              = k  if  u[k][k] .eq. 0.0 .  this is not an error
**                   condition for this subroutine, but it does
**                   indicate that dgesl or dgedi will divide by zero
**                   if called.  use  rcond  in dgeco for a reliable
**                   indication of singularity.
**
**   linpack. this version dated 08/14/78 .
**   cleve moler, university of New Mexico, argonne national lab.
**
**   functions
**
**   blas daxpy,dscal,idamax
**
*/
template <typename T, bool is_rolled>
void dgefa(T *a,int lda,int n,int *ipvt,int *info)

    {
    T t;
    int kp1,l;

    /* gaussian elimination with partial pivoting */

    *info = 0;
    int nm1 = n - 1;
    if (nm1 >=  0)
        for (int k = 0; k < nm1; k++)
            {
            kp1 = k + 1;

            /* find l = pivot index */

            l = idamax(n-k,&a[lda*k+k]) + k;
            ipvt[k] = l;

            /* zero pivot implies this column already
                triangularized */

            if (a[lda*k+l] != ZERO<T>)
                {

                /* interchange if necessary */

                if (l != k)
                    {
                    t = a[lda*k+l];
                    a[lda*k+l] = a[lda*k+k];
                    a[lda*k+k] = t;
                    }

                /* compute multipliers */

                t = -ONE<T>/a[lda*k+k];
                dscal<T,is_rolled>(n-(k+1),t,&a[lda*k+k+1]);

                /* row elimination with column indexing */

                for (int j = kp1; j < n; j++)
                    {
                    t = a[lda*j+l];
                    if (l != k)
                        {
                        a[lda*j+l] = a[lda*j+k];
                        a[lda*j+k] = t;
                        }
                    daxpy<T,is_rolled>(n-(k+1),t,&a[lda*k+k+1],&a[lda*j+k+1]);
                    }
                }
            else
                (*info) = k;
            }
    ipvt[n-1] = n-1;
    if (a[lda*(n-1)+(n-1)] == ZERO<T>)
        (*info) = n-1;
    
    }


/*
**
** DGESL benchmark
**
** We would like to declare a[][lda], but c does not allow it.  In this
** function, references to a[i][j] are written a[lda*i+j].
**
**   dgesl solves the double precision system
**   a * x = b  or  trans(a) * x = b
**   using the factors computed by dgeco or dgefa.
**
**   on entry
**
**      a       double precision[n][lda]
**              the output from dgeco or dgefa.
**
**      lda     integer
**              the leading dimension of the array  a .
**
**      n       integer
**              the order of the matrix  a .
**
**      ipvt    integer[n]
**              the pivot vector from dgeco or dgefa.
**
**      b       double precision[n]
**              the right hand side vector.
**
**      job     integer
**              = 0         to solve  a*x = b ,
**              = nonzero   to solve  trans(a)*x = b  where
**                          trans(a)  is the transpose.
**
**  on return
**
**      b       the solution vector  x .
**
**   error condition
**
**      a division by zero will occur if the input factor contains a
**      zero on the diagonal.  technically this indicates singularity
**      but it is often caused by improper arguments or improper
**      setting of lda .  it will not occur if the subroutines are
**      called correctly and if dgeco has set rcond .gt. 0.0
**      or dgefa has set info .eq. 0 .
**
**   to compute  inverse(a) * c  where  c  is a matrix
**   with  p  columns
**         dgeco(a,lda,n,ipvt,rcond,z)
**         if (!rcond is too small){
**              for (j=0,j<p,j++)
**                      dgesl(a,lda,n,ipvt,c[j][0],0);
**         }
**
**   linpack. this version dated 08/14/78 .
**   cleve moler, university of new mexico, argonne national lab.
**
**   functions
**
**   blas daxpy,ddot
*/
template <typename T, bool is_rolled>
void dgesl(T *a,int lda,int n,int *ipvt,T *b,int job)

    {
    T    t;
    int     l;

    int nm1 = n - 1;
    if (job == 0)
        {

        /* job = 0 , solve  a * x = b   */
        /* first solve  l*y = b         */

        if (nm1 >= 1)
            for (int k = 0; k < nm1; k++)
                {
                l = ipvt[k];
                t = b[l];
                if (l != k)
                    {
                    b[l] = b[k];
                    b[k] = t;
                    }
                daxpy<T,is_rolled>(n-(k+1),t,&a[lda*k+k+1],&b[k+1]);
                }

        /* now solve  u*x = y */

        for (int kb = 0; kb < n; kb++)
            {
            int k = n - (kb + 1);
            b[k] = b[k]/a[lda*k+k];
            t = -b[k];
            daxpy<T,is_rolled>(k,t,&a[lda*k+0],&b[0]);
            }
        }
    else
        {

        /* job = nonzero, solve  trans(a) * x = b  */
        /* first solve  trans(u)*y = b             */

        for (int k = 0; k < n; k++)
            {
            t = ddot<T,is_rolled>(k,&a[lda*k+0],&b[0]);
            b[k] = (b[k] - t)/a[lda*k+k];
            }

        /* now solve trans(l)*x = y     */

        if (nm1 >= 1)
            for (int kb = 1; kb < nm1; kb++)
                {
                int k = n - (kb+1);
                b[k] = b[k] + ddot<T,is_rolled>(n-(k+1),&a[lda*k+k+1],&b[k+1]);
                l = ipvt[k];
                if (l != k)
                    {
                    t = b[l];
                    b[l] = b[k];
                    b[k] = t;
                    }
                }
        }
    }

/*
** For matgen,
** We would like to declare a[][lda], but c does not allow it.  In this
** function, references to a[i][j] are written a[lda*i+j].
*/
template <typename T>
void matgen(T *a,int lda,int n,T *b,T *norma)

    {
    long init = 1325;
    *norma = 0.0;
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++)
            {
            init = 3125L*init % 65536L;
            a[lda*j+i] = (init - 32768.0)/16384.0;
            *norma = (a[lda*j+i] > *norma) ? a[lda*j+i] : *norma;
            }
    for (int i = 0; i < n; i++)
        b[i] = 0.0;
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++)
            b[i] = b[i] + a[lda*j+i];
    }


template <typename T>
static void writeFile (T *a, T *b, int nrows, int ncols){
    FILE *fptra, *fptrb;
    char fileNameA[100], fileNameb[100];
    struct stat st;

    if (stat("./outputFiles", &st) == -1) {
        mkdir("./outputFiles", 0700);
    }   
    
    sprintf(fileNameA, "./outputFiles/output_a-%d.txt", nrows);
    sprintf(fileNameb, "./outputFiles/output_b-%d.txt", nrows);
    if ((fptra = fopen(fileNameA,"w")) == NULL || (fptrb = fopen(fileNameb,"w")) == NULL){
       printf("Error! opening file");

       // Program exits if the file pointer returns NULL.
       exit(1);
    }

    for (int i = 0; i < nrows; ++i){
        for (int j = 0; j < ncols; ++j){
            fprintf(fptra, "%.6f\t", a[ncols*i+j]);
        }
        fprintf(fptra, "\n");
        fprintf(fptrb, "%.6f\n", b[i]);
    }

    fclose(fptra);
    fclose(fptrb);

    snprintf(fileNameA, sizeof(fileNameA), "./outputFiles/output_a-%d.data", nrows);
    snprintf(fileNameb, sizeof(fileNameb), "./outputFiles/output_b-%d.data", nrows);
    if ((fptra = fopen(fileNameA,"wb")) == NULL || (fptrb = fopen(fileNameb,"wb")) == NULL){
        printf("Error! opening binary file");

        // Program exits if the file pointer returns NULL.
        exit(1);
    }

    fwrite(a, sizeof(a[0]), nrows*ncols, fptra);
    fwrite(b, sizeof(b[0]), nrows, fptrb);

    fclose(fptra);
    fclose(fptrb);

}


template <typename T>
T second(void)

    {
    return ((T)clock()/(T)CLOCKS_PER_SEC);
    }

template <typename T>
T linpack (long nreps, long arsize, char wr[5], T *a, T *b, int *ipvt)

    {
    T   norma,t1;
    int   info;

    int lda = arsize;
    int n = arsize/2;
    long arsize2d = arsize*arsize;
    T ops=((2.0*n*n*n)/3.0+2.0*n*n);

    T tdgesl=0;
    T tdgefa=0;
    T totalt=second<T>();
    for (long i=0;i<nreps;i++)
        {
        matgen(a,lda,n,b,&norma);
        t1 = second<T>();
        dgefa<T,true>(a,lda,n,ipvt,&info);
        tdgefa += second<T>()-t1;
        t1 = second<T>();
        dgesl<T,true>(a,lda,n,ipvt,b,0);
        tdgesl += second<T>()-t1;
        }
    for (long i=0;i<nreps;i++)
        {
        matgen(a,lda,n,b,&norma);
        t1 = second<T>();
        dgefa<T,false>(a,lda,n,ipvt,&info);
        tdgefa += second<T>()-t1;
        t1 = second<T>();
        dgesl<T,false>(a,lda,n,ipvt,b,0);
        tdgesl += second<T>()-t1;
        }

    totalt=second<T>()-totalt;

    if (wr[0] == 'Y' || wr[0] == 'y'){
        writeFile(a, b, lda, lda);
    }
    if (totalt<0.5 || tdgefa+tdgesl<0.2) return(0.);
    T kflops=2.*nreps*ops/(1000.*(tdgefa+tdgesl));
    T toverhead=totalt-tdgefa-tdgesl;
    if (tdgefa<0.) tdgefa=0.;
    if (tdgesl<0.) tdgesl=0.;
    if (toverhead<0.) toverhead=0.;

    std::cout << std::fixed;
    std::cout << std::setprecision(2);
    std::cout << nreps << " " << totalt << " " << 100.*tdgefa/totalt << "% "         \
    << 100.*tdgesl/totalt << "% " << 100.*toverhead/totalt << "% " << kflops << '\n';

    return(totalt);
    }


template <typename T>
void callLinpack (long arsize, const char *PREC){

    char wr[5];

    std::cout << "Would you like to save the obtained output into a file? (Y/N): ";
    std::cin >> wr;

    long arsize2d = arsize*arsize;
    
    std::vector<T> a;
    std::vector<T> b;
    std::vector<int> ipvt;
    try
    {
        a.reserve(arsize2d);
        b.reserve(arsize);
        ipvt.reserve(arsize);

        std::cout << "\n\nLINPACK benchmark, " << PREC <<" precision. \n";
        std::cout << "Machine precision: " << BASE10DIG<T> << " digits. \n";
        std::cout << "Array size " << arsize << " X " << arsize << ". \n";
        std::cout << "Average rolled and unrolled performance:\n\n";
        std::cout << "    Reps Time(s) DGEFA   DGESL  OVERHEAD    KFLOPS\n";
        std::cout << "----------------------------------------------------\n";

        long nreps = 1;
        while (linpack<T>(nreps,arsize,wr,a.data(),b.data(),ipvt.data())<10.)
            nreps *= 2;

    }
    catch (const std::bad_alloc& ba) 
    {
        std::cerr <<"ERROR: ";
        std::cerr << ba.what() << "\n";
    }

}


int main (){
    char buf[80], precission[10];
    long arsize;
    const char *PREC;

    while(1){
        std::cout << "Introduce the precission required: \n s - single precission \n d - double precission \n q - Quit \n";
        std::cin >> precission;

        if (precission[0]=='q' || precission[0]=='Q') break;

        std::cout << "Enter array size [200]:  ";
        std::cin >> buf;

        if (buf[0]=='\0' || buf[0]=='\n') arsize=200;
        else arsize=atoi(buf);

        arsize/=2;
        arsize*=2;
        if (arsize<10)
            {
            std::cout << "Too small.\n";
            continue;
            }

        if (precission[0] == 's' || precission[0]=='S'){
            PREC = "single";
            callLinpack <float> (arsize, PREC);
        } else {
            PREC = "double";
            callLinpack <double> (arsize, PREC);
        }

        std::cout << '\n';
    }
    
}