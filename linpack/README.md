# FOLDER DECRIPTION

In this folder, you can find different versions of the LINPACK benchmark. The programming languages used are C and C++.


## linpackc.new

Contains the code in C language from the original LINPACK benchmark. 

This code contains the following functions:
  - *linpack* : call the functions matgen, dgefa and dgesl to measure time.
  - *matgen* : initalizes matrix a and vector b.
  - *dgefa* : factors matrix a by gaussian elimination.
  - *dgesl* : solves the double precision system a * x = b  or  trans(a) * x = b using the factors computed by dgefa.
  - *daxpy_r* : returns the result of computing constant times a vector plus a vector.
  - *ddot_r* : forms the dot product of two vectors.
  - *dscal_r* : scales a vector by a constant.
  - *daxpy_ur* : unrolled version of daxpy_r. The computation is performed in groups of four vector components each time.
  - *ddot_ur* : unrolled version of ddot_r. - idamax: returns the index of the element with the maximum absolute value.The computation is performed in groups of five vector components each time.
  - *dscal_ur* : unrolled version of dscal_r. The computation is performed in groups of five vector components each time.
  - *idamax* : returns the index of the element with the maximum absolute value.
  - *second* : returns the current execution time in seconds.


## linpackcpp.new

Contains the code in C++ language based on _linpackc.new_. 

This code contains the following functions:
  - *linpack* : call the functions matgen, dgefa and dgesl to measure time.
  - *matgen* : initalizes matrix a and vector b.
  - *dgefa* - idamax: returns the index of the element with the maximum absolute value.: factors matrix a by gaussian elimination.
  - *dgesl* : solves the double precision system a * x = b  or  trans(a) * x = b using the factors computed by dgefa.
  - *daxpy_r* : returns the result of computing constant times a vector plus a vector.
  - *ddot_r* : forms the dot product of two vectors.
  - *dscal_r* : scales a vector by a constant.
  - *daxpy_ur* : unrolled version of daxpy_r. The computation is performed in groups of four vector components each time.
  - *ddot_ur* : unrolled version of ddot_r. The computation is performed in groups of five vector components each time.
  - *dscal_ur* : unrolled version of dscal_r. The computation is performed in groups of five vector components each time.
  - *idamax* : returns the index of the element with the maximum absolute value.


## linpackcpp_original

Contains the code in C++ language based on _linpackc.new_ using Linear Algebra library. 

This code contains the following functions:
  - *dgefa* : factors matrix a by gaussian elimination.
  - *dgesl* : solves the double precision system a * x = b  or  trans(a) * x = b using the factors computed by dgefa.
  - *daxpy_r* : returns the result of computing constant times a vector plus a vector.
  - *ddot_r* : forms the dot product of two vectors.
  - *dscal_r* : scales a vector by a constant.
  - *daxpy_ur* : unrolled version of daxpy_r. The computation is performed in groups of four vector components each time.
  - *ddot_ur* : unrolled version of ddot_r. The computation is performed in groups of five vector components each time.
  - *dscal_ur* : unrolled version of dscal_r. The computation is performed in groups of five vector components each time.
  - *idamax* : returns the index of the element with the maximum absolute value.


## linpackcpp_functions

Contains the code in C++ language with some changes with respect to the original script in _linpackc.new_. 

This code contains the following functions:
  - *linpack* : call the functions matgen, dgefa and dgesl to measure time.
  - *matgen* : initalizes matrix a and vector b.
  - *dgefa* : generates LU decomposition of matrix a.
  - *dgesl* : solves the double precision system a * x = b  or  trans(a) * x = b using the factors computed by dgefa.
  - *daxpy_r* : returns the result of computing constant times a vector plus a vector.
  - *ddot_r* : forms the dot product of two vectors.
  - *dscal_r* : scales a vector by a constant.
  - *daxpy_ur* : unrolled version of daxpy_r. The computation is performed in groups of four vector components each time.
  - *ddot_ur* : unrolled version of ddot_r. The computation is performed in groups of fi- idamax: returns the index of the element with the maximum absolute value.ve vector components each time.
  - *dscal_ur* : unrolled version of dscal_r. The computation is performed in groups of five vector components each time.
  - *idamax* : returns the index of the element with the maximum absolute value.


## linpackcpp

Contains the code in C++ language with some changes with respect to the original script in _linpackc.new_. 

This code contains the following functions:
  - *linpack* : call the functions matgen, dgefa and dgesl to measure time.
  - *matgen* : initalizes matrix a and vector b.
  - *dgefa* : generates LU decomposition of matrix a.
  - *dgesl* : solves the double precision system a * x = b  or  trans(a) * x = b using the factors computed by dgefa.
  - *idamax* : returns the index of the element with the maximum absolute value.- idamax: returns the index of the element with the maximum absolute value.

In this version, less functions are defined since they are substituted by predefined functions from Linear Algebra library.