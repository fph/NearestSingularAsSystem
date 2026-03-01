# Computing the nearest singular matrix using a nonlinear system approach

This Matlab package computes the nearest structured singular matrix, for
both sparsity structures and arbitrary (dense) linear structures.

Consider installing also https://github.com/fph/RiemannOracle alongside, 
to compare with those algorithms.

Usage for sparse matrices:

```matlab
>> A = mmread('orani678.mtx');
>> [AplusDelta Delta] = nearest_singular_sparse_newtonlike(A);
Computed initial values from singular values of A
## k = 1 
minres converged at iteration 1610 to a solution with relative residual 0.01.
line search: accepted alpha=1; stopping=5.5231e-05, norm(rhs)=0.00537761
## k = 2 
minres converged at iteration 2160 to a solution with relative residual 0.01.
line search: accepted alpha=1; stopping=5.68615e-07, norm(rhs)=5.5231e-05
## k = 3 
minres converged at iteration 2197 to a solution with relative residual 0.01.
line search: accepted alpha=1; stopping=5.67628e-09, norm(rhs)=5.68615e-07
## k = 4 
minres converged at iteration 2218 to a solution with relative residual 0.0099.
line search: accepted alpha=1; stopping=5.63395e-11, norm(rhs)=5.67628e-09
## k = 5 
minres converged at iteration 2233 to a solution with relative residual 0.01.
line search: accepted alpha=1; stopping=5.62916e-13, norm(rhs)=5.63397e-11
norm(rhs)=5.62916e-13

normDelta =

    0.0268
```

There is also a version for arbitrary dense structures. We show an example 
which uses `autobasis` from https://github.com/fph/RiemannOracle to compute
an orthonormal basis for the perturbation space.

```matlab
>> A = toeplitz([1,2,3], [1 5 6])

A =

     1     5     6
     2     1     5
     3     2     1

>> addpath(genpath('../NearestUnstableMatlab/'));
>> P = autobasis(A)

P(:,:,1) =

    0.5774         0         0
         0    0.5774         0
         0         0    0.5774


P(:,:,2) =

         0         0         0
    0.7071         0         0
         0    0.7071         0


P(:,:,3) =

     0     0     0
     0     0     0
     1     0     0


P(:,:,4) =

         0    0.7071         0
         0         0    0.7071
         0         0         0


P(:,:,5) =

     0     0     1
     0     0     0
     0     0     0

>> [AplusDelta Delta] = nearest_singular_structured_dense_newtonlike(A, P)
Computed initial values from singular values of A
*** k = 1 
System solved
line search: accepted alpha=1; stopping=0.0942635, norm(rhs)=2.72217
*** k = 2 
System solved
line search: accepted alpha=1; stopping=0.000128979, norm(rhs)=0.0942635
*** k = 3 
System solved
line search: accepted alpha=1; stopping=1.17664e-10, norm(rhs)=0.000128979
*** k = 4 
System solved
line search: accepted alpha=1; stopping=2.08296e-15, norm(rhs)=1.17664e-10
*** k = 5 
System solved
line search: accepted alpha=0.5; stopping=1.27555e-15, norm(rhs)=2.08296e-15
*** k = 6 
System solved
line search: accepted alpha=0.5; stopping=9.99201e-16, norm(rhs)=1.27555e-15
*** k = 7 
System solved
Line search failed, cannot improve the solution anymore
*** k = 8 
System solved
Line search failed, cannot improve the solution anymore
*** k = 9 
System solved
Line search failed, cannot improve the solution anymore
*** k = 10 
System solved
Line search failed, cannot improve the solution anymore
Maximum number of iterations reached

normDelta =

    2.8165


AplusDelta =

    2.0636    3.6470    6.6693
    1.5393    2.0636    3.6470
    3.0776    1.5393    2.0636


Delta =

    1.0636   -1.3530    0.6693
   -0.4607    1.0636   -1.3530
    0.0776   -0.4607    1.0636
```
