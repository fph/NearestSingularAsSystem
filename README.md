# Computing the nearest singular matrix using a nonlinear system approach

Work in progress.

Consider installing also https://github.com/fph/RiemannOracle in parallel, to compare with those algorithms.

Usage:

```matlab
>> A = mmread('orani678.mtx');
>> [AplusDelta Delta] = nearest_singular_sparse_newtonlike(A);
Computed initial value from singular values of A
*** k = 1 
line search: accepted alpha=1; stopping=0.000617654, norm(rhs)=0.00340858
*** k = 2 
line search: accepted alpha=1; stopping=2.03735e-06, norm(rhs)=0.000617654
*** k = 3 
line search: accepted alpha=1; stopping=1.20966e-11, norm(rhs)=2.03735e-06
*** k = 4 
line search: accepted alpha=1; stopping=9.04583e-16, norm(rhs)=1.20971e-11
norm(rhs)=1.70706e-14, frobnorm(Delta)=0.0268131
```
