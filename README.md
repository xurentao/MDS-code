# Modified Direct Swapping

This repository contains Fortran 77 and MATLAB implementations in [1] for the swapping problem of 2-by-2 blocks in the real Schur form. For any questions, please contact Rentao Xu, rtxu24@m.fudan.edu.cn.

- **DS:** A direct swapping algorithm for the reordering problem in the real Schur form. See *On swapping diagonal blocks in real Schur form*, Z. Bai and J. Demmel., 1993. The reference implementation was obtained from https://www.netlib.org/lapack/.
- **MDS:** A modifed direct swapping algorithm replacing Householder QR with pivoted Givens QR. 
- **RDS:** A refined direct swapping algorithm . See *Swapping 2x2 blocks in the Schur and generalized Schur form*, D. Camps, N. Mastronardi, R. Vandebril, and P. Van Dooren., 2019. The reference implementation was obtained from https://github.com/campsd

## 1. Project Structure

```
.
├── README.md                           # This file
│
├── LICENSE
│
├── Fortran 77/                         # Fortran 77 implementations
│   ├── dlaexc_g.f                      # Our algorithm
│   ├── dlaexc_h.f                      # LAPACK routine
│   ├── test_g.f                        # Test program for dlaexc_g
│   ├── test_h.f                        # Test program for dlaexc_h
│   └── test matrices/                  # Test matrices
│       ├── test matrix 1.dat
│       ├── test matrix 2.dat
│       ├── test matrix 3.dat
│       └── test matrix 4.dat
│
└── MATLAB/                             # MATLAB implementations
    ├── mds/                            # Our algorithm
    │   ├── gecp.m                      # GECP
    │   ├── mds.m                       # Main function
    │   └── fig_mds.m                   # Reproduce Figures 1 and 2 
    │
    └── rds/                            # Reference algorithm
        ├── cosine_sv.m
        ├── refine_swap_sevp.m          
        ├── swap22_sevp.m               # Main function 
        └── fig_rds.m                   # Reproduce Figures 1 and 2
```

## 2. Requirements

The dependencies and versions listed below reflect the specific environments used by the author to develop and test this project.

**Fortran 77 Environment:**
- **Fortran Compiler:** `gfortran` (version 12.2.0 or later)
- **BLAS Library:** Reference BLAS
- **LAPACK Library:** Reference LAPACK (version 3.11.0 or later)

**MATLAB Environment:**
- **MATLAB:** Version R2024b or later

## 3. Algorithm Implementations

- **DS:** In `dlaexc_h.f`, changes are only made for the case `N1 = 2, N2 = 2`. An output argument
  ```
           TOL = MAX( ABS( D( 3, 1 ) ), ABS( D( 3, 2 ) ), 
       $              ABS( D( 4, 1 ) ), ABS( D( 4, 2 ) ) ) / DNORM
  ```
  is added to report the error. Besides, the orthogonal similarity transformation would still be applied to the test matrix even `INFO = 1`. 
- **MDS:** 
    - In `dlaexc_g.f`, changes are only made for the case `N1 = 2, N2 = 2`. The QR factorization is replaced with proposed pivoted Givens QR. The remaining changes to `dlaexc_g.f` are similar to thoes made for `dlaexc_h`. 
    - In `mds.m`, `gecp.m` is used to solve a linear system with complete pivoting. The standardization of both 2-by-2 diagonal blocks are skipped to keep consistent with `swap22_sevp.m`. Note that `mds.m` only works for swapping 2-by-2 diagonal blocks.
- **RDS:** In `swap22_sevp.m`, the stability test is changed to
  ```
  norm_max = max(max(abs(A)));
  A = Q'*A*Q;
  tol = max(max(abs(A(3:4,1:2))))/norm_max;
  if tol > 10*eps
  ...
  ```
  And a similar change is applied to `refine_swap_sevp.m`. Moreover, `tol` in `swap22_sevp.m` is added as an output argument to report the error. In our test the transformation `Q` is computed via an SVD. 

All modifications in DS and RDS do not affect performances of both algorithms. They are only for numerical experiments.

## 4. Test Matrices

Test matrices 1--3 are provided as test matrix N.dat and match the numbering in [1]. Test matrix 4 is also provided for reference. It was used in Bai and Demmel's paper to demonstrate the superiority of the direct swapping algorithm over the iterative algorithm.

## 5. Compiling and Running

We show how to compile and run our Fortran programs. The following commands are written for the Windows Command Prompt (cmd) environment. Take `test_g.f` as an example. To compile `test_g.f`, type
```
gfortran -std=legacy test_g.f dlaexc_g.f -llapack -lblas -o test_g
```

To run `test_g.exe` with `test matrix N.dat`, type
```
test_g "test matrices/test matrix N.dat"
```

## 6. Reproducing Results

The following shows the expected output when running `test_g` with `test matrix 1.dat`:
```
A before swapping:
  5.270536E-04  2.418000E-12  1.964442E+00  1.963905E+00
 -7.781646E+00  5.270536E-04  1.021202E+00  9.426382E-01
  0.000000E+00  0.000000E+00  1.784755E-08  3.640901E-08
  0.000000E+00  0.000000E+00 -3.640901E-08  1.784755E-08

eigenvalues:
Lambda1 =  5.270536000000000E-04 +  4.337743656326408E-06i

Lambda2 =  5.270536000000000E-04 + -4.337743656326408E-06i

Lambda3 =  1.784755000000000E-08 +  3.640901000000000E-08i

Lambda4 =  1.784755000000000E-08 + -3.640901000000000E-08i

A after swapping by DLAEXC_G:
  1.784755E-08  2.668324E+00 -7.310351E+00  1.387360E+00
 -4.967972E-16  1.784755E-08 -1.881260E-04 -9.337855E-01
  0.000000E+00  0.000000E+00  5.270536E-04  2.616004E+00
  0.000000E+00  0.000000E+00 -7.192637E-12  5.270536E-04

eigenvalues:
Lambda1 =  1.784755012858628E-08 +  3.640901012888724E-08i

Lambda2 =  1.784755012858628E-08 + -3.640901012888724E-08i

Lambda3 =  5.270535999998065E-04 +  4.337737849039124E-06i

Lambda4 =  5.270535999998065E-04 + -4.337737849039124E-06i

tol:
  5.706880136285596E-17

```

## 7. Acknowledgement

DeepSeek was used to generate the initial draft of `gecp.m`, improve the output format for `test_g.f` and `test_h.f`, and suggest the outline of this README. The generated codes were reviewed and tested to ensure correctness.

## 8. License

The MDS implementations (`dlaexc_g.f` and `mds.m`) are released under the [MIT License](LICENSE).

## 9. Reference
[1] Meiyue Shao and Rentao Xu. ``A Modified Direct Eigenvalue Swapping Algorithm with Guaranteed Stability.'' arXiv: [2610.30795](https://arxiv.org/abs/2610.30795).