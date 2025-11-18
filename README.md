✅ README Template: Program 4 — Linear Least Squares via Householder Reflections
📌 Overview

This project implements the Householder transformation method to solve linear least squares problems:

Exact solutions when the system is consistent

Minimum-norm projections when overdetermined

Testing includes several structured cases (square systems, inconsistent systems, polynomial regressions, Chebyshev bases).

🧠 Key Concepts

Householder reflectors

Orthogonal transformations

Numerical stability vs Gaussian elimination

Back substitution on upper-triangular 
𝑅
R

Residual analysis

🧮 Algorithm Summary

For each column 
𝑗
j:

Extract the active subvector 
𝑎
a

Construct the Householder vector

Apply reflector to:

Remaining submatrix of 
𝐴
A

Remaining subvector of 
𝑏
b

After 
𝑘
k reflections:

Extract 
𝑅
R

Solve 
𝑅
𝑥
=
𝑐
Rx=c

Residual is computed from 
∣
∣
𝑑
∣
∣
2
∣∣d∣∣
2
	​

.

🧪 Test Cases
✔️ 1. Square Non-Singular Systems

Solution error ≈ machine precision

Growth factor ≈ 1–2

✔️ 2. Overdetermined, consistent

Residual ≈ 0

Conditioning grows mildly with size

✔️ 3. Overdetermined, inconsistent

Solution error small

Residual reflects model mismatch (~1e-2)

✔️ 4. Growth Factor Stress Test

Householder vastly outperforms LU (Program 3)

✔️ 5. Arithmetic Mean

Perfect conditioning

Residual reflects data variance (expected)

✔️ 6. Polynomial Regression: Monomial Basis

Higher degree ⇒ lower residual

Higher degree ⇒ worse conditioning

✔️ 7. Chebyshev Basis

Chebyshev nodes dramatically improve conditioning

Better accuracy vs uniform nodes

🎯 Key Takeaways

Householder is numerically robust and stable

Performs well even on ill-conditioned problems

Superior to LU for least squares tasks
