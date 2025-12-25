# Linear Least Squares via Householder Transformations

## Overview
This project implements the **Householder transformation algorithm** to solve linear least squares (LLS) problems. The algorithm transforms a matrix $A$ into an upper triangular matrix $R$ through a sequence of orthogonal reflectors, providing a numerically robust method for finding solution vectors in both consistent and inconsistent systems.

The implementation is evaluated on various matrix structures to analyze **numerical stability**, **accuracy**, and **growth behavior**. Performance is assessed using **relative solution error**, **residual error**, **growth factor**, and **condition number** metrics.

All routines are implemented in **MATLAB**.

---

## Methods Implemented

### 1. Householder Reflections
Transforms matrix $A \in \mathbb{R}^{n \times k}$ into an upper triangular matrix $R$.
- Uses a sequence of $k$ orthogonal Householder reflectors $H = I - 2vv^T$.
- Eliminates entries below the diagonal while preserving the 2-norm of the system.
- Maintains higher numerical stability compared to LU factorization methods.

### 2. Back Substitution Solver
Once $A$ is transformed to $R$ and $b$ is transformed to $c$, the system $Rx = c$ is solved.
- Efficiently computes the solution vector $x$.
- Includes a rank-check to ensure the matrix has full column rank before solving.

---

## Implementation Details

- **Efficiency**: Reflectors are applied directly to active sub-matrices and sub-vectors to optimize storage and computational cost.
- **Numerical Robustness**: Includes guards against ill-conditioned problems by checking if diagonal elements of $R$ are near machine precision.
- **Storage**: The algorithm is designed to avoid storing unnecessary data points.

### Computational Complexity
- **Overdetermined systems**: $O(2nk^2)$
- **Square systems ($n=k$)**: $O(\frac{4}{3}n^3)$
- **Flop count**: $2nk^2 - \frac{2}{3}k^3$

---

## Experimental Design

The algorithm was tested on matrices of increasing size using random values from the open interval (0,1). The following metrics were computed using the **standard 2-norm**:

### Metrics

**Relative Solution Error ($E_{sol}$)**
> $||x_{true} - x_{comp}||_2 / ||x_{true}||_2$

**Relative Residual Error ($E_{res}$)**
> $||b - Ax_{comp}||_2 / ||b||_2$

**Growth Factor ($\gamma$)**
> $max|R(:)| / max|A(:)|$

**Condition Number ($K_2(A)$)**
> $||A||_2 ||A^{-1}||_2$

---

## Test Cases and Results

The implementation was validated against three primary scenarios:

* **Square Non-Singular Matrices**: Produced machine precision accuracy with a max solution error of $O(10^{-15})$.
* **Overdetermined Systems ($b \in \mathcal{R}(A)$)**: Solution and residual errors remained near machine precision.
* **Overdetermined Systems ($b \notin \mathcal{R}(A)$)**: Successfully produced the least squares solution; residual error was $O(10^{-2})$, confirming the system's inconsistency while solution error remained low.
* **Polynomial Comparison**: Tested Monomial vs. Chebyshev bases. The Chebyshev matrix demonstrated higher numerical stability for higher-degree fits.

**Key Finding**: Unlike LU factorization, the Householder algorithm yields **linear growth** for growth-factor stress matrices rather than exponential growth, demonstrating superior stability.

---

## Files Included

* `LLSHouse.m` - The core function implementing the Householder LLS solver.
* `Chebyshev.m` - Script for comparing monomial and Chebyshev bases in regression.
* `Program4TestDriver.m` - Runs the full suite of test cases and generates performance metrics.

---

## Key Takeaways
* Householder transformations are highly stable for solving linear least squares.
* The algorithm's stability is independent of the consistency of the system.
* Orthogonal reflectors effectively control the growth factor compared to standard pivoting strategies.
