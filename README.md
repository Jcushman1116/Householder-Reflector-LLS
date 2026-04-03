# Least Squares via Householder Reflectors

## Overview

Implements a linear least squares solver using Householder transformations. The
algorithm handles square non-singular systems, overdetermined systems where a
solution exists, and overdetermined systems where no exact solution exists and a
minimum-norm approximation is computed. Tested across multiple problem types
including polynomial regression with monomial and Chebyshev bases.

## Methodology

The algorithm reduces the matrix $A \in \mathbb{R}^{n \times k}$ to upper
triangular form $R$ by applying a sequence of $k$ Householder reflectors. For
each column $j$, the active sub-column is extracted and used to construct a
unit vector $v$ that reflects the column onto a multiple of $e_1$:

$$
v = a + \text{sign}(a_1)\|a\|_2 e_1, \qquad v \leftarrow \frac{v}{\|v\|_2}
$$

The reflector is then applied in-place to the trailing submatrix and the
right-hand side vector:

$$
A(j{:}n,\, j{:}k) \leftarrow A(j{:}n,\, j{:}k) - 2v(v^T A(j{:}n,\, j{:}k))
\qquad
b(j{:}n) \leftarrow b(j{:}n) - 2v(v^T b(j{:}n))
$$

After $k$ reflectors, $R = A(1{:}k,\, 1{:}k)$ is upper triangular and the
transformed right-hand side splits into $c = b(1{:}k)$ and residual
$d = b(k+1{:}n)$. The solution is obtained by back substitution on $Rx = c$:

$$
x_i = \frac{1}{R_{ii}}\left(c_i - \sum_{j=i+1}^{k} R_{ij} x_j\right),
\qquad i = k, k-1, \ldots, 1
$$

If $d = 0$ the solution is exact. Otherwise the algorithm returns the least
squares minimizer $x = \text{argmin}\|Ax - b\|_2$. A rank check on the diagonal
of $R$ guards against ill-conditioned inputs. Algorithm complexity is
$O(2nk^2 - \frac{2}{3}k^3)$ for overdetermined systems and $O(\frac{4}{3}n^3)$
when $n = k$.

Accuracy is evaluated using the 2-norm across four metrics:

$$
E_{\text{sol}} = \frac{\|x_{\text{true}} - x_{\text{comp}}\|_2}{\|x_{\text{true}}\|_2},
\quad
\gamma = \frac{\max|R(:)|}{\max|A(:)|},
\quad
E_{\text{res}} = \frac{\|b - Ax_{\text{comp}}\|_2}{\|b\|_2},
\quad
K_2(A) = \|A\|_2\|A^{-1}\|_2
$$

## Test Cases

**Case 1 — Square non-singular:** Solution error at machine precision
($\sim 10^{-15}$), growth factor between 1 and 2, condition number grows
slowly with $n$.

**Case 2 — Overdetermined, $b \in \mathfrak{R}(A)$:** Near machine precision
solution and residual error, confirming exact consistency.

**Case 3 — Overdetermined, $b \notin \mathfrak{R}(A)$:** Solution error remains
near machine precision while residual error rises to $\sim 10^{-2}$, confirming
the algorithm correctly computes the projection. Stability is independent of
system consistency.

**Growth factor matrix:** Householder yields linear growth factor vs. exponential
blowup observed under no pivoting and partial pivoting LU from the prior
assignment.

**Arithmetic mean:** Near-zero solution error with perfect conditioning. Large
residuals reflect data variability around the mean, not algorithmic error.

**Monomial basis ($d = 1, 2, 3$):** Solution error near machine precision.
Residual decreases with degree. Conditioning deteriorates as $d$ increases due
to growing correlation between Gram matrix columns.

**Chebyshev basis ($d = 1, 2, 3$):** Lower conditioning than monomial basis
across all degrees. Chebyshev root points produce a near-diagonal Gram matrix
and outperform uniformly spaced points on both solution error and conditioning.

## Language

MATLAB

## How to Run

1. Ensure all three `.m` files are in the same directory: `LLSHouse`,
   `Program4TestDriver`, `Chebyshev`
2. Run `Program4TestDriver` — executes all test cases at once
3. Results are printed as tables showing condition number, mean/max solution
   error, mean/max growth factor, and mean/max residual error per dimension
4. Gram matrices for monomial and Chebyshev bases are printed at each tested
   degree
