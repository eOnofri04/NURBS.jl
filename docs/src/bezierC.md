# Bezier Curves

Bezier Curves are a special kind of NURBS curves build upon a controll polygon `B` and a Bernstein basis `J_{n,i}`.

Mathematically, a parametric Bézier curve is defined by
``P(t) = \sum_{i=0}^n B_iJ_{n,i}(t) \quad 0 \leq t \leq 1``
Where the Bernstein Basis is defined as
``J_{n, i}(t) = \binom ni t^i(1-t)^{n-i}``
where the convention `(0)^0 = 1` and `0! = 1` have been made.

In particular `J_{n,i}` is the `i`-th base function of order `n`, while `n` is also the number of segments of the the controll polygon (number of points minus one).

## Bezier Properties

Whe have a small variety of properties, descending directly from the definition:
 - Base function are often real.
 - The degree of the polynomial curve is one less than the number of point of the controll polygon.
 - The curve generally follows the shape of the control polygon.
 - The edges of the curves are the edges of the control polygon.
 - The curve is located inside the convex hull of the controll polygon.
 - The curve is invariant under Affine transformation.

## Matrix Representation

The equation of a Bézier Curve could also be implemented as a Matrix Multiplication (particulary usefull in GPU computations).
``P(t) = [F][G]``
where
``[F] = [J_{n,0}, \dots, J_{n, n}]	\qquad	[G]^t = [B_0, \dots, B_n]``

Moreover is possible to collect the coefficient of the Basis in a square matrix `N \in \mathcal M_{n,n}(\mathbb R)` obtaining:
``P(t) = [t^n, t^{n-1}, \dots, t, 1] [N] [G] = [T][N][G] ``
where
``[N] = \begin{bmatrix}
	\binom n0\binom nn (-1)^n & \binom n1 \binom{n-1}{n-1}(-1)^{n-1} & \dots & \binom nn \binom{n-n}{n-n}(-1)^0\\
	\binom n0\binom n{n-1} (-1)^{n-1} & \binom n1 \binom{n-1}{n-2}(-1)^{n-2} & \dots & 0\\
	\vdots & \vdots & \iddots & \vdots\\
	\binom n0\binom n1 (-1)^1 & \binom n1 \binom{n-1}0(-1)^0 & \dots & 0\\
	\binom n0\binom n0 (-1)^0 & 0 & \dots & 0\\
\end{bmatrix}``

It is also possible to decompose the matrix `[N]` even further in the product of two matrices:
``[N] = [C][D] \quad \Rightarrow \quad P(t) = [T][C][D][G]``
where
``[N] = \begin{bmatrix}
	\binom nn (-1)^n & \binom{n-1}{n-1}(-1)^{n-1} & \dots & \binom{n-n}{n-n}(-1)^0\\
	\binom n{n-1} (-1)^{n-1} & \binom{n-1}{n-2}(-1)^{n-2} & \dots & 0\\
	\vdots & \vdots & \iddots & \vdots\\
	\binom n1 (-1)^1 & \binom{n-1}0(-1)^0 & \dots & 0\\
	\binom n0 (-1)^0 & 0 & \dots & 0\\
\end{bmatrix}``
``[N] = \begin{bmatrix}
	\binom n0 & 0 & \dots & 0\\
	0 & \binom n1 & \dots & 0\\
	\vdots & \vdots & \iddots & \vdots\\
	0 & 0 & \dots & \binom nn\\
\end{bmatrix}``

## Bezier Derivatives

The two derivatives could be obtained starting from the original function:
``P'(t) = \sum_{i=0}^n B_iJ'_{n,i}(t)``
``P''(t) = \sum_{i=0}^n B_iJ''_{n,i}(t)``
where the two derivatives could be obtained with the following formulas:
``J'_{n,i}(t) = \frac{1-nt}{t(1-t)}J_{n,i}(t)``
``J''_{n,i}(t) = \frac{(i-nt)^2-nt^2-i(1-2t)}{t^2(1-t)^2} J_{n,i}(t)``