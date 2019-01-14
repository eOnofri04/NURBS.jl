# Bezier Curves

Bezier Curves are a special kind of NURBS curves build upon a controll polygon ``B`` and a Bernstein basis ``J_{n,i}``.

Mathematically, a parametric Bézier curve is defined by

```math
P(t) = \sum_{i=0}^n B_iJ_{n,i}(t) \quad 0 \leq t \leq 1
```

Where the Bernstein Basis is defined as

```math
J_{n, i}(t) = \binom ni t^i(1-t)^{n-i}
```

where the convention ``(0)^0 = 1`` and ``0! = 1`` have been made.

In particular ``J_{n,i}`` is the ``i``-th base function of order ``n``, while ``n`` is also the number of segments of the the controll polygon (number of points minus one).

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

```math
P(t) = [F][G]
```

where

```math
[F] = [J_{n,0}, \dots, J_{n, n}]	\qquad	[G]^t = [B_0, \dots, B_n]
```

Moreover is possible to collect the coefficient of the Basis in a square matrix ``N \in \mathcal M_{n,n}(\mathbb R)`` obtaining:

```math
P(t) = [t^n, t^{n-1}, \dots, t, 1] [N] [G] = [T][N][G]
```

where

```math
[N] = \begin{bmatrix}
	\binom n0\binom nn (-1)^n & \binom n1 \binom{n-1}{n-1}(-1)^{n-1} & \dots & \binom nn \binom{n-n}{n-n}(-1)^0\\
	\binom n0\binom n{n-1} (-1)^{n-1} & \binom n1 \binom{n-1}{n-2}(-1)^{n-2} & \dots & 0\\
	\vdots & \vdots & \ddots & \vdots\\
	\binom n0\binom n1 (-1)^1 & \binom n1 \binom{n-1}0(-1)^0 & \dots & 0\\
	\binom n0\binom n0 (-1)^0 & 0 & \dots & 0\\
\end{bmatrix}
```

It is also possible to decompose the matrix ``[N]`` even further in the product of two matrices:

```math
[N] = [C][D] \quad \Rightarrow \quad P(t) = [T][C][D][G]
```

where

```math
[C] = \begin{bmatrix}
	\binom nn (-1)^n & \binom{n-1}{n-1}(-1)^{n-1} & \dots & \binom{n-n}{n-n}(-1)^0\\
	\binom n{n-1} (-1)^{n-1} & \binom{n-1}{n-2}(-1)^{n-2} & \dots & 0\\
	\vdots & \vdots & \ddots & \vdots\\
	\binom n1 (-1)^1 & \binom{n-1}0(-1)^0 & \dots & 0\\
	\binom n0 (-1)^0 & 0 & \dots & 0\\
\end{bmatrix}
```

```math
[D] = \begin{bmatrix}
	\binom n0 & 0 & \dots & 0\\
	0 & \binom n1 & \dots & 0\\
	\vdots & \vdots & \ddots & \vdots\\
	0 & 0 & \dots & \binom nn\\
\end{bmatrix}
```

## Bezier Derivatives

The two derivatives could be obtained starting from the original function:

```math
P'(t) = \sum_{i=0}^n B_iJ'_{n,i}(t)
```
```math
P''(t) = \sum_{i=0}^n B_iJ''_{n,i}(t)
```

where the two derivatives could be obtained with the following formulas:

```math
J'_{n,i}(t) = \frac{1-nt}{t(1-t)}J_{n,i}(t)
```
```math
J''_{n,i}(t) = \frac{(i-nt)^2-nt^2-i(1-2t)}{t^2(1-t)^2} J_{n,i}(t)
```
defined in all the points but the first and the last (where `t = 0` and `t = 1`).

### Derivatives Matrix Representation

In order to enhance the matrix representation it is usefull to represent the derivatives in two matrices:
```math
[Der1] = \begin{bmatrix}
	\frac{-nt_{(1)}}{t_{(1)}*(1-t_{(1)})}	&	\dots	&	\dots	&	\dots	&	\frac{n-nt_{(0)}}{t_{(0)}(1-t_{(0)})}	\\
%
	\vdots	&	\ddots	&					\vdots						&	\vdots	&	\vdots	\\
	\dots	&	\dots	&	\frac{(j-1)-nt_{(i)}}{t_{(i)}(1-t_{(i)})}	&	\dots	&	\dots	\\
	\vdots	&	\vdots	&					\vdots						&	\ddots	&	\vdots	\\
%
	\frac{-nt_{(bpts-1)}}{t_{(bpts-1)}(1-t_{(bpts-1)})}	&	\dots	&	\dots	&	\dots	&	\frac{n-nt_{(bpts-1)}}{t_{(bpts-1)}(1-t_{(bpts-1)})} \\
\end{bmatrix}
```
and
```math
[Der2] = begin{bmatrix}
	\frac{(-nt_{(1)})^2- nt_{(1)}^2}{t_{(1)}^2*(1-t_{(1)})^2}	&	\dots	&	\dots	&	\dots	&	\frac{(n-nt_{(0)})^2-nt_{(0)}^2-n(1-2t_{(0)})}{t_{(0)}^2(1-t_{(0)})^2}	\\
%
	\vdots	&	\ddots	&									\vdots											&	\vdots	&	\vdots	\\
	\dots	&	\dots	&	\frac{((j-1)-nt_{(i)})^2-nt_{(i)}^2-(j-1)(1-2t_{(i)})}{t_{(i)}^2(1-t_{(i)})^2}	&	\dots	&	\dots	\\
	\vdots	&	\vdots	&									\vdots											&	\ddots	&	\vdots	\\
%
	\frac{(-nt_{(bpts-1)})^2-nt_{(bpts-1)}^2}{t_{(bpts-1)}^2(1-t_{(bpts-1)})^2}	&	\dots	&	\dots	&	\dots	&	\frac{(n-nt_{(bpts-1)})^2-nt_{(bpts-1)}^2-n(1-2t_{(bpts-1)})}{t_{(bpts-1)}^2(1-t_{(bpts-1)})^2} \\
\end{bmatrix}
```

which could be scalar multiplied with the base matrix in order to obtain two derivative base matrices.
Of course the `T` matrix used must be reduced by eliminating the first and the last points.

The result of the operation is then:

```math
\begin{split}
	P'(t) =& 	Der1\cdot([T]_1^{dpts-1}[C][D])[G] \\
	P''(t) = &	Der2\cdot([T]_1^{dpts-1}[C][D])[G]
```




