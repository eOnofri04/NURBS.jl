# Bézier Surfaces

Bézier Surfaces are special kind of NURBS build upon a controll net ``B`` and a couple of Bernstein basis ``J_{n,i}(u)`` and ``K_{m,j}(w)``.

Mathematically, a parametric Bézier surface is defined by

```math
Q(u, w) = \sum_{i=0}^n \sum_{j=0}^m B_{i,j} J_{n,i}(u) K_{m,j}(w)
```
where ``u`` and ``w`` are two parametric directions and ``n = npts - 1`` while ``m = mpts - 1``.

We point out that the Bernstein Basis are defined as

```math
\begin{split}
	J_{n, i}(u) = \binom ni u^i(1-u)^{n-i}\\
	K_{m, j}(w) = \binom mj w^j(1-w)^{m-j}
\end{split}
```

## Bézier Properties

Many properties of Bézier Surfaces are known:
 - The degree of the surface in each parametric direction is one less than the number of control net vertices in that direction.
 - The continuity of the surface in each parametric direction is two less than the number of control net vertices in that direction.
 - The surface generally follows the shape of the control net.
 - Only the corner points of the control net and the resulting Bézier surface are coincident.
 - The surface is contained within the convex hull of the control net.
 - The surface does not exhibit the variation-diminishing property. The variation-diminishing property for bivariant surfaces is both undefined and unknown.
 - The surface is invariant under an affine transformation.

## Matrix Representation

Like in Bézier Curves, also Bézier surfaces are very well described by a matrix point of view. Here we could simplify the calculus by writing:
```math
	Q(u, w) = [U][N][B][M]^T[W]^T
```
where
```math
\begin{split}
	[U] =& \begin{bmatrix} 
		u^n & u^{n-1} & \dots & 1
	\end{matrix}\\
	[W] =& \begin{bmatrix} 
		w^n & u^{w-1} & \dots & 1
	\end{matrix}\\
	[B] =& \begin{bmatrix} 
		B_{0,0}	&	\dots	&	B_{0,m}	\\
		\vdots	&	\ddots	&	\vdots	\\
		B_{n,0}	&	\dots	&	B_{n,m}
	\end{matrix}
\end{split}
```

and ``N`` and ``M`` are the coefficient matrix of the base, so given by:
```math
	begin{split}
	[N] = \begin{bmatrix}
		\binom n0\binom nn (-1)^n & \binom n1 \binom{n-1}{n-1}(-1)^{n-1} & \dots & \binom nn \binom{n-n}{n-n}(-1)^0\\
		\binom n0\binom n{n-1} (-1)^{n-1} & \binom n1 \binom{n-1}{n-2}(-1)^{n-2} & \dots & 0\\
		\vdots & \vdots & \ddots & \vdots\\
		\binom n0\binom n1 (-1)^1 & \binom n1 \binom{n-1}0(-1)^0 & \dots & 0\\
		\binom n0\binom n0 (-1)^0 & 0 & \dots & 0\\
	\end{bmatrix} \\
	[M] = \begin{bmatrix}
		\binom m0\binom mm (-1)^m & \binom m1 \binom{m-1}{m-1}(-1)^{m-1} & \dots & \binom mm \binom{m-m}{m-m}(-1)^0\\
		\binom m0\binom m{m-1} (-1)^{m-1} & \binom m1 \binom{m-1}{m-2}(-1)^{m-2} & \dots & 0\\
		\vdots & \vdots & \ddots & \vdots\\
		\binom m0\binom m1 (-1)^1 & \binom m1 \binom{m-1}0(-1)^0 & \dots & 0\\
		\binom m0\binom m0 (-1)^0 & 0 & \dots & 0\\
	\end{bmatrix} \\
```
or, in other worlds:
```math
\begin{split}
	\left(N_{i+1,j+1}\right)_{i,j=0}^n =& \begin{cases}
		\binom nj \binom{n-j}{n-i-j}(-1)^{n-i-j}	&	\mbox{if } 0 \leq i+j leq n\\
		0											&	\mbox{otherwise}
	\end{cases}\\
	\left(M_{i+1,j+1}\right)_{i,j=0}^m =& \begin{cases}
		\binom mj \binom{m-j}{m-i-j}(-1)^{m-i-j}	&	\mbox{if } 0 \leq i+j leq m\\
		0											&	\mbox{otherwise}
	\end{cases}\\
\end{split}
```

of course it is possible to decompose ``[N]`` and ``[M]`` in two more matrices, like for curves, leading to the following scheme:
```math
Q(u, w) = [U]([C][D])[B]([E][F])^T[W]^T = [U][C][D][B][F]^T[E]^T[W]^T
```

Much more interesting is the situation where ``m = n`` as ``[M] = [N]``

## Bézier Surface Derivatives

Of course is simple to determine the derivatives starting from the matrix notation.

We recall th results about the first and the second derivatives but the process could be extended further with no particular changes:
```math
\begin{split}
	\frac {\partial}{\partial u} Q(u,w) = &	[U'][N][B][M]^T[W]^T\\
	\frac {\partial}{\partial w} Q(u,w) = &	[U][N][B][M]^T[W']^T\\
	\frac {\partial^2}{\partial u\partial w} Q(u,w) = &	[U'][N][B][M]^T[W']^T\\
	\frac {\partial^2}{\partial u^2} Q(u,w) = &	[U''][N][B][M]^T[W]^T\\
	\frac {\partial^2}{\partial w^2} Q(u,w) = &	[U][N][B][M]^T[W'']^T
\end{split}
```

where, as expected, we have:

```math
\begin{split}
	[U'] =& [nu^{n-1} \ (n-1)u^{n-2} \ \dots \ 2n \ 1 \ 0]\\
	[W'] =& [mw^{m-1} \ (m-1)w^{m-2} \ \dots \ 2m \ 1 \ 0]\\
	[U''] =& [n(n-1)u^{n-2} \ (n-1)(n-2)u^{n-3} \ \dots \ 1 \ 0 \ 0]\\
	[W''] =& [m(m-1)w^{m-2} \ (m-1)(m-2)w^{m-3} \ \dots \ 1 \ 0 \ 0]\\
\end{split}
```