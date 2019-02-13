# Rational B-Spline Curves

A rational B-spline curve is the projection of a nonrational B-spline curve defined in four-dimensional homogeneous coordinate space back to the three-dimensional physical space.

A rational B-spline curve is defined as

```math
P(t) = \sum_{i=1}^{n+1} B_i^hN_{i,k}(t)
```

where ``B_i^h`` are the four-dimensional homogeneous control polygon vertices for the nonrational four-dimensional B-spline curve and ``N_{i,k}(t)`` is the nonrational B-spline basis function.

Projecting into the three-dimensional space by dividing through by thr homogeneous coordinate yelds

```math
P(t) = \dfrac{\sum_{i=1}^{n+1} B_ih_iN_{i,k}(t)}{\sum_{i=1}^{n+1}h_iN_{i,k}(t)} = \sum_{i=1}^{n+1}B_iR_{i,k}(t)
```

here the ``B_i`` are the three-dimensional control polygon vertices and ``R_{i,k}(t) = \dfrac{h_iN_{i,k}(t)}{\sum_{i=1}^{n+1}h_iN_{i,k}(t)}`` are the rational B-spline basis function where ``h_i>0`` for every value of ``i``.

## RB-Spline Properties

RB-spline are a generalization of nonrational B-spline, thus they carry forward nearly all the analytic and geometric characteristics:
 - Every rational basis function is positive or zero for all parameter values.
 - The sum of the rational B-spline basis functions for any value of ``t`` is one.
 - Each basis function has one maximum, except for the first order basis.
 - The maximum order of the rational B-spline is equal to the number of control polygon vertices.
 - A RB-spline generally follows the shape of the control polygon.

## RB-Spline Basis Functions

To generate RB-spline basis functions and curves are used open uniform, periodic uniform and nonuniform knot vectors.

The homogeneous coordinates ``h_i``, also called homogeneous weight factors, provide additional blendig capability, i.e. as a certain weight ``h_j`` increases the curve is pulled closer to the polygon vertex ``B_j``.

## RB-Spline Derivatives

The derivatives are obtained by formal differantiation of the curve's function

```math
P'(t) = \sum_{i=1}^{n+1}B_iR'_{i,k}(t)
```

where

```math
R'_{i,k}(t)=\dfrac{h_iN'_{i,k}(t)}{\sum_{i=1}^{n+1}h_iN_{i,k}}-\dfrac{h_iN_{i,k}\sum_{i=1}^{n+1}h_iN'_{i,k}}{(\sum_{i=1}^{n+1}h_iN_{i.k})^2}
```
For example, evaluating this result at ``t=0`` and ``t=n-k+2``

```math
P'(0)=(k-1)\dfrac{h_2}{h_1}(B_2-B_1)
```
and
```math
P'(n-k+2)=(k-1)\dfrac{h_n}{h_{n+1}}(B_{n+1}-B_n)
```

Higher order dervatives are obtained in a similar manner.
