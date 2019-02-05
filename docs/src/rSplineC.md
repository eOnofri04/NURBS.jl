# Rational B-Spline Curves

A rational B-spline curve is the projection of a nonrational B-spline curve defined in four-dimensional homogeneous coordinate space back to the three-dimensional physical space. A rational B-spline curve is defined as

```math
P(t) = \sum_{i=1}^{n+1} B_i^hN_{i,k}(t)
```

where ``B_i^h`` are the four-dimensional homogeneous control polygon vertices for the nonrational four-dimensional B-spline curve.