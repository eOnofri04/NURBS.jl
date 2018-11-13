# Algorithms NURBS

Here we present the Algorithm list of the C code in [Rog00].

---

Modifier [-] stands for an algorithm defined in a previoous chapter;

Modifier [+] stands for an optional algorithm;

Modifier [?] stands for a missing pseudocode algorithm.

---

### Chapter 2: Bezier Curves

 - [ ] `bezier___`  Calculates a Bezier curve.
 - [ ] `dbezier__` Calculates a Bezier curve and it's derivatives.


### Chapter 3: B-spline Curves
 
 - [ ] `basis____` Calculates the B-spline basis functions.
 - [ ] `bsplfit__` [+] Subroutine to fit a B-spline curve using an open uniform knot vector. Call `basis, knot, param`.
 - [ ] `bspline__` Calculates a B-spline curve. Call `basis, knot`.
 - [ ] `bsplineu_` Calculates a periodic B-spline curve. Call `basis, knotu`.
 - [ ] `dbasis___` Calculates the B-spline basis functions and derivatives.
 - [ ] `dbasisu__` Calculates the periodic B-spline basis functions.
 - [ ] `dbspline_` Calculates a B-spline curve and derivatives. Call `dbasis, knot`.
 - [ ] `dbsplineu` Calculates a periodic B-spline curve and derivatives. Call `knotu, dbasisu`.
 - [ ] `knot_____` Calculates an open knot vector.
 - [ ] `knotc____` Calculates a chord length approximation open knot vector.
 - [ ] `knotu____` Calculates a periodic knot vector.
 - [ ] `matpbspl_` [+] Generate a B-spline curve using matrix methods and a periodic uniform knot vector. Call `nmatrix`.
 - [ ] `nmatrix__` [+] Calculate the general B-spline periodic basis matrix.
 - [ ] `param____` Calculates the chord length paramter values.


### Chapter 4: Rational B-spline (NURBS) Curves

 - [ ] `knot_____` [-] Calculates an open knot vector.
 - [ ] `knotu____` [-] Calculates a periodic knot vector.
 - [ ] `rbasis___` Calculates the rational B-spline basis functions.
 - [ ] `rbspline_` Calculates an open rational B-spline curve. Call `rbasis, knot`.
 - [ ] `rbsplinu_` Calculates a periodic rational B-spline curve. Call `rbasis, knotu`.


### Chapter 5: Bezier Surface

 - [ ] `bezsurf__` Calculates a Bezier surface.
 - [ ] `mbezsurf_` [+] Subroutine to calculate a Bezier surface using matrix methods.

### Chapter 6: B-spline Surfaces

 - [ ] `basis____` [-] Calculates the B-spline basis functions.
 - [ ] `bsplsurf_` Calculates a B-spline surface. Call `basis, knot`.
 - [ ] `bspsurfu_` Calculates a periodic B-spline surface. Call `knotu, basis`.
 - [ ] `dbsurf___` Calculates the B-spline basis functions and derivatives. Call `dbasis, knot`.
 - [ ] `dbasis___` [-] Calculates a B-spline surface and derivatives.
 - [ ] `knot_____` [-] Calculates an open knot vector.
 - [ ] `knotu____` [-] Calculates a periodic knot vector.


### Chapter 7: Rational B-spline (NURBS) Surfaces

 - [ ] `basis____` [-] Calculates the B-spline basis functions.
 - [ ] `_________` [?] Calculates and test the fast B-spline surface algorithm.
 - [ ] `knot_____` [-] Calculates an open knot vector.
 - [ ] `rbspsurf_` Calculates a rational B-spline (NURBS) surface. Call `knot, basis, sumrbas`.
 - [ ] `sumrbas__` [+] Calculate the sum of the nonrational basis functions.

## Summary

 **Chap** | **No.** | **Done**
----------|---------|----------
 2        |     2   |    0
 3        |  11+3   |    0
 4        |   5-2   |    0
 5        |   1+1   |    0
 6        |   7-4   |    0
 7        | 4-2+1   |    0
 **Tot**  |    27   |    0
