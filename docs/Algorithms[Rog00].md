# Algorithms NURBS

Here we present the Algorithm list of the C code in [Rog00].

---

Modifier [-] stands for an algorithm defined in a previoous chapter;

Modifier [+] stands for an optional algorithm;

Modifier [?] stands for a missing pseudocode algorithm.

---

### Chapter 2: Bezier Curves

 - [ ] [E] `bezier___`  Calculates a Bezier curve.
 - [ ] [E] `dbezier__` Calculates a Bezier curve and it's derivatives.


### Chapter 3: B-spline Curves
 
 - [X] `basis____` Calculates the B-spline basis functions.
 - [ ] [GM] `bsplfit__` [+] Subroutine to fit a B-spline curve using an open uniform knot vector. Call `basis, knot, param`.
 - [X] `bspline__` Calculates a B-spline curve. Call `basis, knot`.
 - [X] `bsplineu_` Calculates a periodic B-spline curve. Call `basis, knotu`.
 - [X] `dbasis___` Calculates the B-spline basis functions and derivatives.
 - [X] `dbasisu__` Calculates the periodic B-spline basis functions.
 - [X] `dbspline_` Calculates a B-spline curve and derivatives. Call `dbasis, knot`.
 - [X] `dbsplineu` Calculates a periodic B-spline curve and derivatives. Call `knotu, dbasisu`.
 - [X] `knot_____` Calculates an open knot vector.
 - [ ] [G] `knotc____` Calculates a chord length approximation open knot vector.
 - [X] `knotu____` Calculates a periodic knot vector.
 - [ ] [P] `matpbspl_` [+] Generate a B-spline curve using matrix methods and a periodic uniform knot vector. Call `nmatrix`.
 - [ ] [P] `nmatrix__` [+] Calculate the general B-spline periodic basis matrix.
 - [ ] [P] `param____` Calculates the chord length paramter values.
 - [ ] [G] `raise23` [+]
 - [ ] [G] `raise12` [+]
 - [ ] [GM] `remknot` [+]
 


### Chapter 4: Rational B-spline (NURBS) Curves

 - [X] `knot_____` [-] Calculates an open knot vector.
 - [X] `knotu____` [-] Calculates a periodic knot vector.
 - [X] `rbasis___` Calculates the rational B-spline basis functions.
 - [X] `rbspline_` Calculates an open rational B-spline curve. Call `rbasis, knot`.
 - [X] `rbsplinu_` Calculates a periodic rational B-spline curve. Call `rbasis, knotu`.


### Chapter 5: Bezier Surface

 - [ ] `bezsurf__` Calculates a Bezier surface.
 - [ ] `mbezsurf_` [+] Subroutine to calculate a Bezier surface using matrix methods.

### Chapter 6: B-spline Surfaces

 - [X] `basis____` [-] Calculates the B-spline basis functions.
 - [ ] [P] `bsplsurf_` Calculates a B-spline surface. Call `basis, knot`.
 - [ ] [P] `bspsurfu_` Calculates a periodic B-spline surface. Call `knotu, basis`.
 - [ ] [P] `dbsurf___` Calculates the B-spline basis functions and derivatives. Call `dbasis, knot`.
 - [X] `dbasis___` [-] Calculates a B-spline surface and derivatives.
 - [X] `knot_____` [-] Calculates an open knot vector.
 - [X] `knotu____` [-] Calculates a periodic knot vector.


### Chapter 7: Rational B-spline (NURBS) Surfaces

 - [X] `basis____` [-] Calculates the B-spline basis functions.
 - [ ] `frbsurf__` [?] Calculates and test the fast B-spline surface algorithm. Call `basis, knot, sumrbas`
 - [X] `knot_____` [-] Calculates an open knot vector.
 - [ ] [P] `rbspsurf_` Calculates a rational B-spline (NURBS) surface. Call `knot, basis, sumrbas`.
 - [ ] [P] `sumrbas__` [+] Calculate the sum of the nonrational basis functions.
 - [ ] [GM] `bsurfnaive`
 - [ ] [GM] `rbsurf`

## Summary

 **Chap** | **No.** | **Done**
----------|---------|----------
 2        |     2   |    0
 3        |    17   |    0
 4        |     5   |    0
 5        |     2   |    0
 6        |     7   |    0
 7        |     6   |    0
 **Tot**  |    39   |    0
