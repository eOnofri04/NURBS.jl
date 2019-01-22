# Algorithms NURBS

Here we present the Algorithm list of the C code in [Rog00].

---

Modifier [-] stands for an algorithm defined in a previoous chapter;

Modifier [+] stands for an optional algorithm;

Modifier [?] stands for a missing pseudocode algorithm.

---

### Chapter 2: Bezier Curves

 - [X] [E] `bezier___` Calculates a Bezier curve using matrix method.
 - [X] [E] `dbezier__` Calculates a Bezier curve and it's derivatives.


### Chapter 3: B-spline Curves
 
 - [X] `basis____` Calculates the B-spline basis functions.
 - [X] [T] [D] `bsplfit__` [+] Subroutine to fit a B-spline curve using an open uniform knot vector. Call `basis, knot, param`.
 - [X] `bspline__` Calculates a B-spline curve. Call `basis, knot`.
 - [X] `bsplineu_` Calculates a periodic B-spline curve. Call `basis, knotu`.
 - [X] `dbasis___` Calculates the B-spline basis functions and derivatives.
 - [X] `dbasisu__` Calculates the periodic B-spline basis functions.
 - [X] `dbspline_` Calculates a B-spline curve and derivatives. Call `dbasis, knot`.
 - [X] `dbsplineu` Calculates a periodic B-spline curve and derivatives. Call `knotu, dbasisu`.
 - [X] `knot_____` Calculates an open knot vector.
 - [X] `knotc____` Calculates a chord length approximation open knot vector.
 - [X] `knotu____` Calculates a periodic knot vector.
 - [X] [T] [D] `matpbspl_` [+] Generate a B-spline curve using matrix methods and a periodic uniform knot vector. Call `nmatrix`.
 - [X] [T] [D] `nmatrix__` [+] Calculate the general B-spline periodic basis matrix.
 - [X] [T] [D] `param____` Calculates the chord length paramter values.
 - [X] `raise23` [+]
 - [X] `raise12` [+]
 - [ ] `remknot` [+]
 


### Chapter 4: Rational B-spline Curves (NURBC)

 - [X] `knot_____` [-] Calculates an open knot vector.
 - [X] `knotu____` [-] Calculates a periodic knot vector.
 - [ ] `rbasis___` Calculates the rational B-spline basis functions.
 - [ ] `rbspline_` Calculates an open rational B-spline curve. Call `rbasis, knot`.
 - [ ] `rbsplinu_` Calculates a periodic rational B-spline curve. Call `rbasis, knotu`.


### Chapter 5: Bezier Surface

 - [X] `bezsurf__` Calculates a Bezier surface using matrix method.

### Chapter 6: B-spline Surfaces

 - [X] `basis____` [-] Calculates the B-spline basis functions.
 - [X] [D] `bsplsurf_` Calculates a B-spline surface. Call `basis, knot`.
 - [X] [T] [D] `bspsurfu_` Calculates a periodic B-spline surface. Call `knotu, basis`.
 - [X] [D] `dbsurf___` Calculates the B-spline basis functions and derivatives. Call `dbasis, knot`.
 - [X] `dbasis___` [-] Calculates a B-spline surface and derivatives.
 - [X] `knot_____` [-] Calculates an open knot vector.
 - [X] `knotu____` [-] Calculates a periodic knot vector.


### Chapter 7: Rational B-spline (NURBS) Surfaces

 - [X] `basis____` [-] Calculates the B-spline basis functions.
 - [ ] `frbsurf__` [?] Calculates and test the fast B-spline surface algorithm. Call `basis, knot, sumrbas`
 - [X] `knot_____` [-] Calculates an open knot vector.
 - [X] `rbspsurf_` Calculates a rational B-spline (NURBS) surface. Call `knot, basis, sumrbas`.
 - [X] `sumrbas__` [+] Calculate the sum of the nonrational basis functions.
 - [ ] [ ] `rbsurf`

## Summary

 **Chap** | **No.** | **Test** | **Doc** 
----------|---------|----------|---------
 2        |     2   |     2    |     2
 3        |    17   |    12    |    12
 4        |     5   |     2    |     2
 5        |     1   |     0    |     1
 6        |     7   |     6    |     4
 7        |     6   |     4    |     4
 **Tot**  |    39   |    26    |    25
