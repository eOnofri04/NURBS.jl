var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#NURBS.jl-1",
    "page": "Home",
    "title": "NURBS.jl",
    "category": "section",
    "text": "NURBS.jl (or Non Uniform Rational B-Spline Surfaces) is a Julia library to build parametrical curves and parametrical surface based on the interpolation of polynomial bases. This library is developed by:Elia Onofri - Ln (elia.onofri4@gmail.com)\nGianmarco Caramitti (g.caramitti@gmail.com)\nPaolo Macciacchera (polmacra@outlook.it)\nGiuseppe Santorelli (giu.santorelli15@gmail.com)and would be maintained by the Computational Visual Design Laboratory (CVDLAB) of Università degli Studi di Roma Tre."
},

{
    "location": "index.html#Dependencies-1",
    "page": "Home",
    "title": "Dependencies",
    "category": "section",
    "text": "NURBS.jl has no dependencies but is build to work with:Plasm"
},

{
    "location": "index.html#Docstrings-conventions-1",
    "page": "Home",
    "title": "Docstrings conventions",
    "category": "section",
    "text": "Bold is used to point out theory concepts.\nMonospace is used for everything code related."
},

{
    "location": "nurbs.html#",
    "page": "NURBS Intro",
    "title": "NURBS Intro",
    "category": "page",
    "text": ""
},

{
    "location": "nurbs.html#A-short-introduction-to-NURBS-1",
    "page": "NURBS Intro",
    "title": "A short introduction to NURBS",
    "category": "section",
    "text": ""
},

{
    "location": "bezierC.html#",
    "page": "Bezier Curves",
    "title": "Bezier Curves",
    "category": "page",
    "text": ""
},

{
    "location": "bezierC.html#Bézier-Curves-1",
    "page": "Bezier Curves",
    "title": "Bézier Curves",
    "category": "section",
    "text": "Bézier Curves are a special kind of NURB curves build upon a controll polygon B and a Bernstein basis J_ni.Mathematically, a parametric Bézier curve is defined byP(t) = sum_i=0^n B_iJ_ni(t) quad 0 leq t leq 1Where the Bernstein Basis is defined asJ_n i(t) = binom ni t^i(1-t)^n-iwhere the convention (0)^0 = 1 and 0 = 1 have been made.In particular J_ni is the i-th base function of order n, while n is also the number of segments of the the controll polygon (number of points minus one)."
},

{
    "location": "bezierC.html#Bézier-Properties-1",
    "page": "Bezier Curves",
    "title": "Bézier Properties",
    "category": "section",
    "text": "Whe have a small variety of properties, descending directly from the definition:Base function are often real.\nThe degree of the polynomial curve is one less than the number of point of the controll polygon.\nThe curve generally follows the shape of the control polygon.\nThe edges of the curves are the edges of the control polygon.\nThe curve is located inside the convex hull of the controll polygon.\nThe curve is invariant under Affine transformation."
},

{
    "location": "bezierC.html#Matrix-Representation-1",
    "page": "Bezier Curves",
    "title": "Matrix Representation",
    "category": "section",
    "text": "The equation of a Bézier Curve could also be implemented as a Matrix Multiplication (particulary usefull in GPU computations).P(t) = FGwhereF = J_n0 dots J_n n	qquad	G^t = B_0 dots B_nMoreover is possible to collect the coefficient of the Basis in a square matrix N in mathcal M_nn(mathbb R) obtaining:P(t) = t^n t^n-1 dots t 1 N G = TNGwhereN = beginbmatrix\n	binom n0binom nn (-1)^n  binom n1 binomn-1n-1(-1)^n-1  dots  binom nn binomn-nn-n(-1)^0\n	binom n0binom nn-1 (-1)^n-1  binom n1 binomn-1n-2(-1)^n-2  dots  0\n	vdots  vdots  ddots  vdots\n	binom n0binom n1 (-1)^1  binom n1 binomn-10(-1)^0  dots  0\n	binom n0binom n0 (-1)^0  0  dots  0\nendbmatrixor, in other words:left(N_i+1j+1right)_ij=0^n = begincases\n		binom nj binomn-jn-i-j(-1)^n-i-j		mboxif  0 leq i+j leq n\n		0												mboxotherwise\n	endcasesIt is also possible to decompose the matrix N even further in the product of two matrices:N = CD quad Rightarrow quad P(t) = TCDGwhereC = beginbmatrix\n	binom nn (-1)^n  binomn-1n-1(-1)^n-1  dots  binomn-nn-n(-1)^0\n	binom nn-1 (-1)^n-1  binomn-1n-2(-1)^n-2  dots  0\n	vdots  vdots  ddots  vdots\n	binom n1 (-1)^1  binomn-10(-1)^0  dots  0\n	binom n0 (-1)^0  0  dots  0\nendbmatrixD = beginbmatrix\n	binom n0  0  dots  0\n	0  binom n1  dots  0\n	vdots  vdots  ddots  vdots\n	0  0  dots  binom nn\nendbmatrix"
},

{
    "location": "bezierC.html#Bézier-Derivatives-1",
    "page": "Bezier Curves",
    "title": "Bézier Derivatives",
    "category": "section",
    "text": "The two derivatives could be obtained starting from the original function:P(t) = sum_i=0^n B_iJ_ni(t)P(t) = sum_i=0^n B_iJ_ni(t)where the two derivatives could be obtained with the following formulas:J_ni(t) = frac1-ntt(1-t)J_ni(t)J_ni(t) = frac(i-nt)^2-nt^2-i(1-2t)t^2(1-t)^2 J_ni(t)defined in all the points but the first and the last (where t = 0 and t = 1)."
},

{
    "location": "bezierC.html#Derivatives-Matrix-Representation-1",
    "page": "Bezier Curves",
    "title": "Derivatives Matrix Representation",
    "category": "section",
    "text": "In order to enhance the matrix representation it is usefull to represent the derivatives in two matrices:Der1 = beginbmatrix\n	frac-nt_(1)t_(1)*(1-t_(1))		dots		dots		dots		fracn-nt_(0)t_(0)(1-t_(0))	\n\n	vdots		ddots						vdots							vdots		vdots	\n	dots		dots		frac(j-1)-nt_(i)t_(i)(1-t_(i))		dots		dots	\n	vdots		vdots						vdots							ddots		vdots	\n\n	frac-nt_(bpts-1)t_(bpts-1)(1-t_(bpts-1))		dots		dots		dots		fracn-nt_(bpts-1)t_(bpts-1)(1-t_(bpts-1)) \nendbmatrixandDer2 = beginbmatrix\n	frac(-nt_(1))^2- nt_(1)^2t_(1)^2*(1-t_(1))^2		dots		dots		dots		frac(n-nt_(0))^2-nt_(0)^2-n(1-2t_(0))t_(0)^2(1-t_(0))^2	\n\n	vdots		ddots										vdots												vdots		vdots	\n	dots		dots		frac((j-1)-nt_(i))^2-nt_(i)^2-(j-1)(1-2t_(i))t_(i)^2(1-t_(i))^2		dots		dots	\n	vdots		vdots										vdots												ddots		vdots	\n\n	frac(-nt_(bpts-1))^2-nt_(bpts-1)^2t_(bpts-1)^2(1-t_(bpts-1))^2		dots		dots		dots		frac(n-nt_(bpts-1))^2-nt_(bpts-1)^2-n(1-2t_(bpts-1))t_(bpts-1)^2(1-t_(bpts-1))^2 \nendbmatrixwhich could be scalar multiplied with the base matrix in order to obtain two derivative base matrices. Of course the T matrix used must be reduced by eliminating the first and the last points.The result of the operation is then:beginsplit\n	P(t) = 	Der1cdot(T_1^dpts-1CD)G \n	P(t) = 	Der2cdot(T_1^dpts-1CD)G"
},

{
    "location": "splineC.html#",
    "page": "B-Spline Curves",
    "title": "B-Spline Curves",
    "category": "page",
    "text": ""
},

{
    "location": "splineC.html#B-Spline-Curves-1",
    "page": "B-Spline Curves",
    "title": "B-Spline Curves",
    "category": "section",
    "text": "A B-Spline is the representation of a curve (i.e a function) build from the interpolation of the elements of a normalized base of the analytic functions space.The interpolation is build between a vector of points (called knots) positioned inside a controll polygon (delimited by a family of controll points B).The positioning vector (i.e. the parametric funtion) of a B-Spline is defined as follows:P(t) = sum_i=1^n+1B_iN_ik(t) qquad 2 leq k leq n+1where:B_i is the i-th point between the n+1 controll points (of the polygon)\nN_ik is the function of the i-th basis of the B-Spline normalized with order k (i.e degree k-1).A method to evaluate N_i k is given by +Cox-de Boor* recursive form.    N_ik(t) = frac(t-x_i)N_ik-1(t)x_i+k-1-x_i + frac(x_i+k-t)N_i+1k-1(t)x_i+k-x_i+1\n    qquad mboxcon quad\n    N_i1(t) = \n    begincases\n        1    mbox if  x_i  t  x_i+1\n        0    mbox otherwise\n    endcases"
},

{
    "location": "splineC.html#B-Spline-Properties-1",
    "page": "B-Spline Curves",
    "title": "B-Spline Properties",
    "category": "section",
    "text": "The following are properties hold by B-Splines:The sum of the basis function in every point t is equal to one (i.e. sum_i=1^n+1 N_ik = 1). \nBasis function are non-negative for each and every point t.\nThe order of a curve is at most equals to the number n+1 of controll points (so the maximum degree is n).\nA curve could be modified by an affine function f by applying the function to the controll points of the polygon.\nThe curve is located inside the convex hull of the controll polygon."
},

{
    "location": "splineC.html#Knots-1",
    "page": "B-Spline Curves",
    "title": "Knots",
    "category": "section",
    "text": "Knots choise is very important.There is an important relation between the number of controll polygon points m, the order of the function k and the number of knots m, wich is:m = k + n +1We have two kind of knots:periodics: The first and the last value has k multiplicity;\nopen: Each value has the same multiplicity.wich could be build in two manners:uniform: Knots are evenly spaced;\nnon uniform: there are different space between knots.We have so four classes of knots:Uniform Periodics: Wich are linked to a base function in the following form: N_ik(t) = N_i-1k(t-1) = N_i+1k(t+1)\nOpen Uniform: They have an even space between knots and the multiplicity is k at the edges. Follow from this definition that if the number of controll polygon points is equals to the order of the curve, than this curve is a Bézier curve: in fact the basis is build just like a Bernstein Basis. An example is given by k = 3 qquad 0 0 0 1 2 3 3 3.\nOpen Non-Uniform\nPeriodic Non-Uniform"
},

{
    "location": "splineC.html#B-Spline-Basis-Functions-1",
    "page": "B-Spline Curves",
    "title": "B-Spline Basis Functions",
    "category": "section",
    "text": "As the function N_ik is defined by Cox-de Boor formulas in a recursive way, the evaluation of a basis set could be optimized by saving the previous evaluation. The dependency tree is:    beginmatrix\n        N_i k\n        N_i k-1  N_i+1 k-1\n        N_i k-2  N_i+1 k-2  N_i+2 k-2\n        vdots  vdots  vdots  ddots\n        N_i 1  N_i+1 1  N_i+2 1  dots  N_i+k-1 1\n    endmatrix\n    beginmatrix\n        N_i-k+1 k  dots  N_i-1 k  N_i k  N_i+1 k  dots  N_i+k-1 k\n                     ddots  vdots      vdots    vdots                         \n                             N_i-1 2  N_i 2  N_i+1 2 \n                                         N_i 1 \n    endmatrix"
},

{
    "location": "rSplineC.html#",
    "page": "Rational B-Spline Curves",
    "title": "Rational B-Spline Curves",
    "category": "page",
    "text": ""
},

{
    "location": "rSplineC.html#Rational-B-Spline-Curves-1",
    "page": "Rational B-Spline Curves",
    "title": "Rational B-Spline Curves",
    "category": "section",
    "text": ""
},

{
    "location": "bezierS.html#",
    "page": "Bezier Surfaces",
    "title": "Bezier Surfaces",
    "category": "page",
    "text": ""
},

{
    "location": "bezierS.html#Bézier-Surfaces-1",
    "page": "Bezier Surfaces",
    "title": "Bézier Surfaces",
    "category": "section",
    "text": "Bézier Surfaces are special kind of NURBS build upon a controll net B and a couple of Bernstein basis J_ni(u) and K_mj(w).Mathematically, a parametric Bézier surface is defined byQ(u w) = sum_i=0^n sum_j=0^m B_ij J_ni(u) K_mj(w)where u and w are two parametric directions and n = npts - 1 while m = mpts - 1.We point out that the Bernstein Basis are defined asbeginsplit\n	J_n i(u) = binom ni u^i(1-u)^n-i\n	K_m j(w) = binom mj w^j(1-w)^m-j\nendsplit"
},

{
    "location": "bezierS.html#Bézier-Properties-1",
    "page": "Bezier Surfaces",
    "title": "Bézier Properties",
    "category": "section",
    "text": "Many properties of Bézier Surfaces are known:The degree of the surface in each parametric direction is one less than the number of control net vertices in that direction.\nThe continuity of the surface in each parametric direction is two less than the number of control net vertices in that direction.\nThe surface generally follows the shape of the control net.\nOnly the corner points of the control net and the resulting Bézier surface are coincident.\nThe surface is contained within the convex hull of the control net.\nThe surface does not exhibit the variation-diminishing property. The variation-diminishing property for bivariant surfaces is both undefined and unknown.\nThe surface is invariant under an affine transformation."
},

{
    "location": "bezierS.html#Matrix-Representation-1",
    "page": "Bezier Surfaces",
    "title": "Matrix Representation",
    "category": "section",
    "text": "Like in Bézier Curves, also Bézier surfaces are very well described by a matrix point of view. Here we could simplify the calculus by writing:	Q(u w) = UNBM^TW^Twherebeginsplit\n	U = beginbmatrix \n		u^n  u^n-1  dots  1\n	endmatrix\n	W = beginbmatrix \n		w^n  u^w-1  dots  1\n	endmatrix\n	B = beginbmatrix \n		B_00		dots		B_0m	\n		vdots		ddots		vdots	\n		B_n0		dots		B_nm\n	endmatrix\nendsplitand N and M are the coefficient matrix of the base, so given by:	beginsplit\n	N = beginbmatrix\n		binom n0binom nn (-1)^n  binom n1 binomn-1n-1(-1)^n-1  dots  binom nn binomn-nn-n(-1)^0\n		binom n0binom nn-1 (-1)^n-1  binom n1 binomn-1n-2(-1)^n-2  dots  0\n		vdots  vdots  ddots  vdots\n		binom n0binom n1 (-1)^1  binom n1 binomn-10(-1)^0  dots  0\n		binom n0binom n0 (-1)^0  0  dots  0\n	endbmatrix \n	M = beginbmatrix\n		binom m0binom mm (-1)^m  binom m1 binomm-1m-1(-1)^m-1  dots  binom mm binomm-mm-m(-1)^0\n		binom m0binom mm-1 (-1)^m-1  binom m1 binomm-1m-2(-1)^m-2  dots  0\n		vdots  vdots  ddots  vdots\n		binom m0binom m1 (-1)^1  binom m1 binomm-10(-1)^0  dots  0\n		binom m0binom m0 (-1)^0  0  dots  0\n	endbmatrix or, in other worlds:beginsplit\n	left(N_i+1j+1right)_ij=0^n = begincases\n		binom nj binomn-jn-i-j(-1)^n-i-j		mboxif  0 leq i+j leq n\n		0												mboxotherwise\n	endcases\n	left(M_i+1j+1right)_ij=0^m = begincases\n		binom mj binomm-jm-i-j(-1)^m-i-j		mboxif  0 leq i+j leq m\n		0												mboxotherwise\n	endcases\nendsplitof course it is possible to decompose N and M in two more matrices, like for curves, leading to the following scheme:Q(u w) = U(CD)B(EF)^TW^T = UCDBF^TE^TW^TMuch more interesting is the situation where m = n as M = N"
},

{
    "location": "bezierS.html#Bézier-Surface-Derivatives-1",
    "page": "Bezier Surfaces",
    "title": "Bézier Surface Derivatives",
    "category": "section",
    "text": "Of course is simple to determine the derivatives starting from the matrix notation.We recall th results about the first and the second derivatives but the process could be extended further with no particular changes:beginsplit\n	frac partialpartial u Q(uw) = 	UNBM^TW^T\n	frac partialpartial w Q(uw) = 	UNBM^TW^T\n	frac partial^2partial upartial w Q(uw) = 	UNBM^TW^T\n	frac partial^2partial u^2 Q(uw) = 	UNBM^TW^T\n	frac partial^2partial w^2 Q(uw) = 	UNBM^TW^T\nendsplitwhere, as expected, we have:beginsplit\n	U = nu^n-1  (n-1)u^n-2  dots  2n  1  0\n	W = mw^m-1  (m-1)w^m-2  dots  2m  1  0\n	U = n(n-1)u^n-2  (n-1)(n-2)u^n-3  dots  1  0  0\n	W = m(m-1)w^m-2  (m-1)(m-2)w^m-3  dots  1  0  0\nendsplit"
},

{
    "location": "splineS.html#",
    "page": "B-Spline Surfaces",
    "title": "B-Spline Surfaces",
    "category": "page",
    "text": ""
},

{
    "location": "splineS.html#B-Spline-Surfaces-1",
    "page": "B-Spline Surfaces",
    "title": "B-Spline Surfaces",
    "category": "section",
    "text": ""
},

{
    "location": "rSplineS.html#",
    "page": "Rational B-Spline Surfaces",
    "title": "Rational B-Spline Surfaces",
    "category": "page",
    "text": ""
},

{
    "location": "rSplineS.html#Rational-B-Spline-Surfaces-1",
    "page": "Rational B-Spline Surfaces",
    "title": "Rational B-Spline Surfaces",
    "category": "section",
    "text": ""
},

]}
