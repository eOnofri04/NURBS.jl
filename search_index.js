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
    "text": "NURBS.jl (or Non Uniform Rational B-Spline) is a Julia library to build parametrical curves and parametrical surface based on the interpolation of polynomial bases. This library is developed by:Elia Onofri - Ln (elia.onofri4@gmail.com)\nGianmarco Caramitti (g.caramitti@gmail.com)\nPaolo Macciacchera (polmacra@outlook.it)\nGiuseppe Santorelli (giu.santorelli15@gmail.com)and would be maintained by the Computational Visual Design Laboratory (CVDLAB) of Università degli Studi di Roma Tre."
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
    "location": "bezierC.html#Bezier-Curves-1",
    "page": "Bezier Curves",
    "title": "Bezier Curves",
    "category": "section",
    "text": "Bezier Curves are a special kind of NURBS curves build upon a controll polygon B and a Bernstein basis J_ni.Mathematically, a parametric Bézier curve is defined byP(t) = sum_i=0^n B_iJ_ni(t) quad 0 leq t leq 1Where the Bernstein Basis is defined asJ_n i(t) = binom ni t^i(1-t)^n-iwhere the convention (0)^0 = 1 and 0 = 1 have been made.In particular J_ni is the i-th base function of order n, while n is also the number of segments of the the controll polygon (number of points minus one)."
},

{
    "location": "bezierC.html#Bezier-Properties-1",
    "page": "Bezier Curves",
    "title": "Bezier Properties",
    "category": "section",
    "text": "Whe have a small variety of properties, descending directly from the definition:Base function are often real.\nThe degree of the polynomial curve is one less than the number of point of the controll polygon.\nThe curve generally follows the shape of the control polygon.\nThe edges of the curves are the edges of the control polygon.\nThe curve is located inside the convex hull of the controll polygon.\nThe curve is invariant under Affine transformation."
},

{
    "location": "bezierC.html#Matrix-Representation-1",
    "page": "Bezier Curves",
    "title": "Matrix Representation",
    "category": "section",
    "text": "The equation of a Bézier Curve could also be implemented as a Matrix Multiplication (particulary usefull in GPU computations).P(t) = FGwhereF = J_n0 dots J_n n	qquad	G^t = B_0 dots B_nMoreover is possible to collect the coefficient of the Basis in a square matrix N in mathcal M_nn(mathbb R) obtaining:P(t) = t^n t^n-1 dots t 1 N G = TNGwhereN = beginbmatrix\n	binom n0binom nn (-1)^n  binom n1 binomn-1n-1(-1)^n-1  dots  binom nn binomn-nn-n(-1)^0\n	binom n0binom nn-1 (-1)^n-1  binom n1 binomn-1n-2(-1)^n-2  dots  0\n	vdots  vdots  ddots  vdots\n	binom n0binom n1 (-1)^1  binom n1 binomn-10(-1)^0  dots  0\n	binom n0binom n0 (-1)^0  0  dots  0\nendbmatrixIt is also possible to decompose the matrix N even further in the product of two matrices:N = CD quad Rightarrow quad P(t) = TCDGwhereN = beginbmatrix\n	binom nn (-1)^n  binomn-1n-1(-1)^n-1  dots  binomn-nn-n(-1)^0\n	binom nn-1 (-1)^n-1  binomn-1n-2(-1)^n-2  dots  0\n	vdots  vdots  ddots  vdots\n	binom n1 (-1)^1  binomn-10(-1)^0  dots  0\n	binom n0 (-1)^0  0  dots  0\nendbmatrixN = beginbmatrix\n	binom n0  0  dots  0\n	0  binom n1  dots  0\n	vdots  vdots  ddots  vdots\n	0  0  dots  binom nn\nendbmatrix"
},

{
    "location": "bezierC.html#Bezier-Derivatives-1",
    "page": "Bezier Curves",
    "title": "Bezier Derivatives",
    "category": "section",
    "text": "The two derivatives could be obtained starting from the original function:P(t) = sum_i=0^n B_iJ_ni(t)P(t) = sum_i=0^n B_iJ_ni(t)where the two derivatives could be obtained with the following formulas:J_ni(t) = frac1-ntt(1-t)J_ni(t)J_ni(t) = frac(i-nt)^2-nt^2-i(1-2t)t^2(1-t)^2 J_ni(t)"
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
    "text": "The following are properties hold by B-Splines:The sum of the basis function in every point t is equal to one:sum_i=1^n+1 N_ik = 1Basis function are non-negative for each and every point t\nThe order of a curve is at most equals to the number n+1 of controll points (so the maximum degree is n).\nA curve could be modified by an affine function f by applying the function to the controll points of the polygon.\nThe curve is located inside the convex hull of the controll polygon."
},

{
    "location": "splineC.html#Knots-1",
    "page": "B-Spline Curves",
    "title": "Knots",
    "category": "section",
    "text": "Knots choise is very important.There is an important relation between the number of controll polygon points m, the order of the function k and the number of knots m, wich is:m = k + n +1We have two kind of knots:periodics: The first and the last value has k multiplicity;\nopen: Each value has the same multiplicity.wich could be build in two manners:uniform: Knots are evenly spaced;\nnon uniform: there are different space between knots.We have so four classes of knots:Uniform Periodics: Wich are linked to a base d N_ik(t) = N_i-1k(t-1) = N_i+1k(t+1)\n \n - **Open Uniform** They have an even space between knots and the multiplicity is k at the edges for example\n math\n k = 3 qquad 0 0 0 1 2 3 3 3\n \n   Follow from this definition that if the number of controll polygon points is equals to the order of the curve than this curve is a Bzier curve in fact the basis is build just like a Bernstein Basis\n - **Open Non-Uniform**\n - **Periodic Non-Uniform**\n\n---\n---\n B-Spline Basis Functions\n\nAs the function N_ik is defined by *Cox-de Boor* formulas in a recursive way the evaluation of a basis set could be optimized by saving the previous evaluation The dependency tree ismath     \\begin{matrix}         N_{i, k}\\\n        N_{i, k-1} & N_{i+1, k-1}\\\n        N_{i, k-2} & N_{i+1, k-2} & N_{i+2, k-2}\\\n        \\vdots & \\vdots & \\vdots & \\ddots\\\n        N_{i, 1} & N_{i+1, 1} & N_{i+2, 1} & \\dots & N_{i+k-1, 1}     \\end{matrix}     \\begin{matrix}         N_{i-k+1, k} & \\dots & N_{i-1, k} & N_{i, k} & N_{i+1, k} & \\dots & N_{i+k-1, k}\\\n                    & \\ddots & \\vdots     & \\vdots   & \\vdots     &       &             \\\n                     &       & N_{i-1, 2} & N_{i, 2} & N_{i+1, 2} \\\n                     &       &            & N_{i, 1}      \\end{matrix} ```"
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
    "location": "bezierS.html#Bezier-Surfaces-1",
    "page": "Bezier Surfaces",
    "title": "Bezier Surfaces",
    "category": "section",
    "text": ""
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
