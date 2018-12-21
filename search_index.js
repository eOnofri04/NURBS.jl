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
    "text": ""
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
    "text": "A B-Spline is the representation of a curve (_i.e_ a function) build from the interpolation of the elements of a normalized base of the analytic functions space.The interpolation is build between a vector of points (called _knots_) positioned inside a controll polygon (delimited by a family of controll points B).The positioning vector (_i.e._ the parametric funtion) of a B-Spline is defined as follows: $     P(t) = \\sum_{i=1}^{n+1}B_iN_{i,k}(t), \\qquad 2 \\leq k \\leq n+1 $ where:B_i is the i-th point between the n+1 controll points (of the polygon)\nN_{i,k} is the function of the i-th basis of the B-Spline normalized with order k (_i.e_ degree k-1).A method to evaluate N_{i, k} is given by _Cox-de Boor_ recursive form. $     N_{i,k}(t) = \\frac{(t-x_i)N_{i,k-1}(t)}{x_{i+k-1}-x_i} + \\frac{(x_{i+k}-t)N_{i+1,k-1}(t)}{x_{i+k}-x_{i+1}},     \\qquad \\mbox{con} \\quad     N_{i,1}(t) =      \\begin{cases}         1  &  \\mbox{ if } x_i < t < x_{i+1}\\\n        0  &  \\mbox{ otherwise}     \\end{cases} $"
},

{
    "location": "splineC.html#B-Spline-Properties-1",
    "page": "B-Spline Curves",
    "title": "B-Spline Properties",
    "category": "section",
    "text": "The following are properties hold by B-Splines:The sum of the basis function in every point t is equal to one:$      \\sum_{i=1}^{n+1} N_{i,k} = 1  $Basis function are non-negative for each and every point t\nThe order of a curve is at most equals to the number n+1 of controll points (so the maximum degree is n).\nA curve could be modified by an affine function f by applying the function to the controll points of the polygon.\nThe curve is located inside the convex hull of the controll polygon."
},

{
    "location": "splineC.html#Knots-1",
    "page": "B-Spline Curves",
    "title": "Knots",
    "category": "section",
    "text": "Knots choise is very important.There is an important relation between the number of controll polygon points m, the order of the function k and the number of knots m, wich is: $     m = k + n +1 $We have two kind of knots:periodics: The first and the last value has k multiplicity;\nopen: Each value has the same multiplicity.wich could be build in two manners:uniform: Knots are evenly spaced;\nnon uniform: there are different space between knots.We have so four classes of knots:Uniform Periodics: Wich are linked to a base d $     N_{i,k}(t) = N_{i-1,k}(t-1) = N_{i+1,k}(t+1) $\nOpen Uniform: They have an even space between knots and the multiplicity is k at the edges, for example: $     k = 3 \\qquad [0\\ 0\\ 0\\ 1\\ 2\\ 3\\ 3\\ 3] $ Follow from this definition that if the number of controll polygon points is equals to the order of the curve, than this curve is a Bézier curve: in fact the basis is build just like a Bernstein Basis.\nOpen Non-Uniform\nPeriodic Non-Uniform"
},

{
    "location": "splineC.html#Funzioni-delle-Basi-per-le-B-Spline-1",
    "page": "B-Spline Curves",
    "title": "Funzioni delle Basi per le B-Spline",
    "category": "section",
    "text": "As the function N_{i,k} is defined by _Cox-de Boor_ formulas in a recursive way, the evaluation of a basis set could be optimized by saving the previous evaluation. The dependency tree is: $     \\begin{matrix}         N_{i, k}\\\n        N_{i, k-1} & N_{i+1, k-1}\\\n        N_{i, k-2} & N_{i+1, k-2} & N_{i+2, k-2}\\\n        \\vdots & \\vdots & \\vdots & \\ddots\\\n        N_{i, 1} & N_{i+1, 1} & N_{i+2, 1} & \\dots & N_{i+k-1, 1}     \\end{matrix}     \\begin{matrix}         N_{i-k+1, k} & \\dots & N_{i-1, k} & N_{i, k} & N_{i+1, k} & \\dots & N_{i+k-1, k}\\\n                    & \\ddots & \\vdots     & \\vdots   & \\vdots     &       &             \\\n                     &       & N_{i-1, 2} & N_{i, 2} & N_{i+1, 2} \\\n                     &       &            & N_{i, 1}      \\end{matrix} $"
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
