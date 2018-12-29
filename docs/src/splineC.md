# B-Spline Curves


A B-Spline is the representation of a curve (*i.e* a function) build from the interpolation of the elements of a normalized base of the analytic functions space.

The interpolation is build between a vector of points (called *knots*) positioned inside a controll polygon (delimited by a family of controll points `B`).

The positioning vector (*i.e.* the parametric funtion) of a B-Spline is defined as follows:
``P(t) = \sum_{i=1}^{n+1}B_iN_{i,k}(t), \qquad 2 \leq k \leq n+1``
where:
 - `B_i` is the `i`-th point between the `n+1` controll points (of the polygon)
 - `N_{i,k}` is the function of the `i`-th basis of the B-Spline normalized with order `k` (*i.e* degree `k-1`).

A method to evaluate `N_{i, k}` is given by +Cox-de Boor* recursive form.
``
    N_{i,k}(t) = \frac{(t-x_i)N_{i,k-1}(t)}{x_{i+k-1}-x_i} + \frac{(x_{i+k}-t)N_{i+1,k-1}(t)}{x_{i+k}-x_{i+1}},
    \qquad \mbox{con} \quad
    N_{i,1}(t) = 
    \begin{cases}
        1  &  \mbox{ if } x_i < t < x_{i+1}\\
        0  &  \mbox{ otherwise}
    \end{cases}
``

---
---
## B-Spline Properties

The following are properties hold by B-Splines:
 - The sum of the basis function in every point `t` is equal to one:
  ``\sum_{i=1}^{n+1} N_{i,k} = 1``
 - Basis function are non-negative for each and every point `t`
 - The order of a curve is at most equals to the number `n+1` of controll points (so the maximum degree is `n`).
 - A curve could be modified by an affine function `f` by applying the function to the controll points of the polygon.
 - The curve is located inside the convex hull of the controll polygon.

---
---
## Knots

Knots choise is very important.

There is an important relation between the number of controll polygon points `m`, the order of the function `k` and the number of knots `m`, wich is:
``m = k + n +1``

We have two kind of knots:
 - **periodics**: The first and the last value has `k` multiplicity;
 - **open**: Each value has the same multiplicity.

wich could be build in two manners:
 - **uniform**: Knots are evenly spaced;
 - **non uniform**: there are different space between knots.

We have so four classes of knots:
 - **Uniform Periodics**: Wich are linked to a base d
   ``N_{i,k}(t) = N_{i-1,k}(t-1) = N_{i+1,k}(t+1)``
 - **Open Uniform**: They have an even space between knots and the multiplicity is `k` at the edges, for example:
   ``k = 3 \qquad [0\ 0\ 0\ 1\ 2\ 3\ 3\ 3]``
   Follow from this definition that if the number of controll polygon points is equals to the order of the curve, than this curve is a Bézier curve: in fact the basis is build just like a Bernstein Basis.
 - **Open Non-Uniform**
 - **Periodic Non-Uniform**

---
---
## B-Spline Basis Functions

As the function `N_{i,k}` is defined by *Cox-de Boor* formulas in a recursive way, the evaluation of a basis set could be optimized by saving the previous evaluation. The dependency tree is:
``
    \begin{matrix}
        N_{i, k}\\
        N_{i, k-1} & N_{i+1, k-1}\\
        N_{i, k-2} & N_{i+1, k-2} & N_{i+2, k-2}\\
        \vdots & \vdots & \vdots & \ddots\\
        N_{i, 1} & N_{i+1, 1} & N_{i+2, 1} & \dots & N_{i+k-1, 1}
    \end{matrix}
    \begin{matrix}
        N_{i-k+1, k} & \dots & N_{i-1, k} & N_{i, k} & N_{i+1, k} & \dots & N_{i+k-1, k}\\
                    & \ddots & \vdots     & \vdots   & \vdots     &       &             \\
                     &       & N_{i-1, 2} & N_{i, 2} & N_{i+1, 2} \\
                     &       &            & N_{i, 1} 
    \end{matrix}
``
