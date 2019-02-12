# B-Spline Surfaces

the exstention of Bézier surface are the B-spline surface defined by

"Q(u,w)=sum_{i=1}^{n+1} sum_{j=1}^{m+1} B_{i,j} N_{i,k}(u) M_{j,l}(w)"

where $N_{i,k}(u)$ and M_{j,l}(w) are the B-spline basis functions in the biparametric u and w directions.
The definition for the basis functions is given by

begin{split}
     
    N_{i,1}(u)= \left{ begin{split}
                            1 \quad if x_i \leq u < x_{i+1}
                            0 otherwise
                        end{split}
    N_{i,k}(u) \frac{(u-x_i)N_{i,k-1}(u)}{x_{i+k-1}-x_i} + \frac{(x_{i+k}-u)N_{i+1,k-1}(u)}{x_{i+k}-x_{i+1}}
    
end{split}
and
begin{split}
     
    M_{j,1}(w)= \left{ begin{split}
                            1 \quad if y_i \leq w < y_{j+1}
                            0 otherwise
                        end{split}
    M_{j,l}(w) \frac{(w-y_j)M_{j,l-1}(w)}{y_{j+l-1}-y_j} + \frac{(y_{j+l}-w)N_{j+1,l-1}(w)}{y_{j+l}-x_{j+1}}
    
end{split}

where the x_i and y_j are elements of knot vectors.
Again, the B_{i,j}s are the vertices of a polygonal control net and
the indices n and m are one less than the number of control polygon
vertices.The indice n is in the parametric direction of u, respectively m in w.

                              

