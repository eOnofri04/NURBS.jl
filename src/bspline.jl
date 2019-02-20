export bspline, bsplineu, dbspline, dbsplineu, matpbspl, nmatrix


"""
	bspline(npts, ord, p1, b[])

Generate a B-spline curve using an open uniform knot vector.

---

Returns a tuple `(P,EV)`. P is the coordinates matrix of the p1 points b-spline curve defined by the control points matrix b. EV is the 1-dimensional cellular complex.

# Arguments

- `npts::Int64`: the number of the contol polygon vertices.
- `ord::Int64`: order of the bspline basis function.
- `p1::Int64`: number of the points to be calculated on the curve.
- `b::Array{Float64, 2}`: 2-dimensional Array containing the control polygon vertices in the x, y, z coordinates for columns.

# Example

```jldoctest
julia> bspline(7,3,20,[0. 1 2 3 4 5 6; 1 2 1 2 1 2 3; 0 0 0 0 0 0 0])[1]
3×20 Array{Float64,2}:
 0.0  0.49169  0.914127  1.26731  1.55263  …  4.73269  5.08587  5.50831  6.0
 1.0  1.42244  1.63712   1.64404  1.45014     1.73269  2.08587  2.50831  3.0
 0.0  0.0      0.0       0.0      0.0         0.0      0.0      0.0      0.0

julia> bspline(7,3,20,[0. 1 2 3 4 5 6; 1 2 1 2 1 2 3; 0 0 0 0 0 0 0])[2]
19-element Array{Array{Int64,1},1}:
 [1, 2]  
 [2, 3]  
 [3, 4]  
 [4, 5]  
 [5, 6]  
 [6, 7]  
 [7, 8]  
 [8, 9]  
 [9, 10] 
 [10, 11]
 [11, 12]
 [12, 13]
 [13, 14]
 [14, 15]
 [15, 16]
 [16, 17]
 [17, 18]
 [18, 19]
 [19, 20]
```
---
_By Giuseppe Santorelli, Elia Onofri_
"""

function bspline(npts::Int64, ord::Int64, p1::Int64, b::Array{Float64,2})::Tuple{Array{Float64,2},Array{Array{Int64,1},1}}
    
    @assert length(b) == 3 * npts ("ERROR: array b not matching with parameters")
    @assert npts >= ord ("ERROR: ord > npts")

    m::Int64 = npts + ord    
    x::Array{Float64} = knot(npts,ord) #creates knot array
    
    G::Array{Float64} = transpose(b)
    
    t::Float64 = 0.0
    step::Float64 = x[m] / (p1-1)
            
    N::Array{Float64} = zeros(p1, npts)
       
    for i = 1 : p1 #computes the value of the B-spline in the p1 points
        if (x[m] - t) < 5e-6
            t = x[m]
        end
        nbasis::Array{Float64} = basis(npts,ord,t,x) #computes the basis value depending on patameter t
        for j = 1 : npts
            N[i,j] = nbasis[j]
        end
        t = t + step
    end
    
    P::Array{Float64} = transpose(N*G)
    
    EV = Array{Int64,1}[] #1-dimensional cellular complex used for plotting the curve with Plasm package

    for i = 1 : p1-1
        C = Array{Int64}(2)
        C[1] = i
        C[2] = i+1
        push!(EV,C) 
    end
    
    return (P,EV)
end

#-------------------------------------------------------------------------------------

"""
	bsplineu(npts, ord, p1, b[])


Generate a B-spline curve using a periodic uniform knot vector.

---

Returns a tuple `(P,EV)`. P is the coordinates matrix of the p1 points b-spline curve defined by the control points matrix b using an open knot vector. EV is the 1-dimensional cellular complex.

# Arguments

- `npts::Int64`: the number of the contol polygon vertices.
- `ord::Int64`: order of the bspline basis function.
- `p1::Int64`: number of the points to be calculated on the curve.
- `b::Array{Float64, 2}`: 2-dimensional Array containing the control polygon vertices in the x, y, z coordinates for columns.

# Example

```jldoctest
julia> bsplineu(4,3,20,[0. 3 6 9; 0 10 3 6; 0 0 0 0])[1]
3×20 Array{Float64,2}:
 1.5  1.81579  2.13158  2.44737  2.76316  …  6.55263  6.86842  7.18421  7.5
 5.0  5.95845  6.72853  7.31025  7.7036      4.05125  4.09003  4.23961  4.5
 0.0  0.0      0.0      0.0      0.0         0.0      0.0      0.0      0.0

julia> bspline(4,3,20,[0. 3 6 9; 0 10 3 6; 0 0 0 0])[2]
19-element Array{Array{Int64,1},1}:
 [1, 2]  
 [2, 3]  
 [3, 4]  
 [4, 5]  
 [5, 6]  
 [6, 7]  
 [7, 8]  
 [8, 9]  
 [9, 10] 
 [10, 11]
 [11, 12]
 [12, 13]
 [13, 14]
 [14, 15]
 [15, 16]
 [16, 17]
 [17, 18]
 [18, 19]
 [19, 20]
```
---
_By Giuseppe Santorelli, Elia Onofri_
"""

function bsplineu(npts::Int64,ord::Int64,p1::Int64,b::Array{Float64,2})::Tuple{Array{Float64,2},Array{Array{Int64,1},1}}
    
    @assert length(b) == 3 * npts ("ERROR: array b not matching with parameters")
    @assert npts >= ord ("ERROR: ord > npts")

    m::Int64 = npts + ord    
    x::Array{Float64} = knotu(npts,ord) #creates knot array
    
    G::Array{Float64} = transpose(b)
    
    t::Float64 = ord - 1
    step::Float64 = (npts - (ord - 1)) / (p1 - 1)
    
    N::Array{Float64} = zeros(p1, npts)
    
    for i = 1 : p1 #computes the value of the B-spline in the p1 points
        if (npts - t) < 5e-6
            t = npts
        end
        
        nbasis::Array{Float64} = basis(npts,ord,t,x) #computes the basis value depending on patameter t
        
        for j = 1 : npts
            N[i,j] = nbasis[j]
        end
        t = t + step
        
    end
    
    P::Array{Float64} = transpose(N*G)
        
    EV = Array{Int64,1}[] #1-dimensional cellular complex used for plotting the curve with Plasm package

    for i = 1 : p1-1
        C = Array{Int64}(2)
        C[1] = i
        C[2] = i+1
        push!(EV,C) 
    end
        
    return (P,EV)
    
end


#------------------------------------------------------------------------------------


"""
	dbspline(npts, k, p1, b[])

Generate a B-spline curve and its derivatives using an open uniform knot vector.

Returns a tuple `(P,D1,D2,EV)`. P is the coordinates matrix of the p1 points b-spline curve and defined by the control points matrix b. D1 and D2 are the first and second derivatives of P. EV is the 1-dimensional cellular complex. 

---

# Arguments

- `npts::Int64`: the number of the contol polygon vertices.
- `ord::Int64`: order of the bspline basis function.
- `p1::Int64`: number of the points to be calculated on the curve.
- `b::Array{Float64, 2}`: 2-dimensional Array containing the control polygon vertices in the x, y, z coordinates for columns.

# Example

```jldoctest
julia> dbspline(7,3,20,[0. 1 2 3 4 5 6; 1 2 1 2 1 2 3; 0 0 0 0 0 0 0])[1]
3×20 Array{Float64,2}:
 0.0  0.49169  0.914127  1.26731  1.55263  …  4.73269  5.08587  5.50831  6.0
 1.0  1.42244  1.63712   1.64404  1.45014     1.73269  2.08587  2.50831  3.0
 0.0  0.0      0.0       0.0      0.0         0.0      0.0      0.0      0.0

julia> dbspline(7,3,20,[0. 1 2 3 4 5 6; 1 2 1 2 1 2 3; 0 0 0 0 0 0 0])[2]
3×20 Array{Float64,2}:
 2.0  1.73684  1.47368    1.21053   …  1.21053  1.47368  1.73684  2.0
 2.0  1.21053  0.421053  -0.368421     1.21053  1.47368  1.73684  2.0
 0.0  0.0      0.0        0.0          0.0      0.0      0.0      0.0

julia> dbspline(7,3,20,[0. 1 2 3 4 5 6; 1 2 1 2 1 2 3; 0 0 0 0 0 0 0])[3]
 -1.0  -1.0  -1.0  -1.0  0.0  0.0  0.0  …  0.0  0.0  0.0  1.0  1.0  1.0  1.0
 -3.0  -3.0  -3.0  -3.0  2.0  2.0  2.0     2.0  2.0  2.0  1.0  1.0  1.0  1.0
  0.0   0.0   0.0   0.0  0.0  0.0  0.0     0.0  0.0  0.0  0.0  0.0  0.0  0.0

julia> dbspline(7,3,20,[0. 1 2 3 4 5 6; 1 2 1 2 1 2 3; 0 0 0 0 0 0 0])[4]
19-element Array{Array{Int64,1},1}:
 [1, 2]  
 [2, 3]  
 [3, 4]  
 [4, 5]  
 [5, 6]  
 [6, 7]  
 [7, 8]  
 [8, 9]  
 [9, 10] 
 [10, 11]
 [11, 12]
 [12, 13]
 [13, 14]
 [14, 15]
 [15, 16]
 [16, 17]
 [17, 18]
 [18, 19]
 [19, 20]
```
---
_By Giuseppe Santorelli, Elia Onofri_
"""

function dbspline(npts::Int64,ord::Int64,p1::Int64,b::Array{Float64,2})::Tuple{Array{Float64},Array{Float64},Array{Float64},Array{Array{Int64,1},1}}
    
    @assert length(b) == 3 * npts ("ERROR: array b not matching with parameters")
    @assert npts >= ord ("ERROR: ord > npts")

    m::Int64 = npts + ord
    x = knot(npts,ord)
    
    G::Array{Float64} = transpose(b)
    
    N::Array{Float64} = zeros(p1, npts)
    N1::Array{Float64} = zeros(p1, npts)
    N2::Array{Float64} = zeros(p1, npts)
    
    t::Float64 = 0.0
    step::Float64 = x[m] / (p1 - 1)
    
    for i = 1 : p1
        if (x[m] - t) < 5e-6
            t = x[m]
        end
        
        nbasis, d1nbasis, d2nbasis = dbasis(npts,ord,t,x)
        
        for j = 1 : npts
            N[i,j] = nbasis[j]
            N1[i,j] = d1nbasis[j]
            N2[i,j] = d2nbasis[j]
        end
        t = t + step
    end
    
    P::Array{Float64} = transpose(N*G)
    D1::Array{Float64} = transpose(N1*G)
    D2::Array{Float64} = transpose(N2*G)
    
    EV = Array{Int64,1}[] #1-dimensional cellular complex used for plotting the curve with Plasm package

    for i = 1 : p1-1
        C = Array{Int64}(2)
        C[1] = i
        C[2] = i+1
        push!(EV,C) 
    end
    
    return (P, D1, D2, EV)
end


#--------------------------------------------------------------------------------------------------------------------------------

"""
    dbsplineu(npts,ord,p1,b)

Returns a tuple `(P,D1,D2,EV)`. P is the coordinates matrix of the p1 points b-spline curve and defined by the control points matrix b using an open knot vector. D1 and D2 are the first and second derivatives of P. EV is the 1-dimensional cellular complex. 

# Arguments

- `npts::Int64`: the number of the contol polygon vertices.
- `ord::Int64`: order of the bspline basis function.
- `p1::Int64`: number of the points to be calculated on the curve.
- `b::Array{Float64, 2}`: 2-dimensional Array containing the control polygon vertices in the x, y, z coordinates for columns

# Example

```jldoctest
julia> dbsplineu(4,3,20,[0. 3 6 9; 0 10 3 6; 0 0 0 0])[1]
3×20 Array{Float64,2}:
 1.5  1.81579  2.13158  2.44737  2.76316  …  6.55263  6.86842  7.18421  7.5
 5.0  5.95845  6.72853  7.31025  7.7036      4.05125  4.09003  4.23961  4.5
 0.0  0.0      0.0      0.0      0.0         0.0      0.0      0.0      0.0

julia> dbsplineu(4,3,20,[0. 3 6 9; 0 10 3 6; 0 0 0 0])[2]
3×20 Array{Float64,2}:
  3.0  3.0      3.0      3.0      …   3.0       3.0       3.0      3.0
 10.0  8.21053  6.42105  4.63158     -0.157895  0.894737  1.94737  3.0
  0.0  0.0      0.0      0.0          0.0       0.0       0.0      0.0

julia> dbsplineu(4,3,20,[0. 3 6 9; 0 10 3 6; 0 0 0 0])[3]
 0.0    0.0    0.0    0.0    0.0  …   0.0   0.0   0.0   0.0   0.0   0.0
 -17.0  -17.0  -17.0  -17.0  -17.0     10.0  10.0  10.0  10.0  10.0  10.0
   0.0    0.0    0.0    0.0    0.0      0.0   0.0   0.0   0.0   0.0   0.0

julia> dbspline(4,3,20,[0. 3 6 9; 0 10 3 6; 0 0 0 0])[4]
19-element Array{Array{Int64,1},1}:
 [1, 2]  
 [2, 3]  
 [3, 4]  
 [4, 5]  
 [5, 6]  
 [6, 7]  
 [7, 8]  
 [8, 9]  
 [9, 10] 
 [10, 11]
 [11, 12]
 [12, 13]
 [13, 14]
 [14, 15]
 [15, 16]
 [16, 17]
 [17, 18]
 [18, 19]
 [19, 20]
```
---
_By Giuseppe Santorelli, Elia Onofri_
"""

function dbsplineu(npts::Int64,ord::Int64,p1::Int64,b::Array{Float64,2})::Tuple{Array{Float64},Array{Float64},Array{Float64},Array{Array{Int64,1},1}}
    
    @assert length(b) == 3 * npts ("ERROR: array b not matching with parameters")
    @assert npts >= ord ("ERROR: ord > npts")

    m::Int64 = npts + ord
    x = knotu(npts,ord)
    
    G::Array{Float64} = transpose(b)
    
    N::Array{Float64} = zeros(p1, npts)
    N1::Array{Float64} = zeros(p1, npts)
    N2::Array{Float64} = zeros(p1, npts)
    
    t::Float64 = ord - 1
    step::Float64 = (npts - (ord - 1)) / (p1 - 1)
    
    for i = 1 : p1
        if (npts - t) < 5e-6
            t = npts
        end
        
        nbasis, d1nbasis, d2nbasis = dbasisu(npts,ord,t,x)
        
        for j = 1 : npts
            N[i,j] = nbasis[j]
            N1[i,j] = d1nbasis[j]
            N2[i,j] = d2nbasis[j]
        end
        t = t + step
        
    end
    
    P::Array{Float64} = transpose(N*G)
    D1::Array{Float64} = transpose(N1*G)
    D2::Array{Float64} = transpose(N2*G)
    
    EV = Array{Int64,1}[] #1-dimensional cellular complex used for plotting the curve with Plasm package

    for i = 1 : p1-1
        C = Array{Int64}(2)
        C[1] = i
        C[2] = i+1
        push!(EV,C) 
    end
    
    return (P, D1, D2, EV)
end


#--------------------------------------------------------------------------------------------------------------------------------

"""
matpbspl(npts,ord,p1,b)

Subroutine to generate a B-spline curve using matrix methods ad a periodic uniform knot vector. Returns a tuple `(P,EV)`. P is the coordinates matrix of the p1 points b-spline curve and defined by the control points matrix b. EV is the 1-dimensional cellular complex. 

# Arguments

- `npts::Int64`: the number of the contol polygon vertices.
- `ord::Int64`: order of the bspline basis function.
- `p1::Int64`: number of the points to be calculated on the curve
- `b::Array{Float64}`: 2-dimensional array containing the control polygon vertices

--------------------------------------------------------------------------------------------------------------------------------
_By Giuseppe Santorelli
"""

function matpbspl(npts::Int64, ord::Int64, p1::Int64, b::Array{Float64})::Tuple{Array{Float64,2},Array{Array{Int64,1},1}}
    
    @assert length(b) == 3 * npts ("ERROR: array b not matching with parameters")
    @assert npts >= ord ("ERROR: ord > npts")
    
    G::Array{Float64,2} = zeros(ord, 3)
    T::Array{Float64,2} = zeros(1, ord)
    N::Array{Float64,2}, fcoeff::Float64 = nmatrix(ord)

    P::Array{Float64} = Float64[]
    
    step::Float64 = 1 / (p1 - 1)
    
    for i = 0 : npts
        #building G matrix
        for j = 0 : ord - 1
            G[j + 1, 1] = b[1, mod(i + j, npts) + 1]
            G[j + 1, 2] = b[2, mod(i + j,npts) + 1]
            G[j + 1, 3] = b[3, mod(i + j, npts) + 1]
        end
               
        t::Float64 = 0.0
        
        #calculating T array for every parameter t
        while t <= 1 - step
            for j = 1 : ord
                T[1,j] = t^(ord - j)
            end
            
            T = fcoeff * T
            
            #matrix multiplication
            P1::Array{Float64} = T*N*G
            
            push!(P,P1[1,1])
            push!(P,P1[1,2])
            push!(P,P1[1,3])
            
            t = t + step
        end
    
    end
    
    P = reshape(P,3,:)

    EV = Array{Int64,1}[] #1-dimensional cellular complex used for plotting the curve with Plasm package

    for i = 1 : length(P)/3
        C = Array{Int64}(2)
        C[1] = i
        C[2] = i+1
        push!(EV,C) 
    end
    
    return P, EV
end


#--------------------------------------------------------------------------------------------------------------------------------

"""
nmatrix(ord)

Subroutine to calculate the general B-spline periodic basis matrix. Returns the matrix `N` and the coefficient `fcoeff`.

# Arguments

- `ord::Int64` : order of the periodic basis function

---------------------------------------------------------------------------------------------------------------------------------
_By Giuseppe Santorelli
"""

function nmatrix(ord::Int64)::Tuple{Array{Float64,2}, Float64}
    
    function Ni(n, i)
       return factorial(n) / (factorial(i) * factorial(n - i))
    end
    
    N::Array{Float64,2} = zeros(ord, ord)
    fcoeff::Float64 = 1 / factorial(ord - 1)
    
    for i = 0 : ord - 1
        t = Ni(ord - 1, i)
        
        for j = 0 : ord - 1
            s = 0
            
            for k = j : ord - 1
                s1 = (ord - (k + 1))^i
                s2 = (-1)^(k-j)
                s3 = Ni(ord, k - j)
                s = s + s1 * s2 * s3
            end
            
            N[i + 1, j + 1] = t*s
        end
    end
    
    return N, fcoeff
    
end

