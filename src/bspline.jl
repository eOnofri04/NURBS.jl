import knot, knotu, basis, dbasis, dbasisu
export bspline, bsplineu, dbspline, dbsplineu

"""
bspline(npts,ord,p1,b)

Returns a tuple `(P,EV)`. P is the coordinates matrix of the p1 points b-spline curve defined by the control points matrix b. EV is the 1-dimensional cellular complex.

# Arguments

- `npts::Int64`: the number of the contol polygon vertices.
- `ord::Int64`: order of the bspline basis function.
- `p1::Int64`: number of the points to be calculated on the curve
- `b::Array{Float64}`: 2-dimensional Array containing the control polygon vertices

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

"""


function bspline(npts::Int64, ord::Int64, p1::Int64, b::Array{Float64})::Tuple{Array{Float64,2},Array{Array{Int64,1},1}}
    
    @assert length(b) == 3 * npts ("ERROR: array b not matching with parameters")
    @assert npts >= ord ("ERROR: ord > npts")

    m::Int64 = npts + ord
    P::Array{Float64} = zeros(3,p1)
    
    x::Array{Float64} = knot(npts,ord) #creates knot array
    
    icount::Int64 = 1
    t::Float64 = 0.0
    
    step::Float64 = x[m] / (p1-1)
    
    for i1 = 1 : p1 #computes the value of the B-spline in the p1 points
        if (x[m] - t) < 5e-6
            t = x[m]
        end
        
        nbasis::Array{Float64} = basis(npts,ord,t,x) #computes the basis value depending on patameter t
        
        for j = 1 : 3 #iteration over three components for every point
            for i = 1 : npts
                temp = nbasis[i] * b[j,i]
                P[j,icount] = P[j,icount] + temp
            end
        end
        
        icount = icount + 1
        t = t + step
    end
    
    EV = Array{Int64,1}[] #1-dimensional cellular complex used for plotting the curve with Plasm package

    for i = 1 : p1-1
        C = Array{Int64}(2)
        C[1] = i
        C[2] = i+1
        push!(EV,C) 
    end
    
    return P,EV
end


#-------------------------------------------------------------------------------------


"""
bsplineu(npts,ord,p1,b)

Returns a tuple `(P,EV)`. P is the coordinates matrix of the p1 points b-spline curve defined by the control points matrix b using an open knot vector. EV is the 1-dimensional cellular complex.

# Arguments

- `npts::Int64`: the number of the contol polygon vertices.
- `ord::Int64`: order of the bspline basis function.
- `p1::Int64`: number of the points to be calculated on the curve
- `b::Array{Float64}`: 2-dimensional array containing the control polygon vertices

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

"""

function bsplineu(npts::Int64,ord::Int64,p1::Int64,b::Array{Float64})::Tuple{Array{Float64,2},Array{Array{Int64,1},1}}
    
    @assert length(b) == 3 * npts ("ERROR: array b not matching with parameters")
    @assert npts >= ord ("ERROR: ord > npts")

    m::Int64 = npts + ord
    P::Array{Float64} = zeros(3,p1)
        
    x::Array{Float64} = knotu(npts,ord) #creates knot array
    
    icount::Int64 = 1
    t::Float64 = ord - 1
    step::Float64 = (npts - (ord - 1)) / (p1 - 1)
    
    for i1 = 1 : p1 #computes the value of the B-spline in the p1 points
        if (npts - t) < 5e-6
            t = npts
        end
        
        nbasis::Array{Float64} = basis(npts,ord,t,x) #computes the basis value depending on patameter t
        
        for j = 1 : 3 #iteration over three components for every point
            for i = 1 : npts
                temp = nbasis[i] * b[j,i]
                P[j,icount] = P[j,icount] + temp
            end
        end
        
        icount = icount + 1
        t = t + step
    end
    
    EV = Array{Int64,1}[] #1-dimensional cellular complex used for plotting the curve with Plasm package

    for i = 1 : p1-1
        C = Array{Int64}(2)
        C[1] = i
        C[2] = i+1
        push!(EV,C) 
    end
        
    return P,EV
    
end


#------------------------------------------------------------------------------------


"""
dbspline(npts,ord,p1,b)

Returns a tuple `(P,D1,D2,EV)`. P is the coordinates matrix of the p1 points b-spline curve and defined by the control points matrix b. D1 and D2 are the first and second derivatives of P. EV is the 1-dimensional cellular complex. 

# Arguments

- `npts::Int64`: the number of the contol polygon vertices.
- `ord::Int64`: order of the bspline basis function.
- `p1::Int64`: number of the points to be calculated on the curve
- `b::Array{Float64}`: 2-dimensional array containing the control polygon vertices

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
"""

function dbspline(npts::Int64,ord::Int64,p1::Int64,b::Array{Float64})::Tuple{Array{Float64},Array{Float64},Array{Float64},Array{Array{Int64,1},1}}
    
    @assert length(b) == 3 * npts ("ERROR: array b not matching with parameters")
    @assert npts >= ord ("ERROR: ord > npts")

    m::Int64 = npts + ord

    P::Array{Float64} = zeros(3,p1)
    D1::Array{Float64} = zeros(3,p1)
    D2::Array{Float64} = zeros(3,p1)

    x = knot(npts,ord)
    
    icount::Int64 = 1
    t::Float64 = 0.0
    step::Float64 = x[m] / (p1 - 1)
    
    for i1 = 1 : p1
        if (x[m] - t) < 5e-6
            t = x[m]
        end
        
        nbasis, d1nbasis, d2nbasis = dbasis(npts,ord,t,x)
        
        for j = 1 : 3       
            for i = 1 : npts
                P[j,icount] = P[j,icount] + (nbasis[i] * b[j,i])
                D1[j,icount] = D1[j,icount] + (d1nbasis[i] * b[j,i])
                D2[j,icount] = D2[j,icount] + (d2nbasis[i] * b[j,i])
            end
        end
        
        icount = icount + 1
        t = t + step
    end
    
     EV = Array{Int64,1}[] #1-dimensional cellular complex used for plotting the curve with Plasm package

    for i = 1 : p1-1
        C = Array{Int64}(2)
        C[1] = i
        C[2] = i+1
        push!(EV,C) 
    end
    
    return P, D1, D2, EV
    
end


#------------------------------------------------------------------------------------

"""
dbsplineu(npts,ord,p1,b)

Returns a tuple `(P,D1,D2,EV)`. P is the coordinates matrix of the p1 points b-spline curve and defined by the control points matrix b using an open knot vector. D1 and D2 are the first and second derivatives of P. EV is the 1-dimensional cellular complex. 

# Arguments

- `npts::Int64`: the number of the contol polygon vertices.
- `ord::Int64`: order of the bspline basis function.
- `p1::Int64`: number of the points to be calculated on the curve
- `b::Array{Float64}`: 2-dimensional array containing the control polygon vertices

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

"""

function dbsplineu(npts::Int64,ord::Int64,p1::Int64,b::Array{Float64})::Tuple{Array{Float64},Array{Float64},Array{Float64},Array{Array{Int64,1},1}}
    
    @assert length(b) == 3 * npts ("ERROR: array b not matching with parameters")
    @assert npts >= ord ("ERROR: ord > npts")

    m::Int64 = npts + ord

    P::Array{Float64} = zeros(3,p1)
    D1::Array{Float64} = zeros(3,p1)
    D2::Array{Float64} = zeros(3,p1)
    
    x = knotu(npts,ord)
    
    icount::Int64 = 1
    t::Float64 = ord - 1
    step::Float64 = (npts - (ord - 1)) / (p1 - 1)
    
    for i1 = 1 : p1
        if (npts - t) < 5e-6
            t = npts
        end
        
        nbasis, d1nbasis, d2nbasis = dbasisu(npts,ord,t,x)
        
        for j = 1 : 3
            for i = 1 : npts
                P[j,icount] = P[j,icount] + (nbasis[i] * b[j,i])
                D1[j,icount] = D1[j,icount] + (d1nbasis[i] * b[j,i])
                D2[j,icount] = D2[j,icount] + (d2nbasis[i] * b[j,i])
            end
        end
        
        icount = icount + 1
        t = t + step
    end
    
    EV = Array{Int64,1}[] #1-dimensional cellular complex used for plotting the curve with Plasm package

    for i = 1 : p1-1
        C = Array{Int64}(2)
        C[1] = i
        C[2] = i+1
        push!(EV,C) 
    end
    
    return P, D1, D2, EV
    
end