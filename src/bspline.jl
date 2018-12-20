import knot, knotu, basis, dbasis, dbasisu
export bspline, bsplineu, dbspline, dbsplineu, matpbspl, nmatrix

"""
	bsplfit(dpts, d[], npts, ord)

Fit a B-spline curve using an open uniform knot vector.

---

# Arguments
- `dpts::Int64`: Number of Data Points.
- `d::Array{Float64}`: Array of Data Points.
- `npts::Int64`: Number of the Control Poligon Vertices.
- `ord::Int64`: the order of the B-Spline (`degree k-1`).

---

# Examples
```jldoctest
julia> 

```
```jldoctest
julia> 

```
```jldoctest
julia> 

```
---
_By Gianmarco Caramitti, Elia Onofri_
"""

function bsplfit(dpts::Int64, d::Array{Float64}, npts::Int64, ord::Int64)::Array{Float64}

	b::Array{Float64} = zeros(dpts, 2) 		# Array containing the control polygon vertices
	nbasis::Array{Float64} = zeros(npts) 		# Array containing the basis ]unctions ]or a single value of t
	n::Array{Float64} = zeros(dpts, npts) 		# Matrix of basis function
	ntrn::Array{Float64} = zeros(npts, dpts) 	# Transpose of the n matrix
	ntemp::Array{Float64} = zeros(npts, npts) 	# Temporary matrix to hold trn(n) x n
	ntmp::Array{Float64} = zeros(npts, 2) 		# Temporary matrix to hold inverse of trn(n) x n x d
	ninv::Array{Float64} = zeros(npts, npts) 	# Inverse of trn(n) x n
	m::Int8 = npts + ord						# Number of knot values

	x::Array{Float64} = knot(npts, ord)
	tpar::Array{Float64} = param(dpts, d)
	

	# Generate the matrix of basis functions
	for i = 1 : dpts

		t = tpar[i] * x[m]
		nbasis = basis(ord, t, npts, x)

		for j = 1 : npts

			n[i,j] = nbasis[j]

		end

	end

	ntrn = n'
	ntemp = ntrn * n
	ninv = inv(ntemp)
	ntmp = ntrn * d
	b = ninv * ntmp

	return b

end


#--------------------------------------------------------------------------------------------------------------------------------

"""
	bspline(npts, ord, p1, b[])

Generate a B-spline curve using an open uniform knot vector.

---

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
---
_By Giuseppe Santorelli, Elia Onofri_
"""

function bspline(npts::Int64, ord::Int64, p1::Int64, b::Array{Float64})::Tuple{Array{Float64,2},Array{Array{Int64,1},1}}
    
    @assert length(b) == 3 * npts ("ERROR: array b not matching with parameters")
    @assert npts >= ord ("ERROR: ord > npts")

    m::Int64 = npts + ord
    P::Array{Float64} = zeros(3, p1)
    
    x::Array{Float64} = knot(npts, ord) #creates knot array
    
    icount::Int64 = 1
    t::Float64 = 0.0
    
    step::Float64 = x[m] / (p1 - 1)
    
    for i1 = 1 : p1 #computes the value of the B-spline in the p1 points
        if (x[m] - t) < 5e-6
            t = x[m]
        end
        
        nbasis::Array{Float64} = basis(npts,ord,t,x) #computes the basis value depending on patameter t
        
        for i = 1 : npts
            for j = 1 : 3 #iteration over three components for every point
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
- `b::Array{Float64}`: 2-dimensional array containing the control polygon vertices.

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
---
_By Giuseppe Santorelli, Elia Onofri_
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
    
    return (P, D1, D2, EV)
end


#--------------------------------------------------------------------------------------------------------------------------------

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
---
_By Giuseppe Santorelli, Elia Onofri_
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
    
    return (P, D1, D2, EV)
end


#--------------------------------------------------------------------------------------------------------------------------------


"""
	matpbspl(ord, npts, p1, D[])

Generate a B-spline curve using matrix methods and a periodic uniform knot vector.

---

# Arguments
- `ord::Int64`: order of the B-spline basis function.
- `npts::Int64`: number of control polygon vertices.
- `p1::Int64`: number of points to be calculated on the span of the curve.
- `D::Array{Float64,2}`: array containing the control polygon vertices.

---

# Examples
```jldoctest
julia> 

```
```jldoctest
julia> 

```
```jldoctest
julia> 

```
---
_By Paolo Macciacchera, Elia Onofri_
"""

function matpbspl(ord::Int64, npts::Int64, p1::Int64, D::Array{Float64,2})::Array{Float64,2}
    b::Array{Float64,2} = zeros(3, ord)
    n::Array{Float64,2} = zeros(ord, ord)
    t::Array{Float64,2} = zeros(1, ord)
    temp::Array{Float64,2} = zeros(1, ord)
    p::Array{Float64,2} = zeros(3, p1 * npts)
    fcoeff::Float64 = 0.0
    
    n, fcoeff = nmatrix(ord)
    icount = 0
    stepl = 1 / p1
    rangel = 0 : stepl : 1 - (1 / p1)
    
    for j = 1 : npts
        for k = 0 : ord-1
            b[1, k + 1] = D[1, mod(j + k, npts) + 1]
            b[2, k + 1] = D[2, mod(j + k, npts) + 1]
            b[3, k + 1] = D[3, mod(j + k, npts) + 1]
        end
        for l in rangel         
            icount = icount + 1               
            for i = 1 : ord
                t[i] = l ^ (ord - i)
            end                      
            t = fcoeff * t
            temp = t * n
            ptemp::Array{Float64} = temp * b'
            p[1, icount] = ptemp[1]
            p[2, icount] = ptemp[2]
            p[3, icount] = ptemp[3]
        end
    end
    return p            
end


#--------------------------------------------------------------------------------------------------------------------------------


"""
	nmatrix(k)

Calculate the general B-spline periodic basis matrix.

This function is used in matpbspl.

---

# Arguments
- `ord::Int64`: order of the periodic basis function.

---

# Examples
```jldoctest
julia> 

```
```jldoctest
julia> 

```
```jldoctest
julia> 

```
---
_By Paolo Macciacchera, Elia Onofri_
"""

function nmatrix(ord::Int64)::Tuple{Array{Float64,2}, Float64}
    n::Array{Float64} = Array{Float64}(zeros(ord, ord))
    fcoeff = 1 / factorial(ord - 1)
    Ni = (n,i) -> factorial(n) / (factorial(i) * factorial(n - i)) 
 
    for i = 0 : ord - 1
        temp = Ni(ord - 1, i)
        for j = 0 : ord - 1
            sum = 0
            for k = j : ord - 1
                sum1 = ( ord - ( k + 1 ) ) ^ i
                sum2 = (-1) ^ (k - j)
                sum3 = Ni(ord, k - j)
                sum = sum + sum1 * sum2 * sum3
            end
            n[i + 1, j + 1] = temp * sum
        end
    end
    return(n, fcoeff)            
end


#--------------------------------------------------------------------------------------------------------------------------------


"""
	param(dpts, d[])

Calculate parameter values based on chord distances.

---

# Arguments
- `dpts::Int64`: Number of Data Points.
- `d::Array{Float64}`: Array of Data Points.

---

# Examples
```jldoctest
julia> 

```
```jldoctest
julia> 

```
```jldoctest
julia> 

```
---
_By Gianmarco Caramitti, Elia Onofri_
"""

function param(dpts::Int64, d::Array{Float64})::Array{Float64}

	tparm::Array{Float64} = zeros(dpts)

    # calculate the chord distances for all the data points
    arc_lengths::Array{Float64} = [sqrt((d[i,1] - d[i-1,1])^2 + (d[i,2] - d[i-1,2])^2) for i = 2 : dpts]

    tparm[1] = 0

    for i = 2 : dpts
        tparm[i] = tparm[i-1] + arc_lengths[i-1] 
    end

    return tparm / sum(arc_lengths)

end

