export rbasis, rbspline, rbsplinu


"""
	rbasis(npts,ord,t,x,h)

Subroutine to generate a rational B-spline basis using an open knot vector.

# Arguments
- `npts::Int64` : the number of the points of the control polygon.
- `ord::Int64` : the order of the B-Spline (`deg = ord-1`).
- `t::Float64` : the parameter value of the parametric curve.
- `x::Array{Float64}` : the knot vector.
- `h::Array{Float64}` : array containing the homogeneous coordinate weighting factors

# Examples

```jldoctest
julia> x = [0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 2.0]

julia> h = [1.0, 1.0, 0.25, 1.0, 1.0]

julia> rbasis(5, 4, 0.5, x, h)
5-element Array{Float64,1}:
 0.153846 
 0.730769 
 0.0769231
 0.0384615
 0.0        
```
_By Giuseppe Santorelli_
"""

function rbasis(npts::Int64, ord::Int64, t::Float64, x::Array{Float64}, h::Array{Float64})::Array{Float64}

	m::Int64 = npts + ord
	r::Array{Float64} = zeros(m)
    
	# first order nonrational basis functions
	for i = 1 : m-1
		if (t >= x[i]) && (t < x[i+1])
			r[i] = 1
		else
			r[i] = 0
		end
	end
    
	# higher order nonrational basis function
	d::Float64 = 0.0
	e::Float64 = 0.0

	for k = 2 : ord

		for i = 1 : m-ord

			if r[i] != 0
				d = ((t - x[i]) * r[i]) / (x[i+k-1] - x[i])
			else
				d = 0
			end
            
			if r[i+1] != 0
				e = ((x[i+k] - t) * r[i+1]) / (x[i+k] - x[i+1])
			else
				e = 0
			end

			r[i] = d + e
		end
	end
    
	# pick up last point
	if t == x[m]
		r[npts] = 1
	end
    
    # calculate sum for denominator of rational basis functions
	sum::Float64 = 0.0

	for i = 1 : npts
		sum = sum + r[i] * h[i]
	end
    
	# form rational basis function
	if sum != 0
		for i = 1 : npts
			r[i] = r[i] * h[i] / sum
		end
	else
		r = zeros(npts)
	end
        
	return r[1:npts]
end

#-----------------------------------------------------------------------

"""
	rbspline(npts, ord, p1, b, h)

Subroutine to generate a rational B-spline curve using an uniform open knot vector. Returns a tuple `(P,EV)`. P is the coordinates matrix of the p1 points b-spline curve defined by the control points matrix b. EV is the 1-dimensional cellular complex.


# Arguments

- `npts::Int64` : the number of the points of the control polygon.
- `ord::Int64` : the order of the B-Spline (`deg = ord-1`).
- `p1::Int64`: number of the points to be calculated on the curve
- `b::Array{Float64}`: 2-dimensional array containing the control polygon vertices
- `h::Array{Float64}` : array containing the homogeneous coordinate weighting factors

# Examples

```jldoctest
julia> b = [0 1 5/2 4 5; 1 2 0 2 0; 0 0 0 0 0]

julia> h = [1, 1, 0.25, 1, 1]

julia> rbspline(5,3,50,b,h)[1]
3×50 Array{Float64,2}:
 0.0  0.118164  0.228377  0.331479  …  4.66852   4.77162  4.88184   5.0
 1.0  1.11652   1.22178   1.31653      0.641603  0.44733  0.233982  0.0
 0.0  0.0       0.0       0.0          0.0       0.0      0.0       0.0

julia> rbspline(5,3,50,b,h)[2]
49-element Array{Array{Int64,1},1}:
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
 ⋮       
 [38, 39]
 [39, 40]
 [40, 41]
 [41, 42]
 [42, 43]
 [43, 44]
 [44, 45]
 [45, 46]
 [46, 47]
 [47, 48]
 [48, 49]
 [49, 50]
```

_By Giuseppe Santorelli_
"""
function rbspline(npts::Int64, ord::Int64, p1::Int64, b::Array{Float64}, h::Array{Float64})::Tuple{Array{Float64,2},Array{Array{Int64,1},1}}
    
	@assert length(b) == 3 * npts ("ERROR: array b not matching with parameters")
	@assert npts >= ord ("ERROR: ord > npts")
    
	m::Int64 = npts + ord
	P::Array{Float64} = zeros(3, p1)
    
	#generate the uniform knot vector
	x::Array{Float64} = knot(npts, ord)
    
	#calculate the points on the rational B-spline curve
	icount::Int64 = 1
	t::Float64 = 0.0
	step::Float64 = x[m] / (p1 - 1)
    
	for i1 = 1 : p1
       
		if (x[m] - t < 5e-6)
			t = x[m]
		end
        
		#generate the basis function
		r::Array{Float64} = rbasis(npts, ord, t, x, h)
        
		#generate a point
		for j = 1 : 3 #iteration over three components for every point
			for i = 1 : npts
				temp = r[i] * b[j,i]
				P[j, icount] = P[j, icount] + temp
			end
		end
        
		icount = icount + 1
		t = t + step
	end
    
	EV = Array{Int64,1}[] #1-dimensional cellular complex used for plotting the curve with Plasm package

	for i = 1 : p1 - 1
		C = Array{Int64}(2)
		C[1] = i
		C[2] = i + 1
		push!(EV, C) 
	end
    
	return (P, EV)
end

#-----------------------------------------------------------------------

"""
	rbsplinu(npts,ord,p1,b,h)

Subroutine to generate a rational B-spline curve using an uniform periodic knot vector. Returns a tuple `(P,EV)`. P is the coordinates matrix of the p1 points b-spline curve defined by the control points matrix b. EV is the 1-dimensional cellular complex.


# Arguments

- `npts::Int64` : the number of the points of the control polygon.
- `ord::Int64` : the order of the B-Spline (`deg = ord-1`).
- `p1::Int64`: number of the points to be calculated on the curve
- `b::Array{Float64}`: 2-dimensional array containing the control polygon vertices
- `h::Array{Float64}` : array containing the homogeneous coordinate weighting factors

# Examples

```jldoctest
julia> b = [0 1 5/2 4 5; 1 2 0 2 0; 0 0 0 0 0]

julia> h = [1, 1, 0.25, 1, 1]

julia> rbspline(5,3,50,b,h)[1]
3×50 Array{Float64,2}:
 0.5  0.559434  0.615602  0.668943  …  4.33106  4.3844   4.44057  4.75
 1.5  1.55779   1.60901   1.65399      1.31653  1.22178  1.11652  0.5 
 0.0  0.0       0.0       0.0          0.0      0.0      0.0      0.0 

julia> rbspline(5,3,50,b,h)[2]
49-element Array{Array{Int64,1},1}:
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
 ⋮       
 [38, 39]
 [39, 40]
 [40, 41]
 [41, 42]
 [42, 43]
 [43, 44]
 [44, 45]
 [45, 46]
 [46, 47]
 [47, 48]
 [48, 49]
 [49, 50]
```
_By Giuseppe Santorelli_
"""
function rbsplinu(npts::Int64, ord::Int64, p1::Int64, b::Array{Float64}, h::Array{Float64})::Tuple{Array{Float64,2},Array{Array{Int64,1},1}}
    
	@assert length(b) == 3 * npts ("ERROR: array b not matching with parameters")
	@assert npts >= ord ("ERROR: ord > npts")
    
	m::Int64 = npts + ord
	P::Array{Float64} = zeros(3, p1)
    
	#generate the uniform periodic knot vector
	x::Array{Float64} = knotu(npts, ord)
    
	#calculate the points on the rational B-spline
	icount::Int64 = 1
	t::Float64 = ord - 1
	step::Float64 = (npts - ord + 1) / (p1 - 1)
    
	for i1 = 1 : p1
       
	if (x[m] - t < 5e-6)
			t = x[m]
		end
        
		#generate the basis function
		r::Array{Float64} = rbasis(npts, ord, t, x, h)
        
		#generate a point
		for j = 1 : 3 #iteration over three components for every point
			for i = 1 : npts
				temp = r[i] * b[j,i]
				P[j,icount] = P[j,icount] + temp
			end
		end
        
		icount = icount + 1
		t = t + step
	end
    
	EV = Array{Int64,1}[] #1-dimensional cellular complex used for plotting the curve with Plasm package

	for i = 1 : p1 - 1
		C = Array{Int64}(2)
		C[1] = i
		C[2] = i + 1
		push!(EV, C)
	end
    
	return (P, EV)
end
