export bspline, bsplineu, dbspline, dbsplineu, matpbspl, nmatrix


"""
	bsplfit(dpts, d, npts, ord)

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
	bspline(npts, k, p1, b)

Generate a B-spline curve using an open uniform knot vector.
"""

function bspline(npts::Int64, k::Int64, p1::Int64, b::Array{Float64})::Array{Float64}
	
end


#--------------------------------------------------------------------------------------------------------------------------------


"""
	bsplineu(npts k, p1, b)

Generate a B-spline curve using a periodic uniform knot vector.
"""

function bsplineu(npts::Int64, k::Int64, p1::Int64, b::Array{Float64})::Array{Float64}
	
end


#--------------------------------------------------------------------------------------------------------------------------------


"""
	dbspline(npts k, p1, b)

Generate a B-spline curve and its derivatives using an open uniform knot vector.
"""

function dbspline(npts::Int64, k::Int64, p1::Int64, b::Array{Float64})::tuple{Array{Float64}, Array{Float64}, Array{Float64}}
	
end


#--------------------------------------------------------------------------------------------------------------------------------


"""
	dbsplineu(npts k, p1, b)

Generate a B-spline curve and its derivatives using an open uniform knot vector .
"""

function dbsplineu(npts::Int64, k::Int64, p1::Int64, b::Array{Float64})::tuple{Array{Float64}, Array{Float64}, Array{Float64}}
	
end


#--------------------------------------------------------------------------------------------------------------------------------


"""
	matpbspl(ord, npts, p1, D)

Generate a B-spline curve using matrix methods and a periodic uniform knot vector.

---

# Arguments
- `ord::Int64`: order of the B-spline basis function.
- `npts::Int64`: number of control polygon vertices.
- `p1::Int64`: number of points to be calculated on the span of the curve.
- `D::Array{Float64,2}`: array containing the control polygon vertices.

---
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
_By Elia Onofri_
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


