export bezier, dbezier


"""
	bezier(npts, b[], dpts)

Calculate a Bezier curve.

---

# Arguments

- `npts::Int64`: the number of the contol polygon vertices.
- `b::Array{Float64, 2}`: 2-dimensional Array containing the control polygon vertices in the x, y, z coordinates for columns.
- `dpts::Int64`: number of data points to be calculated on the curve.

---

# Examples

```jldoctest
julia> bezier(4,[1.0 2 4 3; 1 3 3 1; 0 0 0 0], 7)
3×7 Array{Float64,2}:
 1.0  1.56481  2.18519  2.75  3.14815  3.26852  3.0
 1.0  1.83333  2.33333  2.5   2.33333  1.83333  1.0
 0.0  0.0      0.0      0.0   0.0      0.0      0.0

"""

function bezier(npts::Int64, b::Array{Float64,2}, dpts::Int64)::Array{Float64,2}

	@assert npts == length(b[1,:]) ("ERROR: there are not npts points of the polygon in b[]")

	step::Float64 = 1/(dpts-1)
	n::Int64 = npts-1
	
	C::Array{Float64,2} = zeros(npts, npts)
	for i = 0:n
		for j = 0:n-i
			if (n-j-i)%2==0
				C[i+1, j+1] = binomial(n-j, n-j-i)
			else
				C[i+1, j+1] = 0.0-binomial(n-j, n-j-i)
			end
		end
	end
	
	D::Array{Float64,2} = zeros(npts, npts)
	for i = 1:npts
		D[i,i] = abs(C[i,1])
	end
	
	G::Array{Float64,2} = transpose(b)
	
	
	T::Array{Float64,2} = ones(dpts, npts)
	
	t::Float64 = 0.0
	for i = 1:dpts
		for j = 0:n-1  #Per j=0 il valore è 1 ed è già impostato
			T[i, n-j]=T[i, npts-j]*t
		end
		t = t + step
	end
	return transpose(T*C*D*G)
end

#-----------------------------------------------------------------------

"""
	dbezier(npts, b, cpts)

Calculate a B~zier curve and its first and second derivatives.
"""

function dbezier(npts::Int64, b::Array{Float64}, cpts::Int64)::tuple{Array{Float64}, Array{Float64}, Array{Float64}}
	#[...]
	return ([0.0 0], [0.0 0], [0.0 0])
end

#-----------------------------------------------------------------------

"""
	bern_basis(n, i, t)

Calculate the Bernstein basis.
"""

function bern_basis(n::Int64, i::Int64, t::Float64)::Float64
	return binomial(n, i) * t^i * (1-t)^(n-i)
end
