export bezier, dbezier


"""
	bezier(npts, b[], dpts)

Calculate a Bezier curve over the `npts` control polygon points in `b`.

The method evaluates `dpts` points of the curve and gives back the tuple `(P, EV)`.
Here `P` is the vector of points memorized by columns while `EV` is the `1`-dimensional cellular complex.
This structure as it is formed could be directly imported in Plasm module by typing:
```julia
julia> P, EV = bezier(npts, b[], dpts);
julia> Plasm.view(P, EV)
```

---

# Arguments

- `npts::Int64`: the number of the contol polygon vertices.
- `b::Array{Float64, 2}`: 2-dimensional Array containing the control polygon vertices in the x, y, z coordinates for columns.
- `dpts::Int64`: number of data points to be calculated on the curve.

---

# Examples

```jldoctest
julia> bezier(4,[1.0 2 4 3; 1 3 3 1; 0 0 0 0], 7)[1]
3×7 Array{Float64,2}:
 1.0  1.56481  2.18519  2.75  3.14815  3.26852  3.0
 1.0  1.83333  2.33333  2.5   2.33333  1.83333  1.0
 0.0  0.0      0.0      0.0   0.0      0.0      0.0
```

```jldoctest
julia> bezier(4,[1.0 2 4 3; 1 3 3 1; 0 0 0 0], 7)[2]
6-element Array{Array{Int64,1},1}:
 [1, 2]
 [2, 3]
 [3, 4]
 [4, 5]
 [5, 6]
 [6, 7]
```

---

_By Elia Onofri_
"""

function bezier(npts::Int64, b::Array{Float64,2}, dpts::Int64)::Tuple{Array{Float64,2}, Array{Array{Int64,1},1}}

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

	EV::Array{Array{Int64,1},1} = Array{Int64,1}[] #1-dimensional cellular complex used for plotting the curve with Plasm package

    for i = 1 : dpts-1
        app = Array{Int64}(2)
        app[1] = i
        app[2] = i+1
        push!(EV,app) 
    end

    P::Array{Float64,2} = transpose(T*C*D*G)

	return (P, EV)
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

This method is meant to evaluate the `i`-th Bernstein Basis of order `n` in the parameter `t`.
Of course this method could be implemented in a more quick way, in order to evaluate all the basis of a specific order
but this is out of the scope of this library.

---

# Arguments

- `n::Int64`: the order of the basis (`n = npts-1`).
- `i::Int64`: the `i`-th basis of the given order. `i ∈ [0, …, n]`.
- `t::Float64`: the parameter of the basis `N_{n, i}(t)`.

---

# Examples

```jldoctest
julia> map(x -> bern_basis(3, 0, x), [0, 0.15, 0.35, 0.5, 0.65, 0.85, 1])
7-element Array{Float64,1}:
 1.0     
 0.614125
 0.274625
 0.125   
 0.042875
 0.003375
 0.0   
```

```jldoctest
julia> map(x -> bern_basis(3, 3, x), [0, 0.15, 0.35, 0.5, 0.65, 0.85, 1])
7-element Array{Float64,1}:
 0.0     
 0.003375
 0.042875
 0.125   
 0.274625
 0.614125
 1.0 
```

---

_By Elia Onofri_
"""

function bern_basis(n::Int64, i::Int64, t::Float64)::Float64
	@assert i >= 0 && i<= n ("ERROR: there are $(n+1) basis of order $(n). You cannot choose the $(i)-th one.")
	return binomial(n, i) * t^i * (1-t)^(n-i)
end
