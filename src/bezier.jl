export bezier, dbezier, bern_basis


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

- `npts::Int64`: the number of the control polygon vertices.
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

	@assert npts == length(b[1,:]) ("ERROR: there are not npts = $(npts) points of the polygon in b[]")

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
		D[i,i] = abs(C[npts+1-i,1])
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
	dbezier(npts, b[], dpts)

Calculate Bezier curve and its first and second derivatives over the `npts` control polygon points in `b`.


The method evaluates `dpts` points of the curve and gives back the tuple `(P, D1, D2, EV, DEV)`.

Here `P` is the vector of points memorized by columns while `EV` is the `1`-dimensional cellular complex associated to it.
Moreover `D1` and `D2` are the derivatives vectors and `DEV` is the `1`-dimensional cellular complex associated.
Please note that only `dpts-2` points are evaluated for the derivatives as there are not enought informations about the first and the last points.

This structure as it is formed could be directly imported in Plasm module by typing:
```julia
julia> P, _, _, EV, _ = bezier(npts, b[], dpts);
julia> Plasm.view(P, EV)
```

---

# Arguments

- `npts::Int64`: the number of the control polygon vertices.
- `b::Array{Float64, 2}`: 2-dimensional Array containing the control polygon vertices in the x, y, z coordinates for columns.
- `dpts::Int64`: number of data points to be calculated on the curve.

---

# Examples

```jldoctest
julia> dbezier(4,[1.0 2 4 3; 1 3 3 1; 0 0 0 0], 7)[1]
3×7 Array{Float64,2}:
 1.0  1.56481  2.18519  2.75  3.14815  3.26852  3.0
 1.0  1.83333  2.33333  2.5   2.33333  1.83333  1.0
 0.0  0.0      0.0      0.0   0.0      0.0      0.0
```

```jldoctest
julia> dbezier(4,[1.0 2 4 3; 1 3 3 1; 0 0 0 0], 7)[2]
3×5 Array{Float64,2}:
 3.66667  3.66667  3.0   1.66667  3.26852
 4.0      2.0      0.0  -2.0      1.83333
 0.0      0.0      0.0   0.0      0.0
```

```jldoctest
julia> dbezier(4,[1.0 2 4 3; 1 3 3 1; 0 0 0 0], 7)[3]
3×5 Array{Float64,2}:
   2.0   -2.0   -6.0  -10.0  3.26852
 -12.0  -12.0  -12.0  -12.0  1.83333
   0.0    0.0    0.0    0.0  0.0
```

```jldoctest
julia> dbezier(4,[1.0 2 4 3; 1 3 3 1; 0 0 0 0], 7)[4]
6-element Array{Array{Int64,1},1}:
 [1, 2]
 [2, 3]
 [3, 4]
 [4, 5]
 [5, 6]
 [6, 7]
```

```jldoctest
julia> dbezier(4,[1.0 2 4 3; 1 3 3 1; 0 0 0 0], 7)[5]
6-element Array{Array{Int64,1},1}:
 [1, 2]
 [2, 3]
 [3, 4]
 [4, 5]
```

---

_By Elia Onofri_
"""

function dbezier(npts::Int64, b::Array{Float64,2}, dpts::Int64)::Tuple{Array{Float64,2}, Array{Float64,2}, Array{Float64,2}, Array{Array{Int64,1},1}, Array{Array{Int64,1},1}}

	@assert npts == length(b[1,:]) ("ERROR: there are not npts = $(npts) points of the polygon in b[]")

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
		D[i,i] = abs(C[npts+1-i,1])
	end

	G::Array{Float64,2} = transpose(b)

	T::Array{Float64,2} = ones(dpts, npts)

	t::Float64 = 0.0
	for i = 1:dpts
		for j = 0:n-1  #Per j=0 il valore è 1 ed è già impostato
			T[i, n-j] = T[i, npts-j]*t
		end
		t = t + step
	end

	Der1::Array{Float64,2} = ones(dpts-2, npts)
	Der2::Array{Float64,2} = ones(dpts-2, npts)
    
	t = step
	for i = 2:dpts-2
		den1::Float64 = t * (1 - t)
		den2::Float64 = den1^2
		for j = 0:n
			Der1[i-1, j+1] = (j - n * t) / den1
			Der2[i-1, j+1] = ((j - n * t)^2 - n * t^2 - j * (1 - 2 * t)) / den2
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

	DEV::Array{Array{Int64,1},1} = EV[1:dpts-3] #1-dimensional cellular complex used for plotting the Derivatives with Plasm package

	F::Array{Float64,2} = T*C*D
	F1::Array{Float64,2} = Der1.*(F[range(2,dpts-2),:])
	F2::Array{Float64,2} = Der2.*(F[range(2,dpts-2),:])
	P::Array{Float64,2} = transpose(F*G)
	D1::Array{Float64,2} = transpose(F1*G)
	D2::Array{Float64,2} = transpose(F2*G)

	return (P, D1, D2, EV, DEV)
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
