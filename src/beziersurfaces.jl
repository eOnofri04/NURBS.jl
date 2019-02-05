export bezsurf

"""
	bezsurf(npts, mptsg, b[], udpts, wdpts)

Calculate a Bezier surface over the `npts` by `mpts` control net points in `b`.

The method evaluates `udpts` by `wdpts` points of the curve and gives back the tuple `(P, EV, FV)`.
Here:
 - `P` is the vector of points memorized by columns;
 - `EV` is the `1`-dimensional cellular complex of the edges;
 - `FV` is the `2`-dimensional cellular complex of the faces.

This structure as it is formed could be directly imported in Plasm module by typing:
```julia
julia> P, EV, FV = bezier(npts, b[], dpts);
julia> Plasm.view(Plasm.mkpol(P, FV))
```

---

# Arguments

- `npts::Int64`: the number of the control net rows in `u` direction.
- `mpts::Int64`: the number of the control net rows in `w` direction;
- `b::Array{Float64, 2}`: 2-dimensional Array containing the control polygon vertices in the x, y, z coordinates for columns.
- `udpts::Int64`: number of data net rows in `u` direction.
- `wdpts::Int64`: number of data net rows in `w` direction.

---

# Examples

```jldoctest
julia> b = [
		  0.0 0.0 0.0 1.0 1.0 1.0 2.0 2.0 2.0;
		  0.0 2.0 4.0 0.0 2.0 4.0 0.0 2.0 4.0;
		  2.0 0.0 2.0 0.0 1.0 0.0 2.0 0.0 2.0
		]
3×9 Array{Float64,2}:
 0.0  0.0  0.0  1.0  1.0  1.0  2.0  2.0  2.0
 0.0  2.0  4.0  0.0  2.0  4.0  0.0  2.0  4.0
 2.0  0.0  2.0  0.0  1.0  0.0  2.0  0.0  2.0

julia> bezsurf(3, 3, b, 5, 5)[1]
3×25 Array{Float64,2}:
 0.0  0.5   1.0  1.5   2.0  0.0   0.5       1.0     1.5       2.0   0.0  0.5     1.0   1.5     2.0  0.0   0.5       1.0     1.5       2.0   0.0  0.5   1.0  1.5   2.0
 0.0  0.0   0.0  0.0   0.0  1.0   1.0       1.0     1.0       1.0   2.0  2.0     2.0   2.0     2.0  3.0   3.0       3.0     3.0       3.0   4.0  4.0   4.0  4.0   4.0
 2.0  1.25  1.0  1.25  2.0  1.25  0.921875  0.8125  0.921875  1.25  1.0  0.8125  0.75  0.8125  1.0  1.25  0.921875  0.8125  0.921875  1.25  2.0  1.25  1.0  1.25  2.0

julia> bezsurf(3, 3, b, 5, 5)[2]
44-element Array{Array{Int64,1},1}:
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
 ⋮       
 [12, 17]
 [13, 18]
 [14, 19]
 [15, 20]
 [16, 21]
 [17, 22]
 [18, 23]
 [19, 24]
 [20, 25]

julia> bezsurf(3, 3, b, 5, 5)[3]
16-element Array{Array{Int64,1},1}:
 [1, 2, 6, 7]    
 [2, 3, 7, 8]    
 [3, 4, 8, 9]    
 [4, 5, 9, 10]   
 [6, 7, 11, 12]  
 [7, 8, 12, 13]  
 [8, 9, 13, 14]  
 [9, 10, 14, 15] 
 [11, 12, 16, 17]
 [12, 13, 17, 18]
 [13, 14, 18, 19]
 [14, 15, 19, 20]
 [16, 17, 21, 22]
 [17, 18, 22, 23]
 [18, 19, 23, 24]
 [19, 20, 24, 25]
```

---

_By Elia Onofri_
"""

function bezsurf(npts::Int64, mpts::Int64, b::Array{Float64,2}, udpts::Int64, wdpts::Int64)::Tuple{Array{Float64,2}, Array{Array{Int64,1},1}, Array{Array{Int64,1},1}}

	@assert npts*mpts == length(b[1,:]) ("ERROR: there are not npts * mpts = $(npts) * $(mpts) = $(npts * mpts) points of the net in b[]")

	ustep::Float64 = 1 / (udpts-1)
	wstep::Float64 = 1 / (wdpts - 1)
	n::Int64 = npts - 1
	m::Int64 = mpts - 1

	#=== x, y, z components of B ===#
	x::Array{Float64,2} = zeros(npts, mpts)
	y::Array{Float64,2} = zeros(npts, mpts)
	z::Array{Float64,2} = zeros(npts, mpts)

	for i = 1 : npts
		for j = 1 : mpts
			ij = (i - 1) * mpts + j
			x[i, j] = b[1, ij]
			y[i, j] = b[2, ij]
			z[i, j] = b[3, ij]
		end
	end

	#=== SYNC BLOCK BEGIN ===#
	#=== all the following @async are totally indipendent blocks ===#
	@sync begin
	
		#=== N = C * D ===#
		C::Array{Float64,2} = zeros(2, 2)
		D::Array{Float64,2} = zeros(2, 2)

		@async begin
			C = zeros(npts, npts)
			D = zeros(npts, npts)
			for i = 0 : n
				for j = 0 : (n - i)
					if (n-j-i)%2==0
						C[i+1, j+1] = binomial(n-j, n-j-i)
					else
						C[i+1, j+1] = 0.0-binomial(n-j, n-j-i)
					end
				end
			end

			for i = 1 : npts
				D[i, i] = abs(C[i, 1])
			end
		end

		#=== M = E * F ===#
		E::Array{Float64,2} = zeros(2, 2)
		F::Array{Float64,2} = zeros(2, 2)

		@async begin
			E = zeros(mpts, mpts)
			F = zeros(mpts, mpts)
			for i = 0 : m
				for j = 0 : (m - i)
					if (m - j - i) % 2 == 0
						E[i+1, j+1] = binomial(m-j, m-j-i)
					else
						E[i+1, j+1] = 0.0 - binomial(m-j, m-j-i)
					end
				end
			end

			for i = 1 : mpts
				F[i, i] = abs(E[i, 1])
			end
		end

		#=== U = [u^n ... u 1] ===#
		U::Array{Float64,2} = zeros(2, 2)

		@async begin
			U = ones(udpts, npts)
			u::Float64 = 0.0
			for i = 1 : udpts
				for j = 0 : n-1  #Per j=0 il valore è 1 ed è già impostato
					U[i, n-j] = U[i, npts-j] * u
				end
				u = u + ustep
			end
		end

		#=== W = [w^m ... w 1] ===#
		W::Array{Float64,2} = zeros(2, 2)

		@async begin
			W = ones(wdpts, mpts)
			w::Float64 = 0.0
			for i = 1 : wdpts
				for j = 0 : m-1  #Per j=0 il valore è 1 ed è già impostato
					W[i, m-j] = W[i, mpts-j] * w
				end
				w = w + wstep
			end
		end

		#=== Edge Vector ===#
		EV::Array{Array{Int64,1},1} = Array{Int64,1}[] #1-dimensional cellular complex used for plotting the curve with Plasm package

		@async begin
			for i = 1 : ((udpts * wdpts) - 1)
				app = Array{Int64}(2)
				app[1] = i
				app[2] = i+1
				push!(EV, app) 
			end

			for i = 1 : (wdpts * (udpts - 1))
				app = Array{Int64}(2)
				app[1] = i
				app[2] = i + wdpts
				push!(EV, app) 
			end
		end

		#=== Faces Vector ===#
		FV::Array{Array{Int64,1},1} = Array{Int64,1}[] #1-dimensional cellular complex used for plotting the curve with Plasm package

		@async begin
			for j = 0 : (udpts - 2)
				for i = 1 : wdpts - 1
					app = Array{Int64}(4)
					app[1] = j * wdpts + i
					app[2] = j * wdpts + i + 1
					app[3] = (j + 1) * wdpts + i
					app[4] = (j + 1) * wdpts + i + 1
					push!(FV, app) 
				end
			end
		end
	
	end
	#=== SYNC BLOCK END ===#

	#=== Points Vector ===#
	UCD::Array{Float64,2} = U * C * D
	FEW::Array{Float64,2} = transpose(W * E * F)
	xP::Array{Float64,2} = transpose(UCD * x * FEW)
	yP::Array{Float64,2} = transpose(UCD * y * FEW)
	zP::Array{Float64,2} = transpose(UCD * z * FEW)

	P::Array{Float64,2} = zeros(3, udpts * wdpts)

	for i = 1 : udpts
		for j = 1 : wdpts
			ij::Int64 = (i-1) * udpts + j
			P[1, ij] = xP[i, j]
			P[2, ij] = yP[i, j]
			P[3, ij] = zP[i, j]
		end
	end

	return (P, EV, FV)
end
