export basis, dbasis, dbasisu

"""
    basis(npts, ord, t, x[])

Generate a B-spline basis functions `N[]` for a given open knot vector `x[]`.

A B-Spline basis is a collection of functions of a parameter `t` wich form
 a basis for the vectorial space of functions. The conformation of this set
 higly depends on the choosing of the knots `x[]` the curve is bound to.
The basis is computed with _Cox-de Boor_ recursive function applied to the
 basis dependency tree in order to optimise the computation.

---

# Arguments
- `npts::Int64`: the number of the points of the controll polygon.
- `ord::Int64`: the order of the B-Spline (`deg = ord-1`).
- `t::Float64`: the parameter value of the parametric curve.
- `x::Array{Float64}`: the knot vector.

---

# Examples

```jldoctest
julia> basis(5,3,0.5,[0. 0 0 1 2 3 3 3])
5-element Array{Float64,1}:
 0.25 
 0.625
 0.125
 0.0  
 0.0  
```
_By Elia Onofri, Giuseppe Santorelli_
"""

function basis(npts::Int64, ord::Int64, t::Float64, x::Array{Float64})::Array{Float64}
   
    local m::Int64 = npts + ord	    # Dimension of the knot vector
    
    local N::Array{Float64} = zeros(m)   # Basis progressive vector
    
    # Local check of the knot vector correctness
    @assert length(x) == m ("ERROR: incompatibile knot vector with given parameters n+1 = $(npts), k = $(ord)")
    
    # Eval N_{i,1} for i = 1:max_B
    for i = 1 : m - 1
        if (t >= x[i]) && (t < x[i+1])
            N[i] = 1
        else
            N[i] = 0
        end
    end
    
    # Eval higher basis N_{i,deg} for deg = 2:ord and i = 1:max_B-deg
    for deg = 2 : ord
        for i = 1 : m - deg
            # Eval of the direct dependency
            if N[i] == 0
                ddep = 0.0
            else
                ddep = ((t - x[i]) * N[i]) / (x[i+deg-1] - x[i])
            end
            # Eval of the forward dependency
            if N[i+1] == 0
                fdep = 0.0
            else
                fdep = ((x[i+deg] - t) * N[i+1]) / (x[i+deg] - x[i+1])
            end
            # Collection of the dependencies
            N[i] = ddep + fdep
        end
    end
    
    # Otherwise last point is zero
    if t == x[m]
        N[npts] = 1
    end
    
    return N[1:npts]
end

#--------------------------------------------------------------------------------------------------------------------------------


"""
    dbasis(npts, ord, t, x[])

Generate a B-spline basis functions `N[]` and its derivatives `D1[], D2[]` for a given open knot vector `x[]`.

A B-Spline basis is a collection of functions of a parameter `t` wich form
 a basis for the vectorial space of functions. The conformation of this set
 higly depends on the choosing of the knots `x[]` the curve is bound to.
The basis is computed with _Cox-de Boor_ recursive function applied to the
 basis dependency tree in order to optimise the computation.

---

# Arguments
- `npts::Int64`: the number of the points of the controll polygon.
- `ord::Int64`: the order of the B-Spline (`deg = ord-1`).
- `t::Float64`: the parameter value of the parametric curve.
- `x::Array{Float64}`: the knot vector.

---

# Examples

```jldoctest
julia> dbasis(5,3,0.5,[0. 0 0 1 2 3 3 3])[1]
5-element Array{Float64,1}:
 0.25 
 0.625
 0.125
 0.0  
 0.0  

julia> dbasis(5,3,0.5,[0. 0 0 1 2 3 3 3])[2]
5-element Array{Float64,1}:
 -1.0
  0.5
  0.5
  0.0
  0.0

julia> dbasis(5,3,0.5,[0. 0 0 1 2 3 3 3])[3]
5-element Array{Float64,1}:
  2.0
 -3.0
  1.0
  0.0
  0.0
```

_By Elia Onofri, Giuseppe Santorelli_
"""

function dbasis(npts::Int64, ord::Int64, t::Float64, x::Array{Float64})::Tuple{Array{Float64}, Array{Float64}, Array{Float64}}
   
    local m::Int64 = npts + ord			# Dimension of the knot vector
    
    local N::Array{Float64} = zeros(m)   # Basis progressive vector
    local D1::Array{Float64} = zeros(m)  # First Derivative progressive vector
    local D2::Array{Float64} = zeros(m)  # Second Derivative progressive vector
    
    # Local check of the knot vector correctness
    @assert length(x) == m ("ERROR: incompatibile knot vector with given parameters n+1 = $(npts), k = $(ord)")
    
    # Eval N_{i,1} for i = 1:max_B
    for i = 1 : m - 1
        if (t >= x[i]) && (t < x[i+1])
            N[i] = 1
        else
            N[i] = 0
        end
    end
    
    if t == x[m]    # last(x)
        N[npts] = 1
    end
    
    # Eval higher basis N_{i,deg} for deg = 2:ord and i = 1:max_B-deg
    for deg = 2 : ord
        for i = 1 : m - deg
            
            # ----Eval of the direct dependency----
            if N[i] == 0
                ddep = 0.0
                f1 = 0.0
            else
                ddep = ((t - x[i]) * N[i]) / (x[i+deg-1] - x[i])
                f1 = N[i] / (x[i+deg-1] - x[i])
            end
            
            # ----Eval of the forward dependency----
            if N[i+1] == 0
                fdep = 0.0
                f2 = 0.0
            else
                fdep = ((x[i+deg] - t) * N[i+1]) / (x[i+deg] - x[i+1])
                f2 = -N[i+1] / (x[i+deg] - x[i+1])
            end

            #----Eval of the direct first derivative dependency----
            if D1[i] != 0    # if the lower order basis function is zero skip the calculation
                f3 = ((t - x[i]) * D1[i]) / (x[i+deg-1] - x[i])
                s1 = (2 * D1[i]) / (x[i+deg-1] - x[i])
            else
                f3 = 0.0
                s1 = 0.0
            end
            
            #----Eval of the forward first derivative dependency----
            if D1[i+1] != 0    # if the lower order basis function is zero skip the calculation
                f4 = ((x[i+deg] - t) * D1[i+1]) / (x[i+deg] - x[i+1])
                s2 = ((-2) * D1[i+1]) / (x[i+deg] - x[i+1])
            else
                f4 = 0.0
                s2 = 0.0
            end
		
            #----Eval of the direct second derivative dependency----
            if D2[i] != 0    # if the lower order basis function is zero skip the calculation
                s3 = ((t - x[i]) * D2[i] ) / (x[i+deg-1]-x[i])
            else
                s3 = 0.0
            end
            
            #----Eval of the forward second derivative dependency----
            if D2[i+1] != 0    # if the lower order basis function is zero skip the calculation
                s4 = ((x[i+deg] - t) * D2[i+1]) / (x[i+deg] - x[i+1])
            else
                s4 = 0.0
            end
			
            # Collection of the dependencies
            N[i] = ddep + fdep
            D1[i] = f1 + f2 + f3 + f4
            D2[i] = s1 + s2 + s3 + s4
        end
    end
    
    return (N[1:npts], D1[1:npts], D2[1:npts])
end



#--------------------------------------------------------------------------------------------------------------------------------


"""
    dbasisu(npts, ord, t, x[])

Generate a B-spline basis functions `N[]` and its derivatives `D1[], D2[]` for a given uniform periodic knot vector `x[]`.

A B-Spline basis is a collection of functions of a parameter `t` wich form
 a basis for the vectorial space of functions. The conformation of this set
 higly depends on the choosing of the knots `x[]` the curve is bound to.
The basis is computed with _Cox-de Boor_ recursive function applied to the
 basis dependency tree in order to optimise the computation.

---

# Arguments
- `npts::Int64`: the number of the points of the controll polygon.
- `ord::Int64`: the order of the B-Spline (`deg = ord-1`).
- `t::Float64`: the parameter value of the parametric curve.
- `x::Array{Float64}`: the knot vector.

---

# Examples

```jldoctest
julia> dbasisu(5,3,0.5,[0. 1 2 3 4 5 6 7])[1]
5-element Array{Float64,1}:
 0.125
 0.0  
 0.0  
 0.0  
 0.0  

julia> dbasisu(5,3,0.5,[0. 1 2 3 4 5 6 7])[2]
5-element Array{Float64,1}:
 0.5
 0.0
 0.0
 0.0
 0.0

julia> dbasisu(5,3,0.5,[0. 1 2 3 4 5 6 7])[3]
5-element Array{Float64,1}:
 1.0
 0.0
 0.0
 0.0
 0.0
```

_By Elia Onofri, Giuseppe Santorelli_
"""

function dbasisu(npts::Int64, ord::Int64, t::Float64, x::Array{Float64})::Tuple{Array{Float64}, Array{Float64}, Array{Float64}}
    local m::Int64 = npts+ord			# Dimension of the knot vector
    local N::Array{Float64} = zeros(m)   # Basis progressive vector
    local D1::Array{Float64} = zeros(m)  # First Derivative progressive vector
    local D2::Array{Float64} = zeros(m)  # Second Derivative progressive vector
    
    # Local check of the knot vector correctness
    @assert length(x) == m ("ERROR: incompatibile knot vector with given parameters n+1 = $(npts), k = $(ord)")
    
    # Eval N_{i,1} for i = 1:max_B
    for i = 1 : m-1
        if (t >= x[i]) && (t < x[i+1])
            N[i] = 1
        else
            N[i] = 0
        end
    end
    
    if t == x[npts+1]    # handle the end specially
        N[npts] = 1   	 # resetting the first order basis functions.
        N[npts+1] = 0
    end
    
   
    # Eval higher basis N_{i,deg} for deg = 2:ord and i = 1:max_B-deg
    for deg = 2 : ord
        for i = 1 : m-deg
            
            # ----Eval of the direct dependency----
            if N[i] == 0
                ddep = 0.0
                f1 = 0.0
            else
                ddep = ((t - x[i]) * N[i]) / (x[i+deg-1] - x[i])
                f1 = N[i] / (x[i+deg-1] - x[i])
            end
            
            # ----Eval of the forward dependency----
            if N[i+1] == 0
                fdep = 0.0
                f2 = 0.0
            else
                fdep = ((x[i+deg] - t) * N[i+1]) / (x[i+deg] - x[i+1])
                f2 = -N[i+1] / (x[i+deg] - x[i+1])
            end

            #----Eval of the direct first derivative dependency----
            if D1[i] != 0    # if the lower order basis function is zero skip the calculation
                f3 = ((t - x[i]) * D1[i]) / (x[i+deg-1] - x[i])
                s1 = (2 * D1[i]) / (x[i+deg-1] - x[i])
            else
                f3 = 0.0
                s1 = 0.0
            end
            
            #----Eval of the forward first derivative dependency----
            if D1[i+1] != 0    # if the lower order basis function is zero skip the calculation
                f4 = ((x[i+deg] - t) * D1[i+1]) / (x[i+deg] - x[i+1])
                s2 = ((-2) * D1[i+1]) / (x[i+deg] - x[i+1])
            else
                f4 = 0.0
                s2 = 0.0
            end
		
            #----Eval of the direct second derivative dependency----
            if D2[i] != 0    # if the lower order basis function is zero skip the calculation
                s3 = ((t - x[i]) * D2[i] ) / (x[i+deg-1]-x[i])
            else
                s3 = 0.0
            end
            
            #----Eval of the forward second derivative dependency----
            if D2[i+1] != 0    # if the lower order basis function is zero skip the calculation
                s4 = ((x[i+deg] - t) * D2[i+1]) / (x[i+deg] - x[i+1])
            else
                s4 = 0.0
            end
			
            # Collection of the dependencies
            N[i] = ddep + fdep
            D1[i] = f1 + f2 + f3 + f4
            D2[i] = s1 + s2 + s3 + s4
        end
    end
    
    return (N[1:npts], D1[1:npts], D2[1:npts])
end
