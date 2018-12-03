export basis, dbasis, dbasisu

"""
    basis(ord, t, npts, x[])

Generate a B-spline basis functions `N[]` for a given open knot vector `x[]`.

A B-Spline basis is a collection of functions of a parameter `t` wich form
 a basis for the vectorial space of functions. The conformation of this set
 higly depends on the choosing of the knots `x[]` the curve is bound to.
The basis is computed with _Cox-de Boor_ recursive function applied to the
 basis dependency tree in order to optimise the computation.

---

# Arguments
- `ord::Int64`: the order of the B-Spline (`deg = ord-1`).
- `t::Float64`: the parameter value of the parametric curve.
- `npts::Int64`: the number of the points of the controll polygon.
- `x::Array{Float64}`: the knot vector.

---

_By Elia Onofri, Giuseppe Santorelli_
"""

function basis(ord::Int64, t::Float64, npts::Int64, x::Array{Float64})::Array{Float64}
    local max_N = npts-1+ord        # Needs i+(ord-1)|i=npts = npts-1+ord trivial basis
    local m::Int64 = npts+ord	    # Dimension of the knot vector
    local ddep::Float64 = 0.0       # Direct Dependency partial sum
    local fdep::Float64 = 0.0       # Forward Dependency partial sum
    
    local N::Array{Float64} = Array{Float64}(max_N)   # Basis progressive vector
    
    # Local check of the knot vector correctness
    @assert length(x) == m ("ERROR: incompatibile knot vector with given parameters n+1 = $(npts), k = $(ord)")
    
    # Eval N_{i,1} for i = 1:max_B
    for i = 1 : max_N
        if (t >= x[i]) && (t < x[i+1])
            N[i] = 1
        else
            N[i] = 0
        end
    end
    
    # Eval higher basis N_{i,deg} for deg = 2:ord and i = 1:max_B-deg
    for deg = 2 : ord
        for i = 1 : m-deg
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
    dbasis(ord, t, npts, x[])

Generate a B-spline basis functions `N[]` and its derivatives `d1[], d2[]` for a given open knot vector `x[]`.

A B-Spline basis is a collection of functions of a parameter `t` wich form
 a basis for the vectorial space of functions. The conformation of this set
 higly depends on the choosing of the knots `x[]` the curve is bound to.
The basis is computed with _Cox-de Boor_ recursive function applied to the
 basis dependency tree in order to optimise the computation.

---

# Arguments
- `ord::Int64`: the order of the B-Spline (`deg = ord-1`).
- `t::Float64`: the parameter value of the parametric curve.
- `npts::Int64`: the number of the points of the controll polygon.
- `x::Array{Float64}`: the knot vector.

---

_By Elia Onofri, Giuseppe Santorelli_
"""

function dbasis(ord::Int64, t::Float64, npts::Int64, x::Array{Float64})::Tuple{Array{Float64}, Array{Float64}, Array{Float64}}
    local max_N::Int64 = npts-1+ord		# Needs i+(ord-1)|i=npts = npts-1+ord trivial basis
    local m::Int64 = npts+ord			# Dimension of the knot vector
    local ddep::Float64 = 0.0   		# Direct Dependency partial sum
    local fdep::Float64 = 0.0   		# Forward Dependency partial sum
    
    local N::Array{Float64} = Array{Float64}(max_N)   # Basis progressive vector
    local d1::Array{Float64} = Array{Float64}(max_N)  # First Derivative progressive vector
    local d2::Array{Float64} = Array{Float64}(max_N)  # Second Derivative progressive vector
    
    # Local check of the knot vector correctness
    @assert length(x) == m ("ERROR: incompatibile knot vector with given parameters n+1 = $(npts), k = $(ord)")
    
    # Eval N_{i,1} for i = 1:max_B
    for i = 1 : max_N
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
            if d1[i] != 0    # if the lower order basis function is zero skip the calculation
                f3 = ((t - x[i]) * d1[i]) / (x[i+deg-1] - x[i])
                s1 = (2 * d1[i]) / (x[i+deg-1] - x[i])
            else
                f3 = 0.0
                s1 = 0.0
            end
            
            #----Eval of the forward first derivative dependency----
            if d1[i+1] != 0    # if the lower order basis function is zero skip the calculation
                f4 = ((x[i+deg] - t) * d1[i+1]) / (x[i+deg] - x[i+1])
                s2 = ((-2) * d1[i+1]) / (x[i+deg] - x[i+1])
            else
                f4 = 0.0
                s2 = 0.0
            end
		
            #----Eval of the direct second derivative dependency----
            if d2[i] != 0    # if the lower order basis function is zero skip the calculation
                s3 = ((t - x[i]) * d2[i] ) / (x[i+deg-1]-x[i])
            else
                s3 = 0.0
            end
            
            #----Eval of the forward second derivative dependency----
            if d2[i+1] != 0    # if the lower order basis function is zero skip the calculation
                s4 = ((x[i+deg] - t) * d2[i+1]) / (x[i+deg] - x[i+1])
            else
                s4 = 0.0
            end
			
            # Collection of the dependencies
            N[i] = ddep + fdep
            d1[i] = f1 + f2 + f3 + f4
            d2[i] = s1 + s2 + s3 + s4
        end
    end
    
    return (N[1:npts], d1[1:npts], d2[1:npts])
end


#--------------------------------------------------------------------------------------------------------------------------------


"""
    dbasisu(ord, t, npts, x[])

Generate a B-spline basis functions `N[]` and its derivatives `d1[], d2[]` for a given uniform periodic knot vector `x[]`.

A B-Spline basis is a collection of functions of a parameter `t` wich form
 a basis for the vectorial space of functions. The conformation of this set
 higly depends on the choosing of the knots `x[]` the curve is bound to.
The basis is computed with _Cox-de Boor_ recursive function applied to the
 basis dependency tree in order to optimise the computation.

---

# Arguments
- `ord::Int64`: the order of the B-Spline (`deg = ord-1`).
- `t::Float64`: the parameter value of the parametric curve.
- `npts::Int64`: the number of the points of the controll polygon.
- `x::Array{Float64}`: the knot vector.

---

_By Elia Onofri, Giuseppe Santorelli_
"""

function dbasis(ord::Int64, t::Float64, npts::Int64, x::Array{Float64})::Tuple{Array{Float64}, Array{Float64}, Array{Float64}}
    local max_N::Int64 = npts-1+ord		# Needs i+(ord-1)|i=npts = npts-1+ord trivial basis
    local m::Int64 = npts+ord			# Dimension of the knot vector
    local ddep::Float64 = 0.0   		# Direct Dependency partial sum
    local fdep::Float64 = 0.0   		# Forward Dependency partial sum
    
    local N::Array{Float64} = Array{Float64}(max_N)   # Basis progressive vector
    local d1::Array{Float64} = Array{Float64}(max_N)  # First Derivative progressive vector
    local d2::Array{Float64} = Array{Float64}(max_N)  # Second Derivative progressive vector
    
    # Local check of the knot vector correctness
    @assert length(x) == m ("ERROR: incompatibile knot vector with given parameters n+1 = $(npts), k = $(ord)")
    
    # Eval N_{i,1} for i = 1:max_B
    for i = 1 : max_N
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
            if d1[i] != 0    # if the lower order basis function is zero skip the calculation
                f3 = ((t - x[i]) * d1[i]) / (x[i+deg-1] - x[i])
                s1 = (2 * d1[i]) / (x[i+deg-1] - x[i])
            else
                f3 = 0.0
                s1 = 0.0
            end
            
            #----Eval of the forward first derivative dependency----
            if d1[i+1] != 0    # if the lower order basis function is zero skip the calculation
                f4 = ((x[i+deg] - t) * d1[i+1]) / (x[i+deg] - x[i+1])
                s2 = ((-2) * d1[i+1]) / (x[i+deg] - x[i+1])
            else
                f4 = 0.0
                s2 = 0.0
            end
		
            #----Eval of the direct second derivative dependency----
            if d2[i] != 0    # if the lower order basis function is zero skip the calculation
                s3 = ((t - x[i]) * d2[i] ) / (x[i+deg-1]-x[i])
            else
                s3 = 0.0
            end
            
            #----Eval of the forward second derivative dependency----
            if d2[i+1] != 0    # if the lower order basis function is zero skip the calculation
                s4 = ((x[i+deg] - t) * d2[i+1]) / (x[i+deg] - x[i+1])
            else
                s4 = 0.0
            end
			
            # Collection of the dependencies
            N[i] = ddep + fdep
            d1[i] = f1 + f2 + f3 + f4
            d2[i] = s1 + s2 + s3 + s4
        end
    end
    
    return (N[1:npts], d1[1:npts], d2[1:npts])
end
