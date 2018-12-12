export knot, knotc, knotu


"""
    knot(npts, ord, centre=false, step=1.0)
Generate a `ord` order B-spline open uniform knot vector over `npts = n+1` control vertices.

An open uniform knot vector has multiplicity of knot values
 at the ends equals to the order `ord` of the B-Spline basis function.
Internal knot values are evenly spaced (uniform).
The resulting open uniform basis functions yield curves
 that behave most nearly like Bézier curves.
`centre` Boolean value indicates if the knots are centered in the
 zero value or if them goes from zero to `m-1` where `m=npts+ord`.
 If `true` is choosen, then `m` must be odd.
Necessary condition is that `npts>=ord` in order to correct build the knot vector.

---

# Arguments
- `npts::Int64`: the number of the `n+1` points of the controll polygon.
- `ord::Int64`: the order of the B-Spline (`degree k-1`).
- `centre::Bool=false`: if the knot is centered in zero.
- `step::Float64=1.0`: the step between the knots.

---

# Examples
```jldoctest
julia> knot(5, 2)
7-element Array{Float64,1}:
 0.0
 0.0
 1.0
 2.0
 3.0
 4.0
 4.0
```

```jldoctest
julia> knot(6, 3, false, 0.25)
9-element Array{Float64,1}:
 0.0
 0.0
 0.0
 0.25
 0.5
 0.75
 1.0
 1.0
 1.0
```

```jldoctest
julia> knot(5, 2, true)
7-element Array{Float64,1}:
 -2.0
 -2.0
 -1.0
 0.0
 1.0
 2.0
 2.0
```

---

_By Elia Onofri, Giuseppe Santorelli_
"""

function knot(npts::Int64, ord::Int64, centre::Bool=false, step::Float64=1.0)::Array{Float64}
    # Input Data Check
    @assert npts >= ord ("ERROR: ord > npts")
	
    local m::Int64 = npts+ord      		   # knot vector dimension
    local backstep::Float64 = 0.0   		   # backstep is the affine transition parameter for 0-central knots
    local x::Array{Float64} = Array{Float64}(m)    # knot vector output
	
    # Different preparation for centre knot vector
    if centre
        @assert m % 2 == 1 ("ERROR: no centring points available")
        backstep = step * (npts - ord + 1) / 2
    end
    
    x[1] = 0 - backstep
    for i = 2 : m
        if i > ord && i <= npts + 1
            x[i] = x[i-1] + step
        else
            x[i] = x[i-1]
        end
    end
    return x
end


#--------------------------------------------------------------------------------------------------------------------------------


"""
knotc(npts, ord, b)

Generate a nonuniform open knot vector proportional to the chord lengths between defining polygon vertices.

# Arguments

- `npts::Int64`: the number of the contol polygon vertices.
- `ord::Int64`: order of the basis function.
- `b::Array{Float64}`: 2-dimensional array containing the contol polygon vertices, every column represents a point.

---

# Examples

- Example 3.5 pag. 67

```jldoctest
julia> b = [0. 2 4 6 8; 0 6 3 6 6; 0 0 0 0 0]
3×5 Array{Float64,2}:
 0.0  2.0  4.0  6.0  8.0
 0.0  6.0  3.0  6.0  6.0
 0.0  0.0  0.0  0.0  0.0

julia> knotc(5,3,b)
8-element Array{Float64,1}:
 0.0    
 0.0    
 0.0    
 1.45338
 2.38171
 3.0    
 3.0    
 3.0    
```

---
_Giuseppe Santorelli_
"""

function knotc(npts::Int64,ord::Int64,b::Array{Float64,2})::Array{Float64}
    
    @assert length(b) == 3 * npts ("ERROR: array b not matching with parameters")
    @assert npts >= ord ("ERROR: ord > npts")

    m::Int64 = npts + ord
    x::Array{Float64} = Array{Float64}(m)
    chord::Array{Float64} = Array{Float64}(npts - 1)
    
    maxchord = 0
    icount = 0
    
    # determine chord distance between defining polygon vertices and their sum
    for i = 2 : npts
        icount = icount + 1
        xchord = b[1,i] - b[1,i - 1]
        ychord = b[2,i] - b[2,i - 1]
        zchord = b[3,i] - b[3,i - 1]
        chord[icount] = sqrt(xchord * xchord + ychord * ychord + zchord * zchord)
        maxchord = maxchord + chord[icount]
    end
    
    # multiplicity of c=order zeros at the beginning of the open knot vector
    for i = 1 : ord
        x[i] = 0
    end
    
    # generate the internal knot values
    csum = 0
    for i = 1 : npts - ord
        csum = csum + chord[i]
        numerator = i / (npts - ord + 1) * chord[i + 1] + csum
        x[ord + i] = (numerator / maxchord) * (npts - ord + 1)
    end
    
    # multiplicity of c=order zeros at the end of the open knot vector
    for i = npts + 1 : m
        x[i] = npts - ord + 1
    end
    
    return x
    
end



#--------------------------------------------------------------------------------------------------------------------------------


"""
    knotu(npts, ord, step=1.0)
Generate a B-spline periodic uniform knot vector over `npts = n+1` control vertices.

In a periodic uniform knot vector the knot values are evenly spaced by `step`.
Moreover each base is a translate of the others as it could be defined like:

```math
    N_{i, k}(t) = N_{i-1, k}(t-1) = N_{i+1, k}(t+1)
```

---

# Arguments
- `npts::Int64`: the number of the `n+1` points of the controll polygon.
- `ord::Int64`: the order of the B-Spline (`deg = ord-1`).
- `step::Float64=1.0`: the step between the knots.

---

# Example
```jldoctest
julia> knotu(5, 2)
7-element Array{Float64,1}:
 0.0
 1.0
 2.0
 3.0
 4.0
 5.0
 6.0
```

```
julia> knotu(5, 4, 0.25)
9-element Array{Float64,1}:
 0.0 
 0.25
 0.5 
 0.75
 1.0 
 1.25
 1.5 
 1.75
 2.0
```

---
_By Elia Onofri, Giuseppe Santorelli_
"""

function knotu(npts::Int64, ord::Int64, step::Float64=1.0)::Array{Float64}
    # Input Data Check
    @assert npts >= ord ("ERROR: ord > npts")
	
    local m::Int64 = npts + ord      		   # knot vector dimension
    local x::Array{Float64} = Array{Float64}(m)    # knot vector output

    for i = 1 : m
        x[i] = (i-1) * step
    end

    return x
end
