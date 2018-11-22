export knot, knotc, knotu


"""
    knot(n1, k, centre=false, step=1.0)
Generate a `k` order B-spline open uniform knot vector over `n1 = n+1` control vertices.

An open uniform knot vector has multiplicity of knot values
 at the ends equals to the order `k` of the B-Spline basis function.
Internal knot values are evenly spaced (uniform).
The resulting open uniform basis functions yield curves
 that behave most nearly like Bézier curves.
`centre` Boolean value indicates if the knots are centered in the
 zero value or if them goes from zero to `m-1` where `m=n1+k`.
 If `true` is choosen, then `m` must be odd.
Necessary condition is that `n1>=k` in order to correct build the knot vector.

---

# Arguments
- `n1::Int64`: the number of the `n+1` points of the controll polygon.
- `k::Int64`: the order of the B-Spline (`degree k-1`).
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

_By Elia Onofri_
"""

function knot(n1::Int64, k::Int64, centre::Bool=false, step::Float64=1.0)::Array{Float64}
    local x = Float64[]                # knot vector output
    local m::Int64 = n1+k              # knot vector dimension
    local backstep::Float64 = 0.0      # backstep is the affine transition parameter for 0-central knots
    
    # Different preparation for centre knot vector
    if centre
        if m%2!=1
            print("ERROR: no centring points available")
            return x
        else
            backstep = step*(n1-k+1)/2
        end
    end
    
    # Local data check
    if n1<k
        print("WARNING: k > n1")
    end
    
    push!(x, 0-backstep)
    for i = 2:m
        if i>k && i<=n1+1
            push!(x, last(x)+step)
        else
            push!(x, last(x))
        end
    end
    return x
end


#--------------------------------------------------------------------------------------------------------------------------------


"""
	knotc(n, c)

Generate a nonuniform open knot vector proportional to the chord lengths between defining polygon vertices.
"""

function knotc(n::Int64, c::Int64)::Array{Float64}
    [...]
end


#--------------------------------------------------------------------------------------------------------------------------------


"""
    knotu(n1, k, step=1.0)
Generate a B-spline periodic uniform knot vector over `n` control vertices.

In a periodic uniform knot vector the knot values are evenly spaced by `step`.
Moreover each base is a translate of the others as it could be defined like:

```math
    N_{i, k}(t) = N_{i-1, k}(t-1) = N_{i+1, k}(t+1)
```

---

# Arguments
- `n1::Int64`: the number of the `n+1` points of the controll polygon.
- `k::Int64`: the order of the B-Spline (`degree k-1`).
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
_By Elia Onofri_
"""

function knotu(n1::Int64, k::Int64, step::Float64=1.0)::Array{Float64}
    x=Float64[];          # knot vector output
    m=n1+k;             # knot vector dimension

    # Local data check
    if n1<k
        print("WARNING: k > n1")
    end

    for i=1:m
        push!(x, (i-1)*step);
    end

    return x;
end