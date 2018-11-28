export knot, knotc, knotu


"""
	knot(n, c)

Generate a B-spline open knot vector with multiplicity k at the ends.

# Examples

```jldoctest
julia> knot(5, 2)
7-element Array{Int64,1}:
 0
 0
 1
 2
 3
 4
 4
```
"""

function knot(n,c)
    
    nplusc = n + c
    nplus2 = n + 2

    x::Array{Float64} = Array{Float64}(nplusc)

    x[1] = 0
    
    for i = 2:nplusc
        if i > c  &&  i < nplus2
            
            x[i] = x[i - 1] + 1
        
        else

            x[i] = x[i - 1]
        end
    end

    return x
end

#-----------------------------------------------------------------------

"""
	knotc(n, c)

Generate a nonuniform open knot vector proportional to the chord lengths between defining polygon vertices.
"""

function knotc(n::Int64, c::Int64)::Array{Float64}
    [...]
end

#-----------------------------------------------------------------------

"""
	knotu(n, c)

Generate a B-spline periodic uniform knot vector.
```jldoctest
julia> knotu(5, 2)
7-element Array{Int64,1}:
 0
 1
 2
 3
 4
 5
 6
```
"""

function knotu(n,c)

    nplusc = n + c
    nplus2 = n + 2

    x = Array{Float64}(nplusc)
    
    for i = 1:nplusc
        x[i] = i - 1
    end

    return x
end