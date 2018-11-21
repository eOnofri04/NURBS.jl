export knot, knotc, knotu


"""
	knot(n, c)

Generate a B-spline open knot vector with multiplicity k at the ends.

An open uniform knot vector has multiplicity of knot values
 at the ends equals to the order `k` of the B-Spline basis function.
Internal knot values are evenly spaced.
The resulting open uniform basis functions yield curves
that behave most nearly like Bézier curves.

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

function knot(n::Int64,c::Int64)::Array{Int64}
    x=Int64[];
    nplusc::Int64 = n+c;
    nplus2::Int64 = n+2;
    push!(x, 0);
    for i = 2:nplusc
        if i>c && i<nplus2
            push!(x, last(x)+1);
        else
            push!(x, last(x));
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
Generate a B-spline periodic uniform knot vector over `n` control vertices.

In a periodic uniform knot vector the knot values are evenly spaced.
Moreover each base is a translate of the others as it could be defined like:
```tex
N_{i, k}(t) = N_{i-1, k}(t-1) = N_{i+1, k}(t+1)
```

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

function knotu(n::Int64, c::Int64)::Array{Int64}
    x=Int64[];
    nplusc=n+c;
    nplus2=n+2;
    for i=1:nplusc
        push!(x, i-1);
    end
    return x;
end