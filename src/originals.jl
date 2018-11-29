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


#--------------------------------------------------------------------------------------------------------------------------------

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

function knotu(n::Int64, c::Int64)::Array{Int64}
    x=Int64[];
    nplusc=n+c;
    nplus2=n+2;
    for i=1:nplusc
        push!(x, i-1);
    end
    return x;
end

#--------------------------------------------------------------------------------------------------------------------------------


"""
    basis(c, t, npts, x)

Generate B-spline basis functions for open uniform knot vectors.
"""

function basis(c::Int64, t::Float64, npts::Int64, x::Array{Int64})::Array{Float64}
    tempt = Array{Float64}(36)
    n = Float64[];
    nplusc = npts + c;
    
#    print(knot vector is );
#    for i=1:nplusc
#        print(& &, i, x[i]);
#    end
#    print(t is &, t);
    
    # calculate the first order basis functions n[i][1]
    
    for i=1:nplusc-1
        if (t>=x[i])&&(t<x[i+1])
            temp[i] = 1;
        else
            temp[i] = 0;
        end
    end
    
    # calculate the higher order basis functions
    for k=2:c
        for i=1:nplusc-k
            if temp[i] != 0    # if the lower order basis function is zero skip the calculation
                d = ((t-x[i])*temp[i])/(x[i+k-1]-x[i]);
            else
                d = 0;
            end
            if temp[i+1] != 0    # if the lower order basis function is zero skip the calculation
                e = ((x[i+k]-t)*temp[i+1])/(x[i+k]-x[i+1]);
            else
                e = 0;
            end
            temp[i] = d+e;
        end
    end
    
    if t == x[nplusc]    # pick the last point
        temp[npts] = 1;
    end
    
    # put in n array
    for i=1:npts
        push!(n, temp[i]);
    end
    return n;
end