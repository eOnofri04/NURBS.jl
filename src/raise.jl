export raise23, raise45

"""

raise23(b,x,npts)

Subroutine that increases a second-order B-spline curve to a third-order B-spline curve.
Returns the new control polygon vertices array and the number of its points.

# Arguments

- `b::Array{Float64}`: 2-dimensional array containing the contol polygon vertices
- `x::Array{Float64}`: knot vector
- `npts::Int64`: number of control polygon veritces

# Examples

- Example 3.10 pag 98

```jldoctest
julia> x = knot(4,2)
6-element Array{Float64,1}:
 0.0
 0.0
 1.0
 2.0
 3.0
 3.0

julia> b = [1. 2 3 4;1 2 2 1;1 1 1 1]
3×4 Array{Float64,2}:
 1.0  2.0  3.0  4.0
 1.0  2.0  2.0  1.0
 1.0  1.0  1.0  1.0

julia> raise23(b,x,4)[1]
3×7 Array{Float64,2}:
 1.0  1.5  2.0  2.5  3.0  3.5  4.0
 1.0  1.5  2.0  2.0  2.0  1.5  1.0
 1.0  1.0  1.0  1.0  1.0  1.0  1.0

julia> raise23(b,x,4)[2]
 7
```
"""

function raise23(b::Array{Float64}, x::Array{Float64}, npts::Int64)::Tuple{Array{Float64,2},Int64}

    @assert length(b) == 3*npts ("ERROR: number of control polygon vertices and basis vector length not compatible")
    @assert length(x) == npts + 2 ("ERROR: number of control polygon vertices and knot vector length not compatible")
    
    mpts::Int64 = 1                     #number of vertices in the new array b* 
    bstar::Array{Float64} = Float64[]
    
    if x[1] < x[2]                      #first point of b*
        push!(bstar,b[1,1]/2)
        push!(bstar,b[2,1]/2)
        push!(bstar,b[3,1]/2)
        mpts = mpts + 1
    end
    
    push!(bstar,b[1,1])                   #second point of b* (or first if x1 >= x2) 
    push!(bstar,b[2,1])
    push!(bstar,b[3,1])
    
    for i = 2 : npts                    #filling b* array
        if x[i] < x[i+1]
            mpts = mpts + 1
            push!(bstar,(b[1,i] + b[1,i-1]) / 2)
            push!(bstar,(b[2,i] + b[2,i-1]) / 2)
            push!(bstar,(b[3,i] + b[3,i-1]) / 2)
        end
        
        mpts = mpts + 1
        push!(bstar,b[1,i])
        push!(bstar,b[2,i])
        push!(bstar,b[3,i])
    end
    
    if x[npts + 2] > x[npts + 1]
        mpts = mpts + 1
        push!(bstar,b[1,npts] / 2)
        push!(bstar,b[2,npts] / 2)
        push!(bstar,b[3,npts] / 2)
    end
    
    bstar = reshape(bstar,3,:)       #transforms 1-dimensional bstar array in a 3 x mpts array 
    
    return bstar,mpts
        
end

#------------------------------------------------------------------------------------


"""

raise45(b,x,npts)

Subroutine that increases a cubic B-spline curve to a quartic B-spline curve.
Returns the new control polygon vertices array and the number of its points.

# Arguments

- `b::Array{Float64}`: 2-dimensional array containing the contol polygon vertices
- `x::Array{Float64}`: knot vector
- `npts::Int64`: number of control polygon veritces

# Examples

- Example 3.11 pag 100 
```jldoctest
julia> x = knot(4,4)
8-element Array{Float64,1}:
 0.0
 0.0
 0.0
 0.0
 1.0
 1.0
 1.0
 1.0

julia> b = [1. 2 3 4 ; 1 2 2 1 ; 1 1 1 1]
3×4 Array{Float64,2}:
 1.0  2.0  3.0  4.0
 1.0  2.0  2.0  1.0
 1.0  1.0  1.0  1.0

julia> raise45(b,x,4)[1]
3×5 Array{Float64,2}:
 1.0  1.75  2.5  3.25  4.0
 1.0  1.75  2.0  1.75  1.0
 1.0  1.0   1.0  1.0   1.0
```
- Example 3.12 pag 100
```jldoctest
julia> k = knot(7,4)
11-element Array{Float64,1}:
 0.0
 0.0
 0.0
 0.0
 1.0
 2.0
 3.0
 4.0
 4.0
 4.0
 4.0

julia> b = [1. 2 3 4 5 6 7; 1 2 2 1 0 0 1; 1 1 1 1 1 1 1]
3×7 Array{Float64,2}:
 1.0  2.0  3.0  4.0  5.0  6.0  7.0
 1.0  2.0  2.0  1.0  0.0  0.0  1.0
 1.0  1.0  1.0  1.0  1.0  1.0  1.0

julia> raise45(b,k,7)[1] #the result in position (2,8) is 2/24 and it is correct. On the book it is 20/24 (WRONG!)
3×11 Array{Float64,2}:
 1.0  1.75  2.25  2.95833  3.5  4.0  4.5  5.04167    5.75  6.25  7.0
 1.0  1.75  2.0   1.91667  1.5  1.0  0.5  0.0833333  0.0   0.25  1.0
 1.0  1.0   1.0   1.0      1.0  1.0  1.0  1.0        1.0   1.0   1.0
```
"""
function raise45(b::Array{Float64}, x::Array{Float64}, npts::Int64)::Tuple{Array{Float64,2},Int64}

    @assert length(b) == 3 * npts ("ERROR: number of control polygon vertices and basis vector length not compatible")
    @assert length(x) == npts + 4 ("ERROR: number of control polygon vertices and knot vector length not compatible")
    #mpts = npts + s + ord + 1 = npts + 5 + s ?????
    
    mpts::Int64 = 0
    bstar::Array{Float64} = Float64[]
    
    #coordinates of the first two points in b*
    for i = 1 : 3
        push!(bstar,b[i,1]) #b*[1]
    end
    for  i = 1 : 3
        push!(bstar,(b[i,1] + 3 * b[i,2]) / 4) #[b*[2]]
    end
    mpts = mpts + 2
    
    if npts == 4
        for i = 1 : 3
            push!(bstar,(b[i,2] + b[i,3]) / 2) #b*[3]
        end
        mpts = mpts + 1
    else
        for i = 1 : 3
            push!(bstar,(3 * b[i,2] + b[i,3]) / 4) #b*[3]
        end
        mpts = mpts + 1
        
        if npts == 5
            for i = 1 : 3
                push!(bstar,(b[i,2] + 6 * b[i,3] + b[i,4]) / 8) #b*[4]
            end
            mpts = mpts + 1
        else
            for i = 1 : 3
                push!(bstar,(3 * b[i,2] + 19 * b[i,3] + 2 * b[i,4]) / 24) #b*[4]
            end
            for i = 1 : 3
                push!(bstar,(b[i,3] + b[i,4]) / 2) #b*[5]
            end
            mpts = mpts + 2
            
            if npts > 6
                n4 = npts - 4
                for i = 3 : n4
                    for j = 1 : 3
                        push!(bstar, (b[j,i] + 10 * b[j,i+1] + b[j,i+2]) / 12) #b*[2i]
                    end
                    for j = 1 : 3
                        push!(bstar, (b[j,i+1] + b[j,i+2]) / 2) #b*[2(i+1)]
                    end
                    mpts = mpts + 2
                end
            end
            
            for i = 1 : 3
                push!(bstar, (3 * b[i,npts-1] + 19 * b[i,npts-2] + 2 * b[i,npts-3]) / 24) #b*[mpts-3]
            end
            mpts = mpts + 1
        end
        
        for i = 1 : 3
            push!(bstar, (3 * b[i,npts-1] + b[i,npts-2]) / 4) #b*[mpts-2]
        end
        mpts = mpts + 1
    end
    
    for i = 1 : 3
        push!(bstar, (b[i,npts] + 3 * b[i,npts-1]) / 4) #b*[mpts-1]
    end
    for i = 1 : 3
        push!(bstar, b[i,npts]) #b*[mpts]
    end
    mpts = mpts + 2
    
    bstar = reshape(bstar,3,:)  #transforms 1-dimensional bstar array in a 3 x mpts array 
    
    return bstar, mpts
    
end