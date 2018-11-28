"""
basis(ord,t,npts,x)

# Arguments

- `ord::Int64` : order of the bspline basis function.
- `t::Float64` : parameter value
- `npts::Int64` : the number of the contol polygon vertices.
- `x::Array{Float64}` : knot vector

"""


function basis(ord::Int64,t::Float64,npts::Int64,x::Array{Float64})::Array{Float64}
    
    @assert length(x) == npts + ord ("ERROR: incompatible knot vector")

    m::Int64 = npts + ord
    N::Array{Float64} = zeros(m)
    
    #calculate the first order basis functions
    for i = 1 : m - 1
        if t >= x[i] && t < x[i+1]
            N[i] = 1
        else
            N[i] = 0
        end
    end
    
    #calculate the higher order basis functions
    for k = 2 : ord
        for i = 1 : m - k
            if N[i] != 0
                d = ((t - x[i]) * N[i]) / (x[i+k-1] - x[i])
            else
                d = 0
            end
            if N[i+1] != 0
                e = ((x[i+k] - t) * N[i+1]) / (x[i+k] - x[i+1])
            else
                e = 0
            end
            N[i] = d + e
        end
    end
    
    if t == x[m]
        N[npts] = 1
    end
    
    return N
end


#-------------------------------------------------------------------------------------


"""
dbasis(ord,t,npts,x)

# Arguments

- `ord::Int64` : order of the bspline basis function.
- `t::Float64` : parameter value
- `npts::Int64` : the number of the contol polygon vertices.
- `x::Array{Float64}` : knot vector

"""

function dbasis(ord::Int64,t::Float64,npts::Int64,x::Array{Float64})::Tuple{Array{Float64},Array{Float64},Array{Float64}}
    
    m::Int64 = npts + ord
    N::Array{Float64} = zeros(m)
    D1::Array{Float64} = zeros(m)
    D2::Array{Float64} = zeros(m)
   
    #calculate the first order basis functions
    for i = 1 : m - 1
        if t >= x[i] && t < x[i+1]
            N[i] = 1
        else
            N[i] = 0
        end
    end
    
    if t == x[m]
        N[npts] = 1
    end
        
    #calculate the higher order basis functions, first and second derivatve
    for k = 2 : ord
        for i = 1 : m - k
            if N[i] != 0
                b1 = ((t - x[i]) * N[i]) / (x[i+k-1] - x[i])
                f1 = N[i] / (x[i+k-1] - x[i])
            else
                b1 = 0
                f1 = 0
            end
            if N[i+1] != 0
                b2 = ((x[i+k] - t) * N[i+1]) / (x[i+k] - x[i+1])
                f2 = - N[i+1] / (x[i+k] - x[i+1])
            else
                b2 = 0
                f2 = 0
            end                
            if D1[i] != 0
                f3 = ((t - x[i]) * D1[i]) / (x[i+k-1] - x[i])
                s1 = (2 * D1[i]) / (x[i+k-1] - x[i])
            else
                f3 = 0
                s1 = 0
            end
            if D1[i+1] != 0
                f4 = ((x[i+k] - t) * D1[i+1]) / (x[i+k] - x[i+1])
                s2 = (-2 * D1[i+1]) / (x[i+k] - x[i+1])
            else
                f4 = 0
                s2 = 0
            end
            if D2[i] != 0
                s3 = ((t - x[i]) * D2[i]) / (x[i+k-1] - x[i])
            else
                s3 = 0
            end
            if D2[i+1] != 0
                s4 = ((x[i+k] - t) * D2[i+1]) / (x[i+k] - x[i+1])
            else
                s4 = 0
            end
                    
            N[i] = b1 + b2
            D1[i] = f1 + f2 + f3 + f4
            D2[i] = s1 + s2 + s3 + s4
        end
    end
        
    return N, D1, D2
    
end


#------------------------------------------------------------------------------------


"""
dbasisu(ord,t,npts,x)

# Arguments

- `ord::Int64` : order of the bspline basis function.
- `t::Float64` : parameter value
- `npts::Int64` : the number of the contol polygon vertices.
- `x::Array{Float64}` : knot vector

"""

function dbasisu(ord::Int64,t::Float64,npts::Int64,x::Array{Float64})::Tuple{Array{Float64},Array{Float64},Array{Float64}}
    
    m::Int64 = npts + ord
    N::Array{Float64} = zeros(m)
    D1::Array{Float64} = zeros(m)
    D2::Array{Float64} = zeros(m)
    
    #calculate the first order basis functions
    for i = 1 : m - 1
        if t >= x[i] && t < x[i+1]
            N[i]=1
        else
            N[i]=0
        end
    end
    
    if t == x[npts+1]
        N[npts]=1
        N[npts+1]=0  #the only difference with dbasis function
    end
    
    #calculate the higher order basis function,first and second derivatives
    for k = 2 : ord
        for i = 1 : m - k
            if N[i] != 0
                b1 = ((t - x[i]) * N[i]) / (x[i+k-1] - x[i])
                f1 = N[i] / (x[i+k-1] - x[i])
            else
                b1 = 0
                f1 = 0
            end
            if N[i+1] != 0
                b2 = ((x[i+k] - t) * N[i+1]) / (x[i+k] - x[i+1])
                f2 = -N[i+1] / (x[i+k] - x[i+1])
            else
                b2 = 0
                f2 = 0
            end
            if D1[i] != 0
                f3 = ((t - x[i]) * D1[i]) / (x[i+k-1] - x[i])
                s1 = (2 * D1[i]) / (x[i+k-1] - x[i])
            else
                f3 = 0
                s1 = 0
            end
            if D1[i+1] != 0
                f4 = ((x[i+k] - t) * D1[i+1]) / (x[i+k] - x[i+1])
                s2 = (-2 * D1[i+1]) / (x[i+k] - x[i+1])
            else
                f4 = 0
                s2 = 0
            end
            if D2[i] != 0
                s3 = ((t - x[i]) * D2[i]) / (x[i+k-1] - x[i])
            else
                s3 = 0
            end
            if D2[i+1] != 0
                s4 = ((x[i+k] - t) * D2[i+1]) / (x[i+k] - x[i+1])
            else
                s4 = 0
            end

            N[i] = b1 + b2
            D1[i] = f1 + f2 + f3 + f4
            D2[i] = s1 + s2 + s3 + s4
        end
    end
    
    return N, D1, D2

end