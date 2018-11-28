"""
knot(npts,ord,center=false,step=1.0)

# Arguments

- `npts::Int64` : the number of the contol polygon vertices.
- `ord::Int64` : order of the basis function.
- `center::Bool` : if the vector is centered in zero.
- `step::Float64` : knot distance

"""

function knot(npts::Int64,ord::Int64,center::Bool=false,step::Float64=1.0)::Array{Float64}
    
    @assert npts >= ord ("ERROR: ord > npts")

    m::Int64 = npts + ord
    backstep::Float64 = 0.0
    x::Array{Float64} = Array{Float64}(m)
    
    if center
        @assert m % 2 == 1 ("ERROR: knot vector contains even number of points")
        backstep = step * (npts - ord + 1) / 2
    end

    x[1] = 0 - backstep
    
    for i = 2 : m
        if i > ord && i < npts + 2
            x[i] = x[i-1] + step
        else
            x[i] = x[i-1]
        end
    end
    
    return x
end



#-------------------------------------------------------------------------------------


"""
knotc(npts,ord,b)

Generate a nonuniform open knot vector proportional to the chord lengths between defining polygon vertices.

# Arguments

- `npts::Int64` : the number of the contol polygon vertices.
- `ord::Int64` : order of the basis function.
- `b::Array{Float64}` : array containing the contol polygon vertices.

"""

function knotc(npts::Int64,ord::Int64,b::Array{Float64})::Array{Float64}
    
    @assert length(b) == 3 * npts ("ERROR: array b not matching with parameters")
    @assert npts >= ord ("ERROR: ord > npts")

    chord::Array{Float64} = Array{Float64}(31)
    m::Int64 = npts + ord
    n::Int64 = npts - 1
    x::Array{Float64} = Array{Float64}(m)
        
    maxchord = 0
    icount = 0
    
    # determine chord distance between defining polygon vertices and their sum
    for i = 4 : 3 :3*npts
        icount = icount + 1
        xchord = b[i] - b[i - 3]
        ychord = b[i + 1] - b[i - 2]
        zchord = b[i + 2] - b[i - 1]
        chord[icount] = sqrt(xchord * xchord + ychord * ychord + zchord * zchord)
        maxchord = maxchord + chord[icount]
    end
    
    # multiplicity of c=order zeros at the beginning of the open knot vector
    for i = 1 : ord
        x[i] = 0
    end
    
    # generate the internal knot values
    for i = 1 : npts - ord
        csum = 0
        
        for j = 1 : i
            csum = csum + chord[j]
        end
        
        numerator = i / (npts - ord + 1) * chord[i + 1] + csum
        x[ord + i] = (numerator / maxchord) * (npts - ord + 1)
    end
    
    # multiplicity of c=order zeros at the end of the open knot vector
    for i = npts + 1 : m
        x[i] = npts - ord + 1
    end
    
    return x
    
end


#------------------------------------------------------------------------------------


"""
knotu(npts,ord,step=1.0)

# Arguments

- `npts::Int64` : the number of the contol polygon vertices.
- `ord::Int64` : order of the basis function.
- `step::Float64` : knot distance

"""


function knotu(npts::Int64,ord::Int64,step::Float64=1.0)::Array{Float64}
    
    @assert npts >= ord ("ERROR: ord > npts")

    m::Int64 = npts + ord
    x::Array{Float64} = Array{Int64}(m)
    
    x[1] = 0
    
    for i = 2 : m
        x[i] = x[i - 1] + step
    end

    return x
end