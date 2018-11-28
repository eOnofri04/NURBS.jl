"""
bspline(npts,ord,p1,b)

# Arguments

- `npts::Int64` : the number of the contol polygon vertices.
- `ord::Int64` : order of the bspline basis function.
- `p1::Int64` : number of the points to be calculated on the curve
- `b::Array{Float64}` : Array containing the control polygon vertices

"""


function bspline(npts::Int64,ord::Int64,p1::Int64,b::Array{Float64})::Array{Float64}
    
    @assert length(b) == 3 * npts ("ERROR: array b not matching with parameters")
    @assert npts >= ord ("ERROR: ord > npts")

    m::Int64 = npts + ord
    P::Array{Float64} = zeros(3*p1)
    
    x::Array{Float64} = knot(npts,ord)
    
    icount::Int64 = 0
    t::Float64 = 0.0
    
    step::Float64 = x[m] / (p1-1)
    
    for i1 = 1 : p1
        if x[m] - t < 5e-6
            t = x[m]
        end
        
        nbasis::Array{Float64} = basis(ord,t,npts,x)
        
        for j = 1 : 3
            jcount = j
            P[icount+j] = 0
            for i = 1 : npts
                temp = nbasis[i] * b[jcount]
                P[icount+j] = P[icount+j] + temp
                jcount = jcount+3
            end
        end
        
        icount = icount + 3
        t = t + step
    end
    
    return P
end


#-------------------------------------------------------------------------------------


"""
bsplineu(npts,ord,p1,b)

# Arguments

- `npts::Int64` : the number of the contol polygon vertices.
- `ord::Int64` : order of the bspline basis function.
- `p1::Int64` : number of the points to be calculated on the curve
- `b::Array{Float64}` : Array containing the control polygon vertices

"""

function bsplineu(npts::Int64,ord::Int64,p1::Int64,b::Array{Float64})::Array{Float64}
    
    @assert length(b) == 3 * npts ("ERROR: array b not matching with parameters")
    @assert npts >= ord ("ERROR: ord > npts")

    m::Int64 = npts + ord
    P::Array{Float64} = zeros(3*p1)
    x::Array{Float64} = zeros(m)
    nbasis::Array{Float64} = zeros(npts) 
        
    x = knotu(npts,ord)
    
    icount::Int64 = 0
    t::Float64 = ord - 1
    step::Float64 = (npts - (ord - 1)) / (p1 - 1)
    
    for i1 = 1 : p1
        
        if (npts - t) < 5e-6
            t = npts
        end
        
        nbasis = basis(ord,t,npts,x)
        
        for j = 1 : 3
            jcount = j
            
            for i = 1 : npts
                temp = nbasis[i] * b[jcount]
                P[icount + j] = P[icount + j] + temp
                jcount = jcount + 3
            end
        end
        
        icount = icount + 3
        t = t + step
    end
    
    return P
    
end


#------------------------------------------------------------------------------------


"""
dbspline(npts,ord,p1,b)

# Arguments

- `npts::Int64` : the number of the contol polygon vertices.
- `ord::Int64` : order of the bspline basis function.
- `p1::Int64` : number of the points to be calculated on the curve
- `b::Array{Float64}` : Array containing the control polygon vertices

"""

function dbspline(npts::Int64,ord::Int64,p1::Int64,b::Array{Float64})::Tuple{Array{Float64},Array{Float64},Array{Float64}}
    
    @assert length(b) == 3 * npts ("ERROR: array b not matching with parameters")
    @assert npts >= ord ("ERROR: ord > npts")

    m = npts + ord

    P::Array{Float64} = zeros(3*p1)
    D1::Array{Float64} = zeros(3*p1)
    D2::Array{Float64} = zeros(3*p1)
    x::Array{Float6464} = zeros(m)
    nbasis::Array{Float64} = zeros(npts)
    d1nbasis::Array{Float64} = zeros(npts)
    d2nbasis::Array{Float64} = zeros(npts)
    
    x = knot(npts,ord)
    
    icount::Int64 = 0
    t::Float64 = 0.0
    step::Float64 = x[m] / (p1 - 1)
    
    for i1 = 1 : p1
        if x[nplusc] - t < 5e-6
            t = x[nplusc]
        end
        
        nbasis, d1nbasis, d2nbasis = dbasis(ord,t,npts,x)
        
        for j = 0 : 3
            jcount = j
           
            for i = 1 : npts
                P[icount + j] = P[icount + j] + nbasis[i] * b[jcount]
                D1[icount + j] = D1[icount + j] + d1nbasis[i] * b[jcount]
                D2[icount + j] = D2[icount + j] + d2nbasis[i] * b[jcount]
                jcount = jcount + 3
            end
        end
        
        icount = icount + 3
        t = t + step
    end
    
    return P, D1, D2
    
end


#------------------------------------------------------------------------------------


"""
dbspline(npts,ord,p1,b)

# Arguments

- `npts::Int64` : the number of the contol polygon vertices.
- `ord::Int64` : order of the bspline basis function.
- `p1::Int64` : number of the points to be calculated on the curve
- `b::Array{Float64}` : Array containing the control polygon vertices

"""

function dbsplineu(npts::Int64,ord::Int64,p1::Int64,b::Array{Float64})::Tuple{Array{Float64},Array{Float64},Array{Float64}}
    
    @assert length(b) == 3 * npts ("ERROR: array b not matching with parameters")
    @assert npts >= ord ("ERROR: ord > npts")

    m::Int64 = npts + ord

    P::Array{Float64} = zeros(3*p1)
    D1::Array{Float64} = zeros(3*p1)
    D2::Array{Float64} = zeros(3*p1)
    x::Array{Int64} = zeros(m)
    nbasis::Array{Float64} = zeros(npts)
    d1nbasis::Array{Float64} = zeros(npts)
    d2nbasis::Array{Float64} = zeros(npts)
        
    x = knotu(npts,ord)
    
    icount::Int64 = 0
    t::Float64 = ord - 1
    step::Float64 = (npts - (ord - 1)) / (p1 - 1)
    
    for i1 = 1 : p1
        if npts - t < 5e-6
            t = npts
        end
        
        nbasis, d1nbasis, d2nbasis = dbasisu(ord,t,npts,x)
        
        for j = 0 : 3
            jcount = j
            
            for i = 1 : npts
                P[icount + j] = P[icount + j] + nbasis[i] * b[jcount]
                D1[icount + j] = D1[icount + j] + d1nbasis[i] * b[jcount]
                D2[icount + j] = D2[icount + j] + d2nbasis[i] * b[jcount]
                jcount = jcount + 3
            end
        end
        
        icount = icount + 3
        t = t + step
    end
    
    return P, D1, D2
    
end