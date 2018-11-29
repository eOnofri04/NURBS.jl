export bspline, bsplineu, dbspline, dbsplineu


"""
	bsplfit(dpts, d, npts, k)

Fit a B-spline curve using an open uniform knot vector.
"""

function bsplfit(dpts::Int64, d::Array{Float64}, npts::Int64, k::Int64)::Array{Float64}
	[...]
end


#--------------------------------------------------------------------------------------------------------------------------------


"""
	bspline(npts, k, p1, b)

Generate a B-spline curve using an open uniform knot vector.
"""

function bspline(npts::Int64, k::Int64, p1::Int64, b::Array{Float64})::Array{Float64}
	[...]
end


#--------------------------------------------------------------------------------------------------------------------------------


"""
	bsplineu(npts k, p1, b)

Generate a B-spline curve using a periodic uniform knot vector.
"""

function bsplineu(npts::Int64, k::Int64, p1::Int64, b::Array{Float64})::Array{Float64}
	[...]
end


#--------------------------------------------------------------------------------------------------------------------------------


"""
	dbspline(npts k, p1, b)

Generate a B-spline curve and its derivatives using an open uniform knot vector.
"""

function dbspline(npts::Int64, k::Int64, p1::Int64, b::Array{Float64})::tuple{Array{Float64}, Array{Float64}, Array{Float64}}
	[...]
end


#--------------------------------------------------------------------------------------------------------------------------------


"""
	dbsplineu(npts k, p1, b)

Generate a B-spline curve and its derivatives using an open uniform knot vector .
"""

function dbsplineu(npts::Int64, k::Int64, p1::Int64, b::Array{Float64})::tuple{Array{Float64}, Array{Float64}, Array{Float64}}
	[...]
end


#--------------------------------------------------------------------------------------------------------------------------------


"""
	matpbspl(npts k, p1, b)

Generate a B-spline curve using matrix methods and a periodic uniform knot vector.
"""

function matpbspl(npts::Int64, k::Int64, p1::Int64, b::Array{Float64})::Array{Float64}
	[...]
end


#--------------------------------------------------------------------------------------------------------------------------------


"""
	nmatrix(k)

Calculate the general B-spline periodic basis matrix.

This function is used in matpbspl.
"""

function namtrix(k::Int64)::tuple{Float64, Array{Int64}}
    Ni = (n,i) -> factorial(n) / (factorial(i) * factorial(n-i)) 
    n::Array{Int64} = Array{Int64}(zeros(ord,ord))
    fcoeff = 1 / Factoril(ord-1)
    
    for i = 0:ord-1
        temp = Ni[ord-1,i]
        for j = 0:ord-1
            sum = 0
            for k = j:ord-1
                sum1 = ( ord - ( k + 1 ) ) ^ i
                sum2 = (-1) ^ (k - j)
                sum3 = Ni(ord,k - j)
                sum = sum + sum1 * sum2 * sum3
            end
            n[i+1,j+1] = temp * sum
        end
    end
    return(fcoeff,n)            
end



#--------------------------------------------------------------------------------------------------------------------------------


"""
	param(dpts, d)

Calculate parameter values based on chord distances.
"""

function param(dpts::Int64, d::Array{Float64})::Array{Float64}
    sum::Float64 = 0
    isum::Float64 = 0
    tparm::Array{Float64} = Array{Float64}(zeros(dpts))
    #calculate the sum of the chord distances for all the data points
    for i=2:dpts
        sum = sum + sqrt((d[i,1] - d[i-1,1]) ^ 2 + (d[i,2] - d[i-1,2]) ^ 2)
    end
    #calculate the parameter values
    for i=2:dpts
        isum = isum + sqrt((d[i,1] - d[i-1,1]) ^ 2 + (d[i,2] - d[i-1,2]) ^ 2)
        tparm[i] = isum/sum
    end
    return(tparm)
end


