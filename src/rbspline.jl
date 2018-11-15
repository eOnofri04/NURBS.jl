export rbasis, rbspline, rbsplinu


"""
	rbasis(c, t, npts, x, h)

Generate a rational B-spline basis function using an open knot vector.
"""

function rbasis(c::Int64, t::Float64, npts::Int64, x::Array{Int64}, h::Array{Float64})::Array{Float64}
	[...]
end

#-----------------------------------------------------------------------

"""
	rbspline(npts, k, p1, b, h)

Generate a rational B-spline curve using an open uniform knot vector.
"""

function rbspline(npts::Int64, k::Int64, p1::Int64, b::Array{Float64}, h::Array{Float64})::Array{Float64}
	[...]
end

#-----------------------------------------------------------------------

"""
	rbsplinu(npts, k, p1, b, h)

Generate a rational B-spline curve using a periodic uniform knot vector.
"""

function rbsplinu(npts::Int64, k::Int64, p1::Int64, b::Array{Float64}, h::Array{Float64})::Array{Float64}
	[...]
end

