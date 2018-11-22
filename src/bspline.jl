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
	[...]
end


#--------------------------------------------------------------------------------------------------------------------------------


"""
	param(dpts, d)

Calculate parameter values based on chord distances.
"""

function param(dpts::Int64, d::Array{Float64})::Array{Float64}
	[...]
end


