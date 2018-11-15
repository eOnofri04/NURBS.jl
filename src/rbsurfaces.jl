export frbsurf, rbspsurf


"""
	frbsrf(b, k, l, npts, mpts, p1, p2, p_itest, ibnum, bold, niku, mjlw, rsumij, savrsumij)

Calculates and test the fast B-spline surface algorithm.

Call: basis, knot, sumrbas
"""

function frbsurf(b::Array{Float64}, k::Int64, l::Int64, npts::Int64, mpts::Int64, p1::Int64, p2::Int64, p_itest, ibnum::Int64, bold::Array{Float64}, niku::Array{Float64}, mjlw::Array{Float64}, rsumij::Array{Float64}, savrsumij::Array{Float64})::Array{Float64}
	[...]
end

#-----------------------------------------------------------------------

"""
	rbspsurf(b::Array{Float64}, k::Int64, l::Int64, npts::Int64, mpts::Int64, p1::Int64, p2::Int64)

Calculate a Cartesian product rational B-spline surface (NURBS) using an open uniform knot vector.

Call: basis, knot, sumrbas
"""

function rbspsurf(b, k, l, npts, mpts, p1, p2)::Array{Float64}
	[...]
end

#-----------------------------------------------------------------------

"""
	sumrbas(b, nbasis, mbasis, npts, mpts)

Calculate the sum of the nonrational basis functions.
"""

function sumrbas(b::Array{Float64}, nbasis::Array{Float64}, mbasis::Array{Float64}, npts::Int64, mpts::Int64)::Float64
	[...]
end