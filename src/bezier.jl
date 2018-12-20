export bezier, dbezier


"""
	bezier(npts, Array{Float64}, cpts)

Calculate a Bezier curve.
"""

function bezier(npts::Int64, b::Array{Float64}, cpts::Int64)::Array{Float64}
	#[...]
	return [0.0 0]
}
end

#-----------------------------------------------------------------------

"""
	dbezier(npts, b, cpts)

Calculate a B~zier curve and its first and second derivatives.
"""

function dbezier(npts::Int64, b::Array{Float64}, cpts::Int64)::tuple{Array{Float64}, Array{Float64}, Array{Float64}}
	#[...]
	return ([0.0 0], [0.0 0], [0.0 0])
end

#-----------------------------------------------------------------------

"""
	bern_basis(n, i, t)

Calculate the Bernstein basis.
"""

function bern_basis(n::Int64, i::Int64, t::Float64)::Float64
	#[...]
	return 0.0
end
