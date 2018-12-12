export bspline, bsplineu, dbspline, dbsplineu, bsplfit

"""
	bsplfit(dpts, d, npts, k)

Fit a B-spline curve using an open uniform knot vector.
"""

function bsplfit(dpts::Int64, d::Array{Float64}, npts::Int64, ord::Int64)::Array{Float64}

	b::Array{Float64} = zeros(dpts, 2) 			# Array containing the control polygon vertices
	nbasis::Array{Float64} = zeros(npts) 		# Array containing the basis ]unctions ]or a single value of t
	n::Array{Float64} = zeros(dpts, npts) 		# Matrix of basis function
	ntrn::Array{Float64} = zeros(npts, dpts) 	# Transpose of the n matrix
	ntemp::Array{Float64} = zeros(npts, npts) 	# Temporary matrix to hold trn(n) x n
	ntmp::Array{Float64} = zeros(npts, 2) 		# Temporary matrix to hold inverse of trn(n) x n x d
	ninv::Array{Float64} = zeros(npts, npts) 	# Inverse of trn(n) x n
	m::Int8 = npts + ord						# Number of knot values

	x::Array{Float64} = knot(npts, ord)
	tpar::Array{Float64} = param(dpts, d)
	

	# Generate the matrix of basis functions
	for i = 1 : dpts

		t = tpar[i] * x[m]
		nbasis = basis(ord, t, npts, x)

		for j = 1 : npts

			n[i,j] = nbasis[j]

		end

	end

	ntrn = n'
	ntemp = ntrn * n
	ninv = inv(ntemp)
	ntmp = ntrn * d
	b = ninv * ntmp

	return b

end

#-----------------------------------------------------------------------

"""
	bspline(npts, k, p1, b)

Generate a B-spline curve using an open uniform knot vector.
"""

function bspline(npts::Int64, k::Int64, p1::Int64, b::Array{Float64})::Array{Float64}
	[...]
end

#-----------------------------------------------------------------------

"""
	bsplineu(npts k, p1, b)

Generate a B-spline curve using a periodic uniform knot vector.
"""

function bsplineu(npts::Int64, k::Int64, p1::Int64, b::Array{Float64})::Array{Float64}
	[...]
end

#-----------------------------------------------------------------------

"""
	dbspline(npts k, p1, b)

Generate a B-spline curve and its derivatives using an open uniform knot vector.
"""

function dbspline(npts::Int64, k::Int64, p1::Int64, b::Array{Float64})::tuple{Array{Float64}, Array{Float64}, Array{Float64}}
	[...]
end

#-----------------------------------------------------------------------

"""
	dbsplineu(npts k, p1, b)

Generate a B-spline curve and its derivatives using an open uniform knot vector .
"""

function dbsplineu(npts::Int64, k::Int64, p1::Int64, b::Array{Float64})::tuple{Array{Float64}, Array{Float64}, Array{Float64}}
	[...]
end

#-----------------------------------------------------------------------

"""
	matpbspl(npts k, p1, b)

Generate a B-spline curve using matrix methods and a periodic uniform knot vector.
"""

function matpbspl(npts::Int64, k::Int64, p1::Int64, b::Array{Float64})::Array{Float64}
	[...]
end

#-----------------------------------------------------------------------

"""
	nmatrix(k)

Calculate the general B-spline periodic basis matrix.

This function is used in matpbspl.
"""

function namtrix(k::Int64)::tuple{Float64, Array{Int64}}
	[...]
end

#-----------------------------------------------------------------------

"""
	param(dpts, d)

Calculate parameter values based on chord distances.
"""

function param(dpts::Int64, d::Array{Float64})::Array{Float64}

	tparm::Array{Float64} = zeros(dpts)

    # calculate the chord distances for all the data points
    arc_lengths::Array{Float64} = [sqrt((d[i,1] - d[i-1,1])^2 + (d[i,2] - d[i-1,2])^2) for i = 2 : dpts]

    tparm[1] = 0

    for i = 2 : dpts
        tparm[i] = tparm[i-1] + arc_lengths[i-1] 
    end

    return tparm / sum(arc_lengths)

end


