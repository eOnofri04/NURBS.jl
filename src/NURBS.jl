module NURBS
	export bezier, dbezier
	export knot, knotu, knotc
	export basis, dbasis, dbasisu
	export bspline, bsplineu, dbspline, dbsplineu, matpbspl, nmatrix
	export rbasis, rbspline, rbsplinu
	export bezsurf
	export bsplsurf, bspsurfu, dbsurf, matrixbsplsurf, matrixbsplsurfu
	export rbspsurf, matrixnurbs, constructmatrix
	export raise23, raise45

	export bezsurfj, bezsurf1, bezsurf2, bezsurf3, runn, runner

	include("bezier.jl")
	include("knot.jl")
	include("basis.jl")
	include("bspline.jl")
	include("rbspline.jl")
	include("beziersurfaces.jl")
	include("bsurfaces.jl")
	include("rbsurfaces.jl")
	include("raise.jl")

	include("bezsurftest.jl")
end
