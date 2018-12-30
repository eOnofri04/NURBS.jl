module NURBS
	export bezier, dbezier
	export knot, knotu, knotc
	export basis, dbasis, dbasisu
	export bspline, bsplineu, dbspline, dbsplineu, matpbspl, nmatrix
	#export rbasis, rbspline, rbsplinu
	#export bezsurf, mbezsurf
	export bsplsurf, bspsurfu, dbsurf
	export frbsurf, rbspsurf
	export raise23, raise45

	include("bezier.jl")
	include("knot.jl")
	include("basis.jl")
	include("bspline.jl")
	#include("rbspline.jl")
	#include("beziersurfaces.jl")
	include("bsurfaces.jl")
	include("rbsurfaces.jl")
	include("raise.jl")
end
