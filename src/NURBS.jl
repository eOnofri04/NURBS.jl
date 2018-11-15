module NURBS
	#import ...
	export bezier, dbezier
	export knot, knotu, knotc
	export basis, dbasis, dbasisu
	export bspline, bsplineu, dbspline, dbsplineu
	# more to go

	include("bezier.jl")
	include("knot.jl")
	include("basis.jl")
	include("bspline.jl")
	# more to go
end