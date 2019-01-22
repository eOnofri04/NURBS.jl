push!(LOAD_PATH,"../src/")

using Documenter, NURBS

makedocs(
	format = :html,
	sitename = "NURBS.jl",
	assets = ["assets/nurbs.css", "assets/logo.png"],
	pages = [
		"Home" => "index.md",
		"Getting Started" => "gettingStarted.md"
		"Curves" => [
			"Bezier Curves" => "bezierC.md",
			"B-Spline Curves" => "splineC.md",
			"Rational B-Spline Curves" => "rSplineC.md"
		],
		"Surfaces" => [
			"Bezier Surfaces" => "bezierS.md",
			"B-Spline Surfaces" => "splineS.md",
			"Rational B-Spline Surfaces" => "rSplineS.md"
		],
		"Authors" => "authors.md",
	]
)
