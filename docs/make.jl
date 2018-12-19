push!(LOAD_PATH,"../src/")

using Documenter, NURBS

makedocs(
	format = :html,
	sitename = "NonUniformRationalBSplines.jl",
	assets = ["assets/nurbs.css", "assets/logo.png"],
	pages = [
		"Home" => "index.md",
		"L.A.R. Intro" => "lar.md",
		"Interface" => "interface.md",
		"Arrangement" => "arrangement.md",
		"Parametric primitives" => [
			"Mapper" => "mapper.md",
			"Assemblies" => "struct.md"
		],
		"Grid generation" => [
			"Cuboidal grids" => "largrid.md",
			"Simplicial grids" => "simplexn.md"
		],
		"Domain integration" => "integr.md",
	]
)
