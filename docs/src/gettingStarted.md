# Installation

If you have not done so already, download and install [Julia](https://julialang.org/downloads/).

This package was developed with `julia 0.6.4` but it should be supported also on `julia 1.x.x`.


To install NURBS package first you have to download [this](https://github.com/eOnofri04/NURBS.jl) repository.

Then start julia and navigate to the folder called `NURBS.jl` (with the command `;cd`);
here run
```julia
include("src/NURBS.jl")
using NURBS
```

If you want, you can test if everithyng is working fine by running
```julia
include("test/runtests.jl")
```

this will run all the tests written so far and checking all is working properly.

If you want to use a Graphic Interface to preview the curves and the surfaces you are going to build up
you can also use [Plasm](https://github.com/cvdlab/Plasm.jl) package developed by [CVD-LAB](https://github.com/cvdlab) by running
```julia
pkg.add("Plasm")
```

# First Steps

!ToDo