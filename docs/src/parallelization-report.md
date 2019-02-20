# Parallelization Tests

## Bezier Surfaces

#### Results over 20 runs with:
 - `npts = 50`
 - `mpts = 50`
 - `udpts = 2500`
 - `wdpts = 2500`


  **1 proc**  |   **min**  |   **avg**  |   **max**
--------------|------------|------------|------------
   *serial*   | `3.817236` | `4.302450` | `4.587352`
 *parallel 1* | `3.803109` | `4.229924` | `4.662597`
 *parallel 2* | `3.813949` | `4.231862` | `4.651251`
 *parallel 3* | `3.813311` | `4.274616` | `4.871468`

  **10 proc** |   **min**  |   **avg**  |   **max**
--------------|------------|------------|------------
   *serial*   | `3.839719` | `4.609300` | `6.655391`
 *parallel 1* | `3.826293` | `4.541800` | `5.381051`
 *parallel 2* | `3.843635` | `4.517371` | `5.319194`
 *parallel 3* | `3.834900` | `4.452138` | `5.040315`

the test made are:
 - *serial*: no parallelization esplicitly made.
 - *parallel 1*: with `@sync` and `@async` and matrices defined inside `@sync` block.
 - *parallel 2*: with `@sync` and `@async` and matrices defined before `@sync` block.

#### Result over 50 runs with:
 - `npts = 50`
 - `mpts = 50`
 - `udpts = 2500`
 - `wdpts = 2500`

  **1 proc**  |  **average** 
--------------|--------------
   *serial*   | `4.44101856`
 *parallel 1* | `4.45837448`
 *parallel 2* | `4.47456267`
 *parallel 3* | `4.43402378`
--------------|--------------


 **10 proc**  |  **average** 
--------------|--------------
   *serial*   | `4.61977279`
 *parallel 1* | `4.60410667`
 *parallel 2* | `4.55359172`
 *parallel 3* | `4.51270721`

### The Code

The initialisation of the calculus:
```julia
begin
    npts = 50
    mpts = 50
    udpts = 2500
    wdpts = 2500
    bplus = zeros(3, npts * mpts)

    for i = 1 : npts
        ki = i/2
        for j = 1 : mpts
            kj = j /2
            r = sqrt((ki)^2+(kj)^2)
            if (r == 0)
                r = 0.000001
            end
            ij = (i-1)*npts + j
            bplus[1,ij] = i
            bplus[2,ij] = j# + i/10
            bplus[3,ij] = i+j
        end
    end
end
```

and the execution:
```julia
begin
	minor = 10.0
	medium = 0.0 
	major = 0.0
	for i = 1 : 20
	    x = @elapsed bezsurf(npts, mpts, bplus, 1000, 1000)
	    if (x < minor)
	        minor = x
	    end
	    if (x > major)
	        major = x
	    end
	    medium += x
	end
	medium /= 20
	println("minimum = $(minor)")
	println("average = $(medium)")
	println("maximum = $(major)")
end
```

## Conclusions

New tests have been made on small and medium data sets.

All of them have enhanced the same results: over "not so big" data sets no particular changes have been discovered.

At this state of art Julia parallelism seems to be the best choice as manage memory allocation of Activation Registers in the most efficient way.
Usually Julia's default parallelization seems to be a little slowler. However this virtual time cost is totally overcome by the smaller amount of memory leack that, in user parallelizations, come at the cost of a high number of Garbage Collector calls.

More test should be made in order to have a precise results sets but it seems useless at the current state of art to waste time over its analysis.