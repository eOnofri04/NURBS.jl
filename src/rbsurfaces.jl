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
UNTESTED	rbspsurf(b[], ordx, ordy, npts, mpts, p1, p2)

Calculate a Cartesian product rational B-spline surface (NURBS) using an open uniform knot vector.

Call: basis, knot, sumrbas

---

# Arguments
- `B::Array{Float64}`: .
- `ordx::Int64`: order in the x direction
- `ordy::Int64`: order in the y direction
- `npts::Int64`: the number of control net vertices in x direction
- `mpts::Int64`: the number of control net vertices in y direction
- `p1::Int64`: number of parametric lines in the x direction
- `p2::Int64`: number of parametric lines in the y direction

---

# Examples
```jldoctest
julia> 
```

```jldoctest
julia> 
```

```jldoctest
julia> 
```

---

_By Paolo Macciacchera, Elia Onofri_
"""

function rbspsurf(B::Array{Float64,2}, ordx::Int64, ordy::Int64, npts::Int64, mpts::Int64, p1::Int64, p2::Int64)::Array{Float64,2}
    
    nplusc = npts + ordx
    mplusc = mpts + ordy
    x = zeros(nplusc)
    y = zeros(mplusc)
    nbasis = zeros(1, npts)
    mbasis = zeros(1, mpts)
    q = zeros(3, p1 * p2)
    
    #generate the open uniform knot vectors
    x = knot(npts, ordx)
    y = knot(mpts, ordy)
   
    icount = 0
    #calculate the points on the B-spline surface
    stepu = x[nplusc] / (p1 - 1)
    stepw = y[mplusc] / (p2 - 1)
    for u = 0 : stepu : x[nplusc]
        nbasis = basis(ordx, u, npts, x)
        for w = 0 : stepw : y[mplusc]
            mbasis = basis(ordy, w, mpts, y)
            sum = sumrbas(B, nbasis, mbasis, npts, mpts)
            icount = icount + 1
            for i = 1 : npts
                for j = 1 : mpts
                    for s = 1 : 3
                        j1 = mpts * (i - 1) + j
                        qtemp = (B[4, j1] * B[s, j1] * nbasis[i] * mbasis[j]) / sum
                        q[s, icount] = q[s, icount] + qtemp
                    end
                end
            end
        end
    end
    return(q)
end

#-----------------------------------------------------------------------

"""
UNTESTED	sumrbas(B[], nbasis, mbasis, npts, mpts)

Calculate the sum of the nonrational basis functions.

---

# Arguments
- `B::Array{Float64}`: array containing the control net vertices
- `nbasis::Array{Float64}`: array containing the nonrational basis functions for x
- `mbasis::Array{Float64}`: array containing the nonrational basis functions for y
- `npts::Int64`: the number of control net vertices in x direction
- `mpts::Int64`: the number of control net vertices in y direction

---

# Examples
```jldoctest
julia> 
```

```jldoctest
julia> 
```

```jldoctest
julia> 
```

---

_By Paolo Macciacchera, Elia Onofri_
"""

function sumrbas(B::Array{Float64,2}, nbasis::Array{Float64}, mbasis::Array{Float64}, npts::Int64, mpts::Int64)::Float64
    sum = 0
    for i = 1 : npts
        for j = 1 : mpts
            j1 = mpts * (i - 1) + j
            sum = sum + B[4, j1] * nbasis[i] * mbasis[j]
        end
    end
    return(sum)
end
