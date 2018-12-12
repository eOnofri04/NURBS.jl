export bsplsurf, bspsurfu, dbsurf


"""
UNTESTED	bsplsurf(B[], ordx, ordy, npts, mpts, p1, p2)

Calculate a Cartesian product B-spline surface using open uniform knot vectors.

---

# Arguments
- `B::Array{Float64}`: array containing the control net vertices
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

function bsplsurf(B::Array{Float64}, ordx::Int64, ordy::Int64, npts::Int64, mpts::Int64, p1::Int64, p2::Int64)::Array{Float64}
    
    nplusc = npts + ordx
    mplusc = mpts + ordy
    x = zeros(nplusc)
    y = zeros(mplusc)
    nbasis = zeros(1, npts)
    mbasis = zeros(1, mpts)
    q = zeros(p1*p2, 3)
    
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
            icount = icount + 1
            for i = 1 : npts
                for j = 1 : mpts
                    j1 = mpts * (i - 1) + j
                    for s = 1 : 3
                        q[icount, s] = q[icount, s] + B[j1, s] * nbasis[i] * mbasis[j]
                    end
                end
            end
        end
    end
    return(q)
end

#-----------------------------------------------------------------------

"""
UNTESTED	bspsurfu(B[], ordx, ordy, npts, mpts, p1, p2)

Calculate a Cartesian product B-spline surface using periodic uniform knot vectors.

---

# Arguments
- `B::Array{Float64}`: array containing the control net vertices
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

function bsplsurfu(B::Array{Float64}, ordx::Int64, ordy::Int64, npts::Int64, mpts::Int64, p1::Int64, p2::Int64)::Array{Float64}
    
    nplusc = npts + ordx
    mplusc = mpts + ordy
    x = zeros(nplusc)
    y = zeros(mplusc)
    nbasis = zeros(1, npts)
    mbasis = zeros(1, mpts)
    q = zeros(p1 * p2, 3)
    
    #generate the open uniform knot vectors
    x = knotu(npts, ordx)
    y = knotu(mpts, ordy)
   
    icount = 0
    #calculate the points on the B-spline surface
    stepu = (npts - k + 1) / (p1 - 1)
    stepw = (mpts - k + 1) / (p2 - 1)
    for u = ordx-1 : stepu : npts
        nbasis = basis(ordx, u, npts, x)
        for w = ordy-1 : stepw : mpts
            mbasis = basis(ordy, w, mpts, y)
            icount = icount + 1
            for i = 1 : npts
                for j = 1 : mpts
                    for s = 1 : 3
                        j1 = mpts * (i - 1) + j
                        q[icount, s] = q[icount, s] + B[j1, s] * nbasis[i] * mbasis[j]
                    end
                end
            end
        end
    end
    return(q)
end

#-----------------------------------------------------------------------

"""
UNTESTED	dbsurf(B[], ordx, ordy, npts, mpts, p1, p2)

Calculate a Cartesian product B-spline surface and its derivatives using open uniform knot vectors.

---

# Arguments
- `B::Array{Float64}`: array containing the control net vertices
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

function dbsurf(B::Array{Float64}, ordx::Int64, ordy::Int64, npts::Int64, mpts::Int64, p1::Int64, p2::Int64)::Tuple{Array{Float64}, Array{Float64}, Array{float64}, Array{Float64}, Array{Float64}, Array{float64}}
    nplusc = npts + ordx
    mplusc = mpts + ordy
    x = zeros(nplusc)
    y = zeros(mplusc)
    nbasis = zeros(1, npts)
    mbasis = zeros(1, mpts)
    d1nbasis = zeros(1, npts)
    d1mbasis = zeros(1, mpts)
    d2nbasis = zeros(1, npts)
    d2mbasis = zeros(1, mpts)
    q = zeros(p1 * p2, 3)
    qu = zeros(p1 * p2, 3)
    qw = zeros(p1 * p2, 3)
    quu = zeros(p1 * p2, 3)
    quw = zeros(p1 * p2, 3)
    qww = zeros(p1 * p2, 3)
    
    x = knot(npts, ordx)
    y = knot(mpts, ordy)
    
    icount = 0
    #calculate the points on the B-spline surface
    stepu = x[nplusc] / (p1 - 1)
    stepw = y[mplusc] / (p2 - 1)
    for u = 0 : stepu : x[nplusc]
        nbasis = dbasis(k, u, npts, x, d1nbasis, d2nbasis)
        for w = 0 : stepw : y[mplusc]
            mbasis = dbasis(l, w, mpts, y, d1mbasis, d2mbasis)
            icount = icount + 1
            for i = 1 : npts
                for j = 1 : mpts                    
                    for s = 1 : 3
                        j1 = mpts * (i - 1) + j
                        q[icount, s] = q[icount, s] + B[j1, s] * nbasis[i] * mbasis[j]
                        qu[icount, s] = qu[icount, s] + B[j1, s] * d1nbasis[i] * mbasis[j]
                        qw[icount, s] = qw[icount, s] + B[j1, s] * nbasis[i] * d1mbasis[j]
                        quu[icount, s] = quu[icount, s] + B[j1, s] * d1nbasis[i] * d1mbasis[j]
                        quw[icount, s] = quw[icount, s] + B[j1, s] * d2nbasis[i] * mbasis[j]
                        qww[icount, s] = qww[icount, s] + B[j1, s] * nbasis[i] * d2mbasis[j]
                    end
                end
            end
        end
    end
    return(q, qu, qw, quu, quw, qww)

end
