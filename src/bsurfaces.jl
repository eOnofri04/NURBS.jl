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

function bsplsurf(B::Array{Float64,2}, ordx::Int64, ordy::Int64, npts::Int64, mpts::Int64, p1::Int64, p2::Int64)::Tuple{Array{Float64, 2}, Array{Array{Int64, 1}, 1}, Array{Array{Int64, 1}, 1}}
    
    nplusc = npts + ordx
    mplusc = mpts + ordy
    x = zeros(nplusc)
    y = zeros(mplusc)
    nbasis = zeros(1, npts)
    mbasis = zeros(1, mpts)
    q = zeros(3, p1*p2)
    
    #generate the open uniform knot vectors
    x = knot(npts, ordx)
    y = knot(mpts, ordy)
   
    icount = 0
    #calculate the points on the B-spline surface
    stepu = x[nplusc] / (p1 - 1)
    stepw = y[mplusc] / (p2 - 1)
    for u = 0 : stepu : x[nplusc]
        nbasis = basis(npts, ordx, u, x)
        for w = 0 : stepw : y[mplusc]
            mbasis = basis(mpts, ordy, w, y)
            icount = icount + 1
            for i = 1 : npts
                for j = 1 : mpts
                    j1 = mpts * (i - 1) + j
                    for s = 1 : 3
                        q[s, icount] = q[s, icount] + B[s, j1] * nbasis[i] * mbasis[j]
                    end
                end
            end
        end
    end

    EV = Array{Int64,1}[] #1-dimensional cellular complex used for plotting the curve with Plasm package

    for i = 1 : p1 * p2 - 1
        C = Array{Int64}(2)
        C[1] = i
        C[2] = i+1
        push!(EV,C) 
    end
    
    for i = 1 : p2 * (p1 - 1)
        C = Array{Int64}(2)
        C[1] = i
        C[2] = i+p2
        push!(EV,C) 
    end

    FV = Array{Int64,1}[] #1-dimensional cellular complex used for plotting the curve with Plasm package

    for j = 0 : p1 - 2
        for i = 1 : p2 - 1
            C = Array{Int64}(4)
            C[1] = j * p2 + i
            C[2] = j * p2 + i + 1
            C[3] = (j + 1) * p2 + i
            C[4] = (j + 1) * p2 + i + 1
            push!(FV,C) 
        end
    end    
    return(q, EV, FV)
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

function bsplsurfu(B::Array{Float64,2}, ordx::Int64, ordy::Int64, npts::Int64, mpts::Int64, p1::Int64, p2::Int64)::Tuple{Array{Float64, 2}, Array{Array{Int64, 1}, 1}, Array{Array{Int64, 1}, 1}}
    
    nplusc = npts + ordx
    mplusc = mpts + ordy
    x = zeros(nplusc)
    y = zeros(mplusc)
    nbasis = zeros(1, npts)
    mbasis = zeros(1, mpts)
    q = zeros(3, p1 * p2)
    
    #generate the open uniform knot vectors
    x = knotu(npts, ordx)
    y = knotu(mpts, ordy)
   
    icount = 0
    #calculate the points on the B-spline surface
    stepu = (npts - ordx + 1) / (p1 - 1)
    stepw = (mpts - ordy + 1) / (p2 - 1)
    for u = ordx-1 : stepu : npts
        nbasis = basis(npts, ordx, u, x)
        for w = ordy-1 : stepw : mpts
            mbasis = basis(mpts, ordy, w, y)
            icount = icount + 1
            for i = 1 : npts
                for j = 1 : mpts
                    for s = 1 : 3
                        j1 = mpts * (i - 1) + j
                        q[s, icount] = q[s, icount] + B[s, j1] * nbasis[i] * mbasis[j]
                    end
                end
            end
        end
    end

    EV = Array{Int64,1}[] #1-dimensional cellular complex used for plotting the curve with Plasm package

    for i = 1 : p1 * p2 - 1
        C = Array{Int64}(2)
        C[1] = i
        C[2] = i+1
        push!(EV,C) 
    end
    
    for i = 1 : p2 * (p1 - 1)
        C = Array{Int64}(2)
        C[1] = i
        C[2] = i+p2
        push!(EV,C) 
    end

    FV = Array{Int64,1}[] #1-dimensional cellular complex used for plotting the curve with Plasm package

    for j = 0 : p1 - 2
        for i = 1 : p2 - 1
            C = Array{Int64}(4)
            C[1] = j * p2 + i
            C[2] = j * p2 + i + 1
            C[3] = (j + 1) * p2 + i
            C[4] = (j + 1) * p2 + i + 1
            push!(FV,C) 
        end
    end    
    return(q, EV, FV)
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

function dbsurf(B::Array{Float64,2}, ordx::Int64, ordy::Int64, npts::Int64, mpts::Int64, p1::Int64, p2::Int64)::Tuple{Array{Float64,2}, Array{Float64,2}, Array{Float64,2}, Array{Float64,2}, Array{Float64,2}, Array{Float64,2}, Array{Array{Int64, 1}, 1}, Array{Array{Int64, 1}, 1}}
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
    q = zeros(3, p1 * p2)
    qu = zeros(3, p1 * p2)
    qw = zeros(3, p1 * p2)
    quu = zeros(3, p1 * p2)
    quw = zeros(3, p1 * p2)
    qww = zeros(3, p1 * p2)
    
    x = knot(npts, ordx)
    y = knot(mpts, ordy)
    
    icount = 0
    #calculate the points on the B-spline surface
    stepu = x[nplusc] / (p1 - 1)
    stepw = y[mplusc] / (p2 - 1)
    for u = 0 : stepu : x[nplusc]
        nbasis, d1nbasis, d2nbasis = dbasis(npts, ordx, u, x)
        for w = 0 : stepw : y[mplusc]
            mbasis, d1mbasis, d2mbasis = dbasis(mpts, ordy, w, y)
            icount = icount + 1
            for i = 1 : npts
                for j = 1 : mpts                    
                    for s = 1 : 3
                        j1 = mpts * (i - 1) + j
                        q[s, icount] = q[s, icount] + B[s, j1] * nbasis[i] * mbasis[j]
                        qu[s, icount] = qu[s, icount] + B[s, j1] * d1nbasis[i] * mbasis[j]
                        qw[s, icount] = qw[s, icount] + B[s, j1] * nbasis[i] * d1mbasis[j]
                        quu[s, icount] = quu[s, icount] + B[s, j1] * d1nbasis[i] * d1mbasis[j]
                        quw[s, icount] = quw[s, icount] + B[s, j1] * d2nbasis[i] * mbasis[j]
                        qww[s, icount] = qww[s, icount] + B[s, j1] * nbasis[i] * d2mbasis[j]
                    end
                end
            end
        end
    end

    EV = Array{Int64,1}[] #1-dimensional cellular complex used for plotting the curve with Plasm package

    for i = 1 : p1 * p2 - 1
        C = Array{Int64}(2)
        C[1] = i
        C[2] = i+1
        push!(EV,C) 
    end
    
    for i = 1 : p2 * (p1 - 1)
        C = Array{Int64}(2)
        C[1] = i
        C[2] = i+p2
        push!(EV,C) 
    end

    FV = Array{Int64,1}[] #1-dimensional cellular complex used for plotting the curve with Plasm package

    for j = 0 : p1 - 2
        for i = 1 : p2 - 1
            C = Array{Int64}(4)
            C[1] = j * p2 + i
            C[2] = j * p2 + i + 1
            C[3] = (j + 1) * p2 + i
            C[4] = (j + 1) * p2 + i + 1
            push!(FV,C) 
        end
    end    

    return(q, qu, qw, quu, quw, qww, EV, FV)

end
