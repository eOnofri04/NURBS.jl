export basis, dbasis, dbasisu

"""
    basis(ord, t, npts, x[])

Generate a B-spline basis functions `N[]` for a given open knot vectors `x[]`.

A B-Spline basis is a collection of functions of a parameter `t` wich form
 a basis for the vectorial space of functions. The conformation of this set
 higly depends on the choosing of the knots `x[]` the curve is bound to.
The basis is computed with _Cox-de Boor_ recursive function applied to the
 basis dependency tree in order to optimise the computation.

---

# Arguments
- `ord::Int64`: the order of the B-Spline (`deg = ord-1`).
- `t::Float64`: the parameter value of the parametric curve.
- `npts::Int64`: the number of the points of the controll polygon.
- `x::Array{Float64}`: the knot vector.

---

_By Elia Onofri_
"""

function basis(ord::Int64, t::Float64, npts::Int64, x::Array{Float64})::Array{Float64}
    local tmp = Float64[]       # Basis progressive vector
    local N = Float64[]         # Output vector
    local max_N = npts-1+ord    # Needs i+(ord-1)|i=npts = npts-1+ord trivial basis
    local ddep::Float64 = 0.0   # Direct Dependency partial sum
    local fdep::Float64 = 0.0   # Forward Dependency partial sum
    
    # Local check of the knot vector correctness
    @assert length(x) == npts+ord ("ERROR: incompatibile knot vector with given parameters n+1 = $(npts), k = $(ord)")
    
    # Eval N_{i,1} for i = 1:max_B
    for i = 1:max_N
        if (t>=x[i]) && (t<x[i+1])
            append!(tmp, 1)
        else
            append!(tmp, 0)
        end
    end
    
    # Eval higher basis N_{i,deg} for deg = 2:ord and i = 1:max_B-deg
    for deg = 2:ord
        for i = 1:max_N+1-deg
            # Eval of the direct dependency
            if tmp[i]==0
                ddep = 0.0
            else
                ddep = ((t-x[i])*tmp[i])/(x[i+deg-1]-x[i])
            end
            # Eval of the forward dependency
            if tmp[i+1]==0
                fdep = 0.0
            else
                fdep = ((x[i+deg]-t)*tmp[i+1])/(x[i+deg]-x[i+1])
            end
            # Collection of the dependencies
            tmp[i] = ddep+fdep
        end
        temp[i]=d+e
    end
end

if t==x[nplusc]  #pick up last point
    temp[npts]=1
end

    
    # Otherwise last point is zero
    if t == x[npts+ord]
        tmp[npts] = 1
    end
    
    # Collect N{1,ord} to N{npts,ord} in B
    for i=1:npts
        push!(N, tmp[i]);
    end
    return N;

end


#--------------------------------------------------------------------------------------------------------------------------------


"""
	dbasis(c, t, npts, x)

Generate B-spline basis functions and their derivatives for uniform open knot vectors.

---

_By Elia Onofri_
"""

function dbasis(c::Int64,t::Float64,npts::Int64,x::Array{Float64},n::Array{Float64},d1::Array{Float64},d2::Array{Float64})
    nplusc=npts+c
    #zero the temporary arrays
    temp=zeros(nplusc)
    temp1=zeros(nplusc)
    temp2=zeros(nplusc)
    #calculate the first-order basis functions n(i,1)
    for i=1:nplusc-1
        if t>=x[i] && t<x[i+1]
            temp[i]=1
        else
            temp[i]=0
        end
    end
    
    if t==x[nplusc]  temp[npts]=1 end
    
    #calculate higher-order basis functions and their derivatives
    
    for k=2:c
        for i=1:nplusc-k
            
            #calculate basis function
            if temp[i]!=0
                b1=((t-x[i])*temp[i])/(x[i+k-1]-x[i])
            else
                b1=0
            end
            if temp[i+1]!=0
                b2=((x[i+k]-t)*temp[i+1])/(x[i+k]-x[i+1])
            else
                b2=0
            end
            
            #calculate first derivative
            if temp[i]!=0
                f1=temp[i]/(x[i+k-1]-x[i])
            else
                f1=0
            end
            if temp[i+1]!=0
                f2=-temp[i+1]/(x[i+k]-x[i+1])
            else
                f2=0
            end
            if temp1[i]!=0
                f3=((t-x[i])*temp1[i])/(x[i+k-1]-x[i])
            else
                f3=0
            end
            if temp1[i+1]!=0
                f4=((x[i+k]-t)*temp1[i+1])/(x[i+k]-x[i+1])
            else
                f4=0
            end
            
            #calculate second derivative
            if temp1[i]!=0
                s1=(2*temp1[i])/(x[i+k-1]-x[i])
            else
                s1=0
            end
            if temp1[i+1]!=0
                s2=(-2*temp1[i+1])/(x[i+k]-x[i+1])
            else
                s2=0
            end
            if temp2[i]!=0
                s3=((t-x[i])*temp2[i])/(x[i+k-1]-x[i])
            else
                s3=0
            end
            if temp2[i+1]!=0
                s4=((x[i+k]-t)*temp2[i+1])/(x[i+k]-x[i+1])
            else
                s4=0
            end
            
            temp[i]=b1+b2
            temp1[i]=f1+f2+f3+f4
            temp2[i]=s2+s2+s3+s4
        end
    end
    
    #put in arrays
    for i=1:npts
        n[i]=temp[i]
        d1[i]=temp1[i]
        d2[i]=temp2[i]
    end
            
    
end

#--------------------------------------------------------------------------------------------------------------------------------


"""
	dbasisu(c, t, npts, x)

Generate B-spline basis functions and their derivatives for uniform periodic knot vectors.

---

_By Elia Onofri_
"""

function dbasisu(c::Int64,t::Float64,npts::Int64,x::Array{Float64},n::Array{Float64},d1::Array{Float64},d2::Array{Float64})
    nplusc=npts+c
    #zero the temporary arrays
    temp=zeros(nplusc)
    temp1=zeros(nplusc)
    temp2=zeros(nplusc)
    #calculate the first-order basis functions n(i,1)
    for i=1:nplusc-1
        if t>=x[i] && t<x[i+1]
            temp[i]=1
        else
            temp[i]=0
        end
    end
    
    if t==x[npts+1]
        temp[npts]=1
        temp[npts+1]=0
    end
    
    #calculate higher-order basis functions and their derivatives
    
    for k=2:c
        for i=1:nplusc-k
            
            #calculate basis function
            if temp[i]!=0
                b1=((t-x[i])*temp[i])/(x[i+k-1]-x[i])
            else
                b1=0
            end
            if temp[i+1]!=0
                b2=((x[i+k]-t)*temp[i+1])/(x[i+k]-x[i+1])
            else
                b2=0
            end
            
            #calculate first derivative
            if temp[i]!=0
                f1=temp[i]/(x[i+k-1]-x[i])
            else
                f1=0
            end
            if temp[i+1]!=0
                f2=-temp[i+1]/(x[i+k]-x[i+1])
            else
                f2=0
            end
            if temp1[i]!=0
                f3=((t-x[i])*temp1[i])/(x[i+k-1]-x[i])
            else
                f3=0
            end
            if temp1[i+1]!=0
                f4=((x[i+k]-t)*temp1[i+1])/(x[i+k]-x[i+1])
            else
                f4=0
            end
            
            #calculate second derivative
            if temp1[i]!=0
                s1=(2*temp1[i])/(x[i+k-1]-x[i])
            else
                s1=0
            end
            if temp1[i+1]!=0
                s2=(-2*temp1[i+1])/(x[i+k]-x[i+1])
            else
                s2=0
            end
            if temp2[i]!=0
                s3=((t-x[i])*temp2[i])/(x[i+k-1]-x[i])
            else
                s3=0
            end
            if temp2[i+1]!=0
                s4=((x[i+k]-t)*temp2[i+1])/(x[i+k]-x[i+1])
            else
                s4=0
            end
            
            temp[i]=b1+b2
            temp1[i]=f1+f2+f3+f4
            temp2[i]=s2+s2+s3+s4
        end
    end
    
    #put in arrays
    for i=1:npts
        n[i]=temp[i]
        d1[i]=temp1[i]
        d2[i]=temp2[i]
    end
    return (n, d1, d2);
end
