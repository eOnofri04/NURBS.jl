export basis, dbasis, dbasisu


"""
	basis(c, t, npts, x)

Generate B-spline basis functions for open uniform knot vectors.
"""

function basis(c::Int64, t::Float64, npts::Int64, x::Array{Int64})::Array{Float64}
    tempt = Array{Float64}(36)
    n = Float64[];
    nplusc = npts + c;
    
#    print(knot vector is );
#    for i=1:nplusc
#        print(& &, i, x[i]);
#    end
#    print(t is &, t);
    
    # calculate the first order basis functions n[i][1]
    
    for i=1:nplusc-1
        if (t>=x[i])&&(t<x[i+1])
            temp[i] = 1;
        else
            temp[i] = 0;
        end
    end
    
    # calculate the higher order basis functions
    for k=2:c
        for i=1:nplusc-k
            if temp[i] != 0    # if the lower order basis function is zero skip the calculation
                d = ((t-x[i])*temp[i])/(x[i+k-1]-x[i]);
            else
                d = 0;
            end
            if temp[i+1] != 0    # if the lower order basis function is zero skip the calculation
                e = ((x[i+k]-t)*temp[i+1])/(x[i+k]-x[i+1]);
            else
                e = 0;
            end
            temp[i] = d+e;
        end
    end
    
    if t == x[nplusc]    # pick the last point
        temp[npts] = 1;
    end
    
    # put in n array
    for i=1:npts
        push!(n, temp[i]);
    end
    return n;
end


"""
	dbasis(c, t, npts, x)

Generate B-spline basis functions and their derivatives for uniform open knot vectors.
"""

function dbasis(c::Int64, t::Float64, npts::Int64, x::Array{Int64})::Tuple{Array{Float64}, Array{Float64}, Array{Float64}}
    
    #inizialization
    
    temp = zeros(36);    # allows for 35 defining polygon vertices
    temp1 = zeros(36);
    temp2 = zeros(36);
    nplusc = npts+c;
    
    # calculate the first order basis functions n[i]
    
    for i=1:nplusc-1
        if (t>=x[i])&&(t<x[i+1])
            temp[i]=1;
        else
            temp[i]=0;
        end
    end
    
    if t==x[nplusc]    # last(x)
        temp[npts] = 1;
    end
    
    # calculate the higher order basis functions
    for k=2:c
        for i=1:nplusc-k
            if temp[i] != 0    # if the lower order basis function is zero skip the calculation
                b1 = ((t-x[i])*temp[i])/(x[i+k-1]-x[i]);
            else
                b1 = 0;
            end
            if temp[i+1] != 0    # if the lower order basis function is zero skip the calculation
                b2 = ((x[i+k]-t)*temp[i+1])/(x[i+k]-x[i+1]);
            else
                b2 = 0;
            end
            
            # calculate first derivative
            if temp[i] != 0    # if the lower order basis function is zero skip the calculation
                f1 = temp[i]/(x[i+k-1]-x[i]);
            else
                f1 = 0;
            end
            if temp[i+1] != 0    # if the lower order basis function is zero skip the calculation
                f2 = -temp[i+1]/(x[i+k]-x[i+1]);
            else
                f2 = 0;
            end
            if temp1[i] != 0    # if the lower order basis function is zero skip the calculation
                f3 = ((t-x[i])*temp1[i])/(x[i+k-1]-x[i]);
            else
                f3 = 0;
            end
            if temp1[i+1] != 0    # if the lower order basis function is zero skip the calculation
                f4 = ((x[i+k]-t)*temp1[i+1])/(x[i+k]-x[i+1]);
            else
                f4 = 0;
            end
            
            # calculate second derivative
            if temp1[i] != 0    # if the lower order basis function is zero skip the calculation
                s1 = (2*temp1[i])/(x[i+k-1]-x[i]);
            else
                s1 = 0;
            end
            if temp1[i+1] != 0    # if the lower order basis function is zero skip the calculation
                s2 = (-2*temp1[i+1])/(x[i+k]-x[i+1]);
            else
                s2 = 0;
            end
            if temp2[i] != 0    # if the lower order basis function is zero skip the calculation
                s3 = ((t-x[i])*temp2[i])/(x[i+k-1]-x[i]);
            else
                s3 = 0;
            end
            if temp2[i+1] != 0    # if the lower order basis function is zero skip the calculation
                s4 = ((x[i+k]-t)*temp2[i+1])/(x[i+k]-x[i+1]);
            else
                s4 = 0;
            end
            
            temp[i] = b1 + b2;
            temp1[i] = f1 + f2 + f3 + f4;
            temp2[i] = s1 + s2 + s3 + s4;
        end
    end
    
    # prepare output
    for i=1:npts
        push!(n, temp[i]);
        push!(d1, temp1[i]);
        push!(d2, temp2[i]);
    end
    return (n, d1, d2);
end


"""
	dbasisu(c, t, npts, x)

Generate B-spline basis functions and their derivatives for uniform periodic knot vectors.
"""

function dbasisu(c::Int64, t::Float64, npts::Int64, x::Array{Int64})::tuple{Array{Float64}, Array{Float64}, Array{Float64}}

    #inizialization
    
    temp = zeros(36);    # allows for 35 defining polygon vertices
    temp1 = zeros(36);
    temp2 = zeros(36);
    nplusc = npts+c;
    
    # calculate the first order basis functions n[i]
    
    for i=1:nplusc-1
        if (t>=x[i])&&(t<x[i+1])
            temp[i]=1;
        else
            temp[i]=0;
        end
    end
    
    if t==x[npts+1]    # handle the end specially
        temp[npts] = 1;    # resetting the first order basis functions.
        temp[npts+1]=0;
    end
    
    # calculate the higher order basis functions
    for k=2:c
        for i=1:nplusc-k
            if temp[i] != 0    # if the lower order basis function is zero skip the calculation
                b1 = ((t-x[i])*temp[i])/(x[i+k-1]-x[i]);
            else
                b1 = 0;
            end
            if temp[i+1] != 0    # if the lower order basis function is zero skip the calculation
                b2 = ((x[i+k]-t)*temp[i+1])/(x[i+k]-x[i+1]);
            else
                b2 = 0;
            end
            
            # calculate first derivative
            if temp[i] != 0    # if the lower order basis function is zero skip the calculation
                f1 = temp[i]/(x[i+k-1]-x[i]);
            else
                f1 = 0;
            end
            if temp[i+1] != 0    # if the lower order basis function is zero skip the calculation
                f2 = -temp[i+1]/(x[i+k]-x[i+1]);
            else
                f2 = 0;
            end
            if temp1[i] != 0    # if the lower order basis function is zero skip the calculation
                f3 = ((t-x[i])*temp1[i])/(x[i+k-1]-x[i]);
            else
                f3 = 0;
            end
            if temp1[i+1] != 0    # if the lower order basis function is zero skip the calculation
                f4 = ((x[i+k]-t)*temp1[i+1])/(x[i+k]-x[i+1]);
            else
                f4 = 0;
            end
            
            # calculate second derivative
            if temp1[i] != 0    # if the lower order basis function is zero skip the calculation
                s1 = (2*temp1[i])/(x[i+k-1]-x[i]);
            else
                s1 = 0;
            end
            if temp1[i+1] != 0    # if the lower order basis function is zero skip the calculation
                s2 = (-2*temp1[i+1])/(x[i+k]-x[i+1]);
            else
                s2 = 0;
            end
            if temp2[i] != 0    # if the lower order basis function is zero skip the calculation
                s3 = ((t-x[i])*temp2[i])/(x[i+k-1]-x[i]);
            else
                s3 = 0;
            end
            if temp2[i+1] != 0    # if the lower order basis function is zero skip the calculation
                s4 = ((x[i+k]-t)*temp2[i+1])/(x[i+k]-x[i+1]);
            else
                s4 = 0;
            end
            
            temp[i] = b1 + b2;
            temp1[i] = f1 + f2 + f3 + f4;
            temp2[i] = s1 + s2 + s3 + s4;
        end
    end
    
    # prepare output
    for i=1:npts
        push!(n, temp[i]);
        push!(d1, temp1[i]);
        push!(d2, temp2[i]);
    end
    return (n, d1, d2);
end