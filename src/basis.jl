export basis, dbasis, dbasisu


"""
	basis(c, t, npts, x)

Generate B-spline basis functions for open uniform knot vectors.
"""

function basis(c::Int64,t::Float64,npts::Int64,x::Array{Float64},n::Array{Float64})
    nplusc=npts+c    
    temp=zeros(nplusc)
    
    
#calculate the first-order basis functions Ni,1
for i=1:nplusc-1
    if t>=x[i] && t<x[i+1]
        temp[i]=1
    else
        temp[i]=0
    end
end

    
##calculate the higher-order basis functions
for k=2:c
    for i=1:nplusc-k
        if temp[i]!=0   #if basis function is zero skip the calculation
            d=((t-x[i])*temp[i])/(x[i+k-1]-x[i])
        else
            d=0
        end
        if temp[i+1]!=0   #if basis function is zero skip the calculation
            e=((x[i+k]-t)*temp[i+1])/(x[i+k]-x[i+1])
        else
            e=0
        end
        temp[i]=d+e
    end
end

if t==x[nplusc]  #pick up last point
    temp[npts]=1
end

    
#put in n array
for i=1:npts
    n[i]=temp[i]
end


if t==x[nplusc]    #pick up last point
    n[npts]=1
end

end

#-----------------------------------------------------------------------

"""
	dbasis(c, t, npts, x)

Generate B-spline basis functions and their derivatives for uniform open knot vectors.
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
#-----------------------------------------------------------------------

"""
	dbasisu(c, t, npts, x)

Generate B-spline basis functions and their derivatives for uniform periodic knot vectors.
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
            
    
end
