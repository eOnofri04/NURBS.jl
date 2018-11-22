export rbasis, rbspline, rbsplinu


"""
	rbasis(c, t, npts, x, h)

Generate a rational B-spline basis function using an open knot vector.
"""

function rbasis(c::Int64,t::Float64,npts::Int64,x::Array{Float64},h::Array{Float64},r::Array{Float64})
    nplusc=npts+k
    temp=zeros(nplusc)
    basis(c,t,npts,x,temp)
    sum=0
    for i=1:npts
        sum=sum+temp[i]*h[i]
    end
    for i=1:npts
        if sum!=0
            r[i]=(temp[i]*h[i])/sum
        else
            r[i]=0
        end
    end
end

#-----------------------------------------------------------------------

"""
	rbspline(npts, k, p1, b, h)

Generate a rational B-spline curve using an open uniform knot vector.
"""

function rbspline(npts::Int64,k::Int64,p1::Int64,B::Array{Float64},h::Array{Float64})
    nplusc=npts+k
    nbasis=zeros(1,npts)
    x=zeros(nplusc)
    temp=zeros(1,3)
    p=zeros(p1,3)
    
    #generate the open uniform knot vector
    knot(npts,k,x)
    icount=0
    
    #calculate the points on the rational B-spline curve
    step=x[nplusc]/(p1-1)
    for t=0:step:x[nplusc]
        icount=icount+1
        rbasis(k,t,npts,x,h,nbasis)
        temp=nbasis*B
        for j=1:3
            p[icount,j]=temp[j]
        end
    end
    return(p)
end

#-----------------------------------------------------------------------

"""
	rbsplinu(npts, k, p1, b, h)

Generate a rational B-spline curve using a periodic uniform knot vector.
"""

function rbsplineu(npts::Int64,k::Int64,p1::Int64,B::Array{Float64},h::Array{Float64})
    nplusc=npts+k
    nbasis=zeros(1,npts)
    x=zeros(nplusc)
    temp=zeros(1,3)
    p=zeros(p1,3)
    
    #generate the open uniform knot vector
    knotu(npts,k,x)
    icount=0
    
    #calculate the points on the rational B-spline curve
    step=(npts-k+1)/(p1-1)
    for t=(k-1):step:npts
        icount=icount+1
        rbasis(k,t,npts,x,h,nbasis)
        temp=nbasis*B
        for j=1:3
            p[icount,j]=temp[j]
        end
    end
    return(p)
end

