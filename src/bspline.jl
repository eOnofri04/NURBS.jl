export bspline, bsplineu, dbspline, dbsplineu


"""
	bsplfit(dpts, d, npts, k)

Fit a B-spline curve using an open uniform knot vector.
"""

function bsplfit(dpts::Int64, d::Array{Float64}, npts::Int64, k::Int64)::Array{Float64}
	[...]
end

#-----------------------------------------------------------------------

"""
	bspline(npts, k, p1, b)

Generate a B-spline curve using an open uniform knot vector.
"""

function bspline(npts::Int64,k::Int64,p1::Int64,B::Array{Float64})
         
         nplusc=npts+k
         
         nbasis=zeros(1,npts)
         x=zeros(nplusc)
         
         knot(npts,k,x)
         icount=0
         temp=zeros(1,3)
         p=zeros(p1,3)
         step=x[nplusc]/(p1-1)
         for t=0:step:x[nplusc]
             basis(k,t,npts,x,nbasis)
             temp=nbasis*B
             icount=icount+1
             for j=1:3
                 p[icount,j]=temp[j]
             end
         end
         return(p)
end

#-----------------------------------------------------------------------

"""
	bsplineu(npts k, p1, b)

Generate a B-spline curve using a periodic uniform knot vector.
"""

function bsplineu(npts::Int64,k::Int64,p1::Int64,B::Array{Float64})
         
         nplusc=npts+k
         
         nbasis=zeros(1,npts)
         x=zeros(nplusc)
         
         knotu(npts,k,x)
         
         icount=0
         temp=zeros(1,3)
         p=zeros(p1,3)
         step=(npts-k+1)/(p1-1)
         
         for t=k-1:step:npts
             basis(k,t,npts,x,nbasis)
             temp=nbasis*B
             icount=icount+1
             for j=1:3
                 p[icount,j]=temp[j]
             end
         end
         return(p)
end

#-----------------------------------------------------------------------

"""
	dbspline(npts k, p1, b)

Generate a B-spline curve and its derivatives using an open uniform knot vector.
"""

function dbspline(npts::Int64,k::Int64,p1::Int64,B::Array{Float64})
    #zero and redimension the knot vector, basis, curve and derivative arrays
    p=zeros(p1,3)
    d1=zeros(p1,3)
    d2=zeros(p1,3)
    nbasis=zeros(1,npts)
    d1nbasis=zeros(1,npts)
    d2nbasis=zeros(1,npts)
    
    nplusc=npts+k
    x=zeros(nplusc)
    
    #generate the uniform open knot vector
    knot(npts,k,x)
    icount=0
    
    #calculate the points on the B-spline curve and their first and second derivatives
    step=x[nplusc]/(p1-1)
    for t=0:step:x[nplusc]
        icount=icount+1
        if icount==p1 t=x[nplusc] end
        dbasis(k,t,npts,x,nbasis,d1nbasis,d2nbasis)
        temp=nbasis*B
        temp1=d1nbasis*B
        temp2=d2nbasis*B
        for j=1:3
            p[icount,j]=temp[j]
            d1[icount,j]=temp[j]
            d2[icount,j]=temp[j]
        end
    end
    return(p)
end

#-----------------------------------------------------------------------

"""
	dbsplineu(npts k, p1, b)

Generate a B-spline curve and its derivatives using an open uniform knot vector .
"""

function dbsplineu(npts::Int64,k::Int64,p1::Int64,B::Array{Float64})
    #zero and redimension the knot vector, basis, curve and derivative arrays
    p=zeros(p1,3)
    d1=zeros(p1,3)
    d2=zeros(p1,3)
    nbasis=zeros(1,npts)
    d1nbasis=zeros(1,npts)
    d2nbasis=zeros(1,npts)
    
    nplusc=npts+k
    x=zeros(nplusc)
    
    #generate the uniform open knot vector
    knotu(npts,k,x)
    icount=0
    
    #calculate the points on the B-spline curve and their first and second derivatives
    step=(npts-k+1)/(p1-1)
    for t=(k-1):step:npts
        icount=icount+1
        if icount==p1 t=Float64(npts) end
        dbasisu(k,t,npts,x,nbasis,d1nbasis,d2nbasis)
        temp=nbasis*B
        temp1=d1nbasis*B
        temp2=d2nbasis*B
        for j=1:3
            p[icount,j]=temp[j]
            d1[icount,j]=temp[j]
            d2[icount,j]=temp[j]
        end
    end
    return(p)
end

#-----------------------------------------------------------------------

"""
	matpbspl(npts k, p1, b)

Generate a B-spline curve using matrix methods and a periodic uniform knot vector.
"""

function matpbspl(npts::Int64, k::Int64, p1::Int64, b::Array{Float64})::Array{Float64}
	[...]
end

#-----------------------------------------------------------------------

"""
	nmatrix(k)

Calculate the general B-spline periodic basis matrix.

This function is used in matpbspl.
"""

function namtrix(k::Int64)::tuple{Float64, Array{Int64}}
	[...]
end

#-----------------------------------------------------------------------

"""
	param(dpts, d)

Calculate parameter values based on chord distances.
"""

function param(dpts::Int64, d::Array{Float64})::Array{Float64}
	[...]
end


