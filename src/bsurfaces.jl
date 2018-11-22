export bsplsurf, bspsurfu, dbsurf


"""
	bsplsurf()

Calculate a Cartesian product B-spline surface using open uniform knot vectors.
"""

function bsplsurf(B::Array{Float64},k::Int64,l::Int64,npts::Int64,mpts::Int64,p1::Int64,p2::Int64)::Arrayy{Float64}
    
    nplusc = npts + k
    mplusc = mpts+l
    x = zeros(nplusc)
    y = zeros(mplusc)
    nbasis = zeros(1,npts)
    mbasis = zeros(1,mpts)
    q = zeros(p1*p2,3)
    
    #generate the open uniform knot vectors
    knot(npts,k,x)
    knot(mpts,l,y)
   
    icount = 0
    #calculate the points on the B-spline surface
    stepu = x[nplusc]/(p1-1)
    stepw = y[mplusc]/(p2-1)
    for u = 0:stepu:x[nplusc]
        basis(k,u,npts,x,nbasis)
        for w = 0:stepw:y[mplusc]
            basis(l,w,mpts,y,mbasis)
            icount = icount+1
            for i = 1:npts
                for j = 1:mpts
                    j1 = mpts*(i-1)+j
                    for s = 1:3
                        q[icount,s] = q[icount,s]+B[j1,s]*nbasis[i]*mbasis[j]
                    end
                end
            end
        end
    end
    return(q)
end

#-----------------------------------------------------------------------

"""
	bspsurfu()

Calculate a Cartesian product B-spline surface using periodic uniform knot vectors.
"""

function bsplsurfu(B::Array{Float64},k::Int64,l::Int64,npts::Int64,mpts::Int64,p1::Int64,p2::Int64)::Arrayy{Float64}
    
    nplusc = npts+k
    mplusc = mpts+l
    x = zeros(nplusc)
    y = zeros(mplusc)
    nbasis = zeros(1,npts)
    mbasis = zeros(1,mpts)
    q = zeros(p1*p2,3)
    
    #generate the open uniform knot vectors
    knotu(npts,k,x)
    knotu(mpts,l,y)
   
    icount = 0
    #calculate the points on the B-spline surface
    stepu = (npts-k+1) / (p1-1)
    stepw = (mpts-k+1) / (p2-1)
    for u = k-1:stepu:npts
        basis(k,u,npts,x,nbasis)
        for w = k-1:stepw:mpts
            basis(l,w,mpts,y,mbasis)
            icount = icount+1
            for i = 1:npts
                for j = 1:mpts
                    for s = 1:3
                        j1 = mpts*(i-1)+j
                        q[icount,s] = q[icount,s]+B[j1,s]*nbasis[i]*mbasis[j]
                    end
                end
            end
        end
    end
    return(q)
end

#-----------------------------------------------------------------------

"""
	dbsurf()

Calculate a Cartesian product B-spline surface and its derivatives using open uniform knot vectors.
"""

function dbsurf(B::Array{Float64},k::Int64,l::Int64,npts::Int64,mpts::Int64,p1::Int64,p2::Int64)::Array{Float64}
    nplusc = npts + k
    mplusc = mpts+l
    x = zeros(nplusc)
    y = zeros(mplusc)
    nbasis = zeros(1,npts)
    mbasis = zeros(1,mpts)
    d1nbasis = zeros(1,npts)
    d1mbasis = zeros(1,mpts)
    d2nbasis = zeros(1,npts)
    d2mbasis = zeros(1,mpts)
    q = zeros(p1*p2,3)
    qu = zeros(p1*p2,3)
    qw = zeros(p1*p2,3)
    quu = zeros(p1*p2,3)
    quw = zeros(p1*p2,3)
    qww = zeros(p1*p2,3)
    
    knot(npts,k,x)
    knot(mpts,l,y)
    
    icount = 0
    #calculate the points on the B-spline surface
    stepu = x[nplusc]/(p1-1)
    stepw = y[mplusc]/(p2-1)
    for u = 0:stepu:x[nplusc]
        dbasis(k,u,npts,x,nbasis,d1nbasis,d2nbasis)
        for w = 0:stepw:y[mplusc]
            dbasis(l,w,mpts,y,mbasis,d1mbasis,d2mbasis)
            icount = icount+1
            for i = 1:npts
                for j = 1:mpts                    
                    for s = 1:3
                        j1 = mpts*(i-1)+j
                        q[icount,s] = q[icount,s] + B[j1,s] * nbasis[i] * mbasis[j]
                        qu[icount,s] = qu[icount,s] + B[j1,s] * d1nbasis[i] * mbasis[j]
                        qw[icount,s] = qw[icount,s] + B[j1,s] * nbasis[i] * d1mbasis[j]
                        quu[icount,s] = quu[icount,s] + B[j1,s] * d1nbasis[i] * d1mbasis[j]
                        quw[icount,s] = quw[icount,s] + B[j1,s] * d2nbasis[i] * mbasis[j]
                        qww[icount,s] = qww[icount,s] + B[j1,s] * nbasis[i] * d2mbasis[j]
                    end
                end
            end
        end
    end
    return(q)

end
