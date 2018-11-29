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
	rbspsurf(b::Array{Float64}, k::Int64, l::Int64, npts::Int64, mpts::Int64, p1::Int64, p2::Int64)

Calculate a Cartesian product rational B-spline surface (NURBS) using an open uniform knot vector.

Call: basis, knot, sumrbas
"""

function rbspsurf(B::Array{Float64},ordx::Int64,ordy::Int64,npts::Int64,mpts::Int64,p1::Int64,p2::Int64)::Array{Float64}
    
    nplusc = npts + ordx
    mplusc = mpts + ordy
    x = zeros(nplusc)
    y = zeros(mplusc)
    nbasis = zeros(1,npts)
    mbasis = zeros(1,mpts)
    q = zeros(p1*p2,3)
    
    #generate the open uniform knot vectors
    x = knot(npts,ordx)
    y = knot(mpts,ordy)
   
    icount = 0
    #calculate the points on the B-spline surface
    stepu = x[nplusc]/(p1-1)
    stepw = y[mplusc]/(p2-1)
    for u = 0:stepu:x[nplusc]
        nbasis = basis(ordx,u,npts,x)
        for w = 0:stepw:y[mplusc]
            mbasis = basis(ordy,w,mpts,y)
            sum = sumrbas(B,nbasis,mbasis,npts,mpts)
            icount = icount+1
            for i = 1:npts
                for j = 1:mpts
                    for s = 1:3
                        j1 = mpts * (i-1) + j
                        qtemp = (B[j1,4] * B[j1,s] * nbasis[i] * mbasis[j]) / sum
                        q[icount,s] = q[icount,s] + qtemp
                    end
                end
            end
        end
    end
    return(q)
end

#-----------------------------------------------------------------------

"""
	sumrbas(b, nbasis, mbasis, npts, mpts)

Calculate the sum of the nonrational basis functions.
"""

function sumrbas(B::Array{Float64},nbasis::Array{Float64},mbasis::Array{Float64},npts::Int64,mpts::Int64)::Float64
    sum = 0
    for i = 1:npts
        for j = 1:mpts
            j1 = mpts * (i-1) + j
            sum = sum + B[j1,n] * nbasis[i] * mbasis[j]
        end
    end
    return(sum)
end
