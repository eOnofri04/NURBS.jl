export bezsurfj, bezsurf1, bezsurf2, bezsurf3, runn, runner


#=== RUN ===#

function runner()
	runn()
	addprocs(9)
	runn()
end

function runn()
	npts = 50
    mpts = 50
    udpts = 2500
    wdpts = 2500
    bplus = zeros(3, npts * mpts)

    for i = 1 : npts
        ki = i/2
        for j = 1 : mpts
            kj = j /2
            r = sqrt((ki)^2+(kj)^2)
            if (r == 0)
                r = 0.000001
            end
            ij = (i-1)*npts + j
            bplus[1,ij] = i
            bplus[2,ij] = j
            bplus[3,ij] = i+j
        end
    end

	avgj = 0.0
	avg1 = 0.0
	avg2 = 0.0
	avg3 = 0.0

	for i = 1 : 50
	    xj = @elapsed bezsurfj(npts, mpts, bplus, udpts, wdpts)
	    x1 = @elapsed bezsurf1(npts, mpts, bplus, udpts, wdpts)
	    x2 = @elapsed bezsurf2(npts, mpts, bplus, udpts, wdpts)
	    x3 = @elapsed bezsurf3(npts, mpts, bplus, udpts, wdpts)
	    avgj = avgj + xj
	    avg1 = avg1 + x1
	    avg2 = avg2 + x2
	    avg3 = avg3 + x3
	end
	avg1 /= 50
	avg2 /= 50
	avg3 /= 50
	avgj /= 50

	println("average julia = $(avgj)")
	println("average par 1 = $(avg1)")
	println("average par 2 = $(avg2)")
	println("average par 3 = $(avg3)")

end

#=== DEFAULT ===#
function bezsurfj(npts::Int64, mpts::Int64, b::Array{Float64,2}, udpts::Int64, wdpts::Int64)::Tuple{Array{Float64,2}, Array{Array{Int64,1},1}, Array{Array{Int64,1},1}}

	@assert npts*mpts == length(b[1,:]) ("ERROR: there are not npts * mpts = $(npts) * $(mpts) = $(npts * mpts) points of the net in b[]")

	ustep::Float64 = 1 / (udpts-1)
	wstep::Float64 = 1 / (wdpts - 1)
	n::Int64 = npts - 1
	m::Int64 = mpts - 1

	#=== x, y, z components of B ===#
	x::Array{Float64,2} = zeros(npts, mpts)
	y::Array{Float64,2} = zeros(npts, mpts)
	z::Array{Float64,2} = zeros(npts, mpts)

	for i = 1 : npts
		for j = 1 : mpts
			ij = (i - 1) * mpts + j
			x[i, j] = b[1, ij]
			y[i, j] = b[2, ij]
			z[i, j] = b[3, ij]
		end
	end

	#=== N = C * D ===#
	C::Array{Float64,2} = zeros(npts, npts)
	for i = 0:n
		for j = 0 : (n - i)
			if (n-j-i)%2==0
				C[i+1, j+1] = binomial(n-j, n-j-i)
			else
				C[i+1, j+1] = 0.0-binomial(n-j, n-j-i)
			end
		end
	end

	D::Array{Float64,2} = zeros(npts, npts)
	for i = 1 : npts
		D[i, i] = abs(C[i, 1])
	end

	#=== M = E * F ===#
	E::Array{Float64,2} = zeros(mpts, mpts)
	for i = 0 : m
		for j = 0 : (m - i)
			if (m - j - i) % 2 == 0
				E[i+1, j+1] = binomial(m-j, m-j-i)
			else
				E[i+1, j+1] = 0.0 - binomial(m-j, m-j-i)
			end
		end
	end

	F::Array{Float64,2} = zeros(mpts, mpts)
	for i = 1 : mpts
		F[i, i] = abs(E[i, 1])
	end

	#=== U = [u^n ... u 1] ===#
	U::Array{Float64,2} = ones(udpts, npts)

	u::Float64 = 0.0
	for i = 1 : udpts
		for j = 0 : n-1  #Per j=0 il valore è 1 ed è già impostato
			U[i, n-j] = U[i, npts-j] * u
		end
		u = u + ustep
	end

	#=== W = [w^m ... w 1] ===#
	W::Array{Float64,2} = ones(wdpts, mpts)

	w::Float64 = 0.0
	for i = 1 : wdpts
		for j = 0 : m-1  #Per j=0 il valore è 1 ed è già impostato
			W[i, m-j] = W[i, mpts-j] * w
		end
		w = w + wstep
	end

	#=== Edge Vector ===#
	EV::Array{Array{Int64,1},1} = Array{Int64,1}[] #1-dimensional cellular complex used for plotting the curve with Plasm package

	for i = 1 : ((udpts * wdpts) - 1)
		app = Array{Int64}(2)
		app[1] = i
		app[2] = i+1
		push!(EV, app) 
	end

	for i = 1 : (wdpts * (udpts - 1))
		app = Array{Int64}(2)
		app[1] = i
		app[2] = i + wdpts
		push!(EV, app) 
	end

	#=== Faces Vector ===#
	FV::Array{Array{Int64,1},1} = Array{Int64,1}[] #1-dimensional cellular complex used for plotting the curve with Plasm package

    for j = 0 : (udpts - 2)
        for i = 1 : wdpts - 1
            app = Array{Int64}(4)
            app[1] = j * wdpts + i
            app[2] = j * wdpts + i + 1
            app[3] = (j + 1) * wdpts + i
            app[4] = (j + 1) * wdpts + i + 1
            push!(FV, app) 
        end
    end

    #=== Points Vector ===#
    UCD::Array{Float64,2} = U * C * D
    FEW::Array{Float64,2} = transpose(W * E * F)
	xP::Array{Float64,2} = transpose(UCD * x * FEW)
	yP::Array{Float64,2} = transpose(UCD * y * FEW)
	zP::Array{Float64,2} = transpose(UCD * z * FEW)

	P::Array{Float64,2} = zeros(3, udpts * wdpts)

	for i = 1 : udpts
		for j = 1 : wdpts
			ij::Int64 = (i-1) * udpts + j
			P[1, ij] = xP[i, j]
			P[2, ij] = yP[i, j]
			P[3, ij] = zP[i, j]
		end
	end

	return (P, EV, FV)
end













#=== PARALLEL 1 ===#
function bezsurf1(npts::Int64, mpts::Int64, b::Array{Float64,2}, udpts::Int64, wdpts::Int64)::Tuple{Array{Float64,2}, Array{Array{Int64,1},1}, Array{Array{Int64,1},1}}

	@assert npts*mpts == length(b[1,:]) ("ERROR: there are not npts * mpts = $(npts) * $(mpts) = $(npts * mpts) points of the net in b[]")

	ustep::Float64 = 1 / (udpts-1)
	wstep::Float64 = 1 / (wdpts - 1)
	n::Int64 = npts - 1
	m::Int64 = mpts - 1

	#=== x, y, z components of B ===#
	x::Array{Float64,2} = zeros(npts, mpts)
	y::Array{Float64,2} = zeros(npts, mpts)
	z::Array{Float64,2} = zeros(npts, mpts)

	for i = 1 : npts
		for j = 1 : mpts
			ij = (i - 1) * mpts + j
			x[i, j] = b[1, ij]
			y[i, j] = b[2, ij]
			z[i, j] = b[3, ij]
		end
	end

	#=== SYNC BLOCK BEGIN ===#
	#=== all the following @async are totally indipendent blocks ===#
	@sync begin
	
		#=== N = C * D ===#
		C::Array{Float64,2} = zeros(npts, npts)
		D::Array{Float64,2} = zeros(npts, npts)

		@async begin
			for i = 0 : n
				for j = 0 : (n - i)
					if (n-j-i)%2==0
						C[i+1, j+1] = binomial(n-j, n-j-i)
					else
						C[i+1, j+1] = 0.0-binomial(n-j, n-j-i)
					end
				end
			end

			for i = 1 : npts
				D[i, i] = abs(C[i, 1])
			end
		end

		#=== M = E * F ===#
		E::Array{Float64,2} = zeros(mpts, mpts)
		F::Array{Float64,2} = zeros(mpts, mpts)

		@async begin
			for i = 0 : m
				for j = 0 : (m - i)
					if (m - j - i) % 2 == 0
						E[i+1, j+1] = binomial(m-j, m-j-i)
					else
						E[i+1, j+1] = 0.0 - binomial(m-j, m-j-i)
					end
				end
			end

			for i = 1 : mpts
				F[i, i] = abs(E[i, 1])
			end
		end

		#=== U = [u^n ... u 1] ===#
		U::Array{Float64,2} = ones(udpts, npts)

		@async begin
			u::Float64 = 0.0
			for i = 1 : udpts
				for j = 0 : n-1  #Per j=0 il valore è 1 ed è già impostato
					U[i, n-j] = U[i, npts-j] * u
				end
				u = u + ustep
			end
		end

		#=== W = [w^m ... w 1] ===#
		W::Array{Float64,2} = ones(wdpts, mpts)

		@async begin
			w::Float64 = 0.0
			for i = 1 : wdpts
				for j = 0 : m-1  #Per j=0 il valore è 1 ed è già impostato
					W[i, m-j] = W[i, mpts-j] * w
				end
				w = w + wstep
			end
		end

		#=== Edge Vector ===#
		EV::Array{Array{Int64,1},1} = Array{Int64,1}[] #1-dimensional cellular complex used for plotting the curve with Plasm package

		@async begin
			for i = 1 : ((udpts * wdpts) - 1)
				app = Array{Int64}(2)
				app[1] = i
				app[2] = i+1
				push!(EV, app) 
			end

			for i = 1 : (wdpts * (udpts - 1))
				app = Array{Int64}(2)
				app[1] = i
				app[2] = i + wdpts
				push!(EV, app) 
			end
		end

		#=== Faces Vector ===#
		FV::Array{Array{Int64,1},1} = Array{Int64,1}[] #1-dimensional cellular complex used for plotting the curve with Plasm package

		@async begin
			for j = 0 : (udpts - 2)
				for i = 1 : wdpts - 1
					app = Array{Int64}(4)
					app[1] = j * wdpts + i
					app[2] = j * wdpts + i + 1
					app[3] = (j + 1) * wdpts + i
					app[4] = (j + 1) * wdpts + i + 1
					push!(FV, app) 
				end
			end
		end
	
	end
	#=== SYNC BLOCK END ===#

	#=== Points Vector ===#
	UCD::Array{Float64,2} = U * C * D
	FEW::Array{Float64,2} = transpose(W * E * F)
	xP::Array{Float64,2} = transpose(UCD * x * FEW)
	yP::Array{Float64,2} = transpose(UCD * y * FEW)
	zP::Array{Float64,2} = transpose(UCD * z * FEW)

	P::Array{Float64,2} = zeros(3, udpts * wdpts)

	for i = 1 : udpts
		for j = 1 : wdpts
			ij::Int64 = (i-1) * udpts + j
			P[1, ij] = xP[i, j]
			P[2, ij] = yP[i, j]
			P[3, ij] = zP[i, j]
		end
	end

	return (P, EV, FV)
end











#=== PARALLEL 2 ===#
function bezsurf2(npts::Int64, mpts::Int64, b::Array{Float64,2}, udpts::Int64, wdpts::Int64)::Tuple{Array{Float64,2}, Array{Array{Int64,1},1}, Array{Array{Int64,1},1}}

	@assert npts*mpts == length(b[1,:]) ("ERROR: there are not npts * mpts = $(npts) * $(mpts) = $(npts * mpts) points of the net in b[]")

	ustep::Float64 = 1 / (udpts-1)
	wstep::Float64 = 1 / (wdpts - 1)
	n::Int64 = npts - 1
	m::Int64 = mpts - 1

	#=== x, y, z components of B ===#
	x::Array{Float64,2} = zeros(npts, mpts)
	y::Array{Float64,2} = zeros(npts, mpts)
	z::Array{Float64,2} = zeros(npts, mpts)

	for i = 1 : npts
		for j = 1 : mpts
			ij = (i - 1) * mpts + j
			x[i, j] = b[1, ij]
			y[i, j] = b[2, ij]
			z[i, j] = b[3, ij]
		end
	end
    
    
    C::Array{Float64,2} = zeros(npts, npts)
    D::Array{Float64,2} = zeros(npts, npts)
    E::Array{Float64,2} = zeros(mpts, mpts)
    F::Array{Float64,2} = zeros(mpts, mpts)
    U::Array{Float64,2} = ones(udpts, npts)
    W::Array{Float64,2} = ones(wdpts, mpts)
    EV::Array{Array{Int64,1},1} = Array{Int64,1}[] #1-dimensional cellular complex used for plotting the curve with Plasm package
    FV::Array{Array{Int64,1},1} = Array{Int64,1}[] #1-dimensional cellular complex used for plotting the curve with Plasm package

	#=== SYNC BLOCK BEGIN ===#
	#=== all the following @async are totally indipendent blocks ===#
	@sync begin
	
		#=== N = C * D ===#
		@async begin
			for i = 0 : n
				for j = 0 : (n - i)
					if (n-j-i)%2==0
						C[i+1, j+1] = binomial(n-j, n-j-i)
					else
						C[i+1, j+1] = 0.0-binomial(n-j, n-j-i)
					end
				end
			end

			for i = 1 : npts
				D[i, i] = abs(C[i, 1])
			end
		end

		#=== M = E * F ===#
        @async begin
			for i = 0 : m
				for j = 0 : (m - i)
					if (m - j - i) % 2 == 0
						E[i+1, j+1] = binomial(m-j, m-j-i)
					else
						E[i+1, j+1] = 0.0 - binomial(m-j, m-j-i)
					end
				end
			end

			for i = 1 : mpts
				F[i, i] = abs(E[i, 1])
			end
		end

		#=== U = [u^n ... u 1] ===#
		@async begin
			u::Float64 = 0.0
			for i = 1 : udpts
				for j = 0 : n-1  #Per j=0 il valore è 1 ed è già impostato
					U[i, n-j] = U[i, npts-j] * u
				end
				u = u + ustep
			end
		end

		#=== W = [w^m ... w 1] ===#
		@async begin
			w::Float64 = 0.0
			for i = 1 : wdpts
				for j = 0 : m-1  #Per j=0 il valore è 1 ed è già impostato
					W[i, m-j] = W[i, mpts-j] * w
				end
				w = w + wstep
			end
		end

		#=== Edge Vector ===#
		@async begin
			for i = 1 : ((udpts * wdpts) - 1)
				app = Array{Int64}(2)
				app[1] = i
				app[2] = i+1
				push!(EV, app) 
			end

			for i = 1 : (wdpts * (udpts - 1))
				app = Array{Int64}(2)
				app[1] = i
				app[2] = i + wdpts
				push!(EV, app) 
			end
		end

		#=== Faces Vector ===#
		
		@async begin
			for j = 0 : (udpts - 2)
				for i = 1 : wdpts - 1
					app = Array{Int64}(4)
					app[1] = j * wdpts + i
					app[2] = j * wdpts + i + 1
					app[3] = (j + 1) * wdpts + i
					app[4] = (j + 1) * wdpts + i + 1
					push!(FV, app) 
				end
			end
		end
	
	end
	#=== SYNC BLOCK END ===#

	#=== Points Vector ===#
	UCD::Array{Float64,2} = U * C * D
	FEW::Array{Float64,2} = transpose(W * E * F)
	xP::Array{Float64,2} = transpose(UCD * x * FEW)
	yP::Array{Float64,2} = transpose(UCD * y * FEW)
	zP::Array{Float64,2} = transpose(UCD * z * FEW)

	P::Array{Float64,2} = zeros(3, udpts * wdpts)

	for i = 1 : udpts
		for j = 1 : wdpts
			ij::Int64 = (i-1) * udpts + j
			P[1, ij] = xP[i, j]
			P[2, ij] = yP[i, j]
			P[3, ij] = zP[i, j]
		end
	end

	return (P, EV, FV)
end














#=== PARALLEL 3 ===#
function bezsurf3(npts::Int64, mpts::Int64, b::Array{Float64,2}, udpts::Int64, wdpts::Int64)::Tuple{Array{Float64,2}, Array{Array{Int64,1},1}, Array{Array{Int64,1},1}}

	@assert npts*mpts == length(b[1,:]) ("ERROR: there are not npts * mpts = $(npts) * $(mpts) = $(npts * mpts) points of the net in b[]")

	ustep::Float64 = 1 / (udpts-1)
	wstep::Float64 = 1 / (wdpts - 1)
	n::Int64 = npts - 1
	m::Int64 = mpts - 1

	#=== x, y, z components of B ===#
	x::Array{Float64,2} = zeros(npts, mpts)
	y::Array{Float64,2} = zeros(npts, mpts)
	z::Array{Float64,2} = zeros(npts, mpts)

	for i = 1 : npts
		for j = 1 : mpts
			ij = (i - 1) * mpts + j
			x[i, j] = b[1, ij]
			y[i, j] = b[2, ij]
			z[i, j] = b[3, ij]
		end
	end

	#=== SYNC BLOCK BEGIN ===#
	#=== all the following @async are totally indipendent blocks ===#
	@sync begin
	
		#=== N = C * D ===#
		C::Array{Float64,2} = zeros(2, 2)
		D::Array{Float64,2} = zeros(2, 2)

		@async begin
			C = zeros(npts, npts)
			D = zeros(npts, npts)
			for i = 0 : n
				for j = 0 : (n - i)
					if (n-j-i)%2==0
						C[i+1, j+1] = binomial(n-j, n-j-i)
					else
						C[i+1, j+1] = 0.0-binomial(n-j, n-j-i)
					end
				end
			end

			for i = 1 : npts
				D[i, i] = abs(C[i, 1])
			end
		end

		#=== M = E * F ===#
		E::Array{Float64,2} = zeros(2, 2)
		F::Array{Float64,2} = zeros(2, 2)

		@async begin
			E = zeros(mpts, mpts)
			F = zeros(mpts, mpts)
			for i = 0 : m
				for j = 0 : (m - i)
					if (m - j - i) % 2 == 0
						E[i+1, j+1] = binomial(m-j, m-j-i)
					else
						E[i+1, j+1] = 0.0 - binomial(m-j, m-j-i)
					end
				end
			end

			for i = 1 : mpts
				F[i, i] = abs(E[i, 1])
			end
		end

		#=== U = [u^n ... u 1] ===#
		U::Array{Float64,2} = zeros(2, 2)

		@async begin
			U = ones(udpts, npts)
			u::Float64 = 0.0
			for i = 1 : udpts
				for j = 0 : n-1  #Per j=0 il valore è 1 ed è già impostato
					U[i, n-j] = U[i, npts-j] * u
				end
				u = u + ustep
			end
		end

		#=== W = [w^m ... w 1] ===#
		W::Array{Float64,2} = zeros(2, 2)

		@async begin
			W = ones(wdpts, mpts)
			w::Float64 = 0.0
			for i = 1 : wdpts
				for j = 0 : m-1  #Per j=0 il valore è 1 ed è già impostato
					W[i, m-j] = W[i, mpts-j] * w
				end
				w = w + wstep
			end
		end

		#=== Edge Vector ===#
		EV::Array{Array{Int64,1},1} = Array{Int64,1}[] #1-dimensional cellular complex used for plotting the curve with Plasm package

		@async begin
			for i = 1 : ((udpts * wdpts) - 1)
				app = Array{Int64}(2)
				app[1] = i
				app[2] = i+1
				push!(EV, app) 
			end

			for i = 1 : (wdpts * (udpts - 1))
				app = Array{Int64}(2)
				app[1] = i
				app[2] = i + wdpts
				push!(EV, app) 
			end
		end

		#=== Faces Vector ===#
		FV::Array{Array{Int64,1},1} = Array{Int64,1}[] #1-dimensional cellular complex used for plotting the curve with Plasm package

		@async begin
			for j = 0 : (udpts - 2)
				for i = 1 : wdpts - 1
					app = Array{Int64}(4)
					app[1] = j * wdpts + i
					app[2] = j * wdpts + i + 1
					app[3] = (j + 1) * wdpts + i
					app[4] = (j + 1) * wdpts + i + 1
					push!(FV, app) 
				end
			end
		end
	
	end
	#=== SYNC BLOCK END ===#

	#=== Points Vector ===#
	UCD::Array{Float64,2} = U * C * D
	FEW::Array{Float64,2} = transpose(W * E * F)
	xP::Array{Float64,2} = transpose(UCD * x * FEW)
	yP::Array{Float64,2} = transpose(UCD * y * FEW)
	zP::Array{Float64,2} = transpose(UCD * z * FEW)

	P::Array{Float64,2} = zeros(3, udpts * wdpts)

	for i = 1 : udpts
		for j = 1 : wdpts
			ij::Int64 = (i-1) * udpts + j
			P[1, ij] = xP[i, j]
			P[2, ij] = yP[i, j]
			P[3, ij] = zP[i, j]
		end
	end

	return (P, EV, FV)
end