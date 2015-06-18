


function hpfilter(y::Vector{Float64}, lambda::Float64)

#=
Source : Sebastien Villemot
=#

         n=length(y)
	 @assert n>=4

	 diag2=lambda*ones(n-2)
	 diag1=[-2lambda; -4lambda*ones(n-3); -2lambda]
	 diag0=[1+lambda; 1+5lambda; (1+6lambda)*ones(n-4);
	        1+5lambda;1+lambda]
	 D=spdiagm((diag2, diag1, diag0, diag1, diag2), (-2,-1,0,1,2))
	 D\y
end


# Test 
#y=cumsum(randn(40))
#tau=hpfilter(y, 1600.)









