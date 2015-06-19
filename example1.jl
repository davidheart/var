#= This line adds functions to take 
   an AR(2) model for US inflation 
=#
#load inflation datatable
using DataFrames
include("lag0.jl")
# "LOAD_PATH" is a system name. 
#@everywhere push!(LOAD_PATH, "/home/sungguanyun/var/")
df=readtable("inflation.csv")
# way to DataArray to Array
Y=eltype(df[1,2])[df[i,2] for i in [1:size(df,1)]]
# or use Y=df[2] instead
T=size(Y,1)
X=[ones(T,1) lag0(Y,1) lag0(Y,2)]
#Remove missing obs.
Y=Y[3:end]
X=X[3:end,:]
T=size(X,1)
#Step 1: set priors and starting values 
# priors for B
B0=[0;0;0]
sigma0=eye(3)
# priors for sigma2
T0=1
D0=0.1
# starting values
B=B0
sigma2=1.0

reps=5000
burn=4000

# Still wondering how to express in julia equivalently to x=[ ] in Matlab
out1=zeros(3,3)
out2=zeros(1,1)




for i in [1:reps]

  # Step 2: sample B conditional on sigma N(M*, Y*)
  M=inv(inv(sigma0)+(1/sigma2)*(X'*X))*(inv(sigma0)*B0+(1/sigma2)*(X'*Y))
  V=inv(inv(sigma0)+(1/sigma2)*(X'*X))
  
  chck=-1   
     while chck<0
     	B=M+(randn(1,3)*chol(V))'
     	b=[B[2] B[3];1 0]
     	pnum=abs(eigvals(b))
     	ee=max(pnum[1], pnum[2])
     	if ee<=1
     		chck=1
     	end
     end

  # Step 3 sample sigma2 conditional on B from IG(T1,D1)
  # compute residuals
  resids=Y-X*B
  # compute posterior df and scale matrix
  T1=T0+T
  D1=D0+resids'*resids
   # Draw from IG
   z0=randn(T1,1)
   z0z0=z0'*z0
   # Julia issue: typeof(sigma2) now turns into Array
   # which is not acceptable for the operator, "/"
   # sound though quriky, use maximum to convert one-by-one
   # arrary to scalar
   sigma2=maximum(D1/z0z0)

   if i>burn
   	out1=vcat(out1, B')
	out2=vcat(out2,sigma2)
   end

end
