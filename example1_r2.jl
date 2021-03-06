
#= This line adds functions to take 
   an AR(2) model for US inflation 
=#

using DataFrames


function lag0(x,p)
     R::Int32=size(x,1)
     C::Int32=size(x,2)

# Take the first R-p rows of matrix x
     x1=x[1:(R-p),:]
     return out=[zeros(p,C); x1]

end




#load inflation datatable
df=readtable("inflation.csv")

# Way of converting DataArray to Array
Y=[df[i,2] for i in [1:size(df,1)]]
T=size(Y,1)
X=[ones(T,1) lag0(Y,1) lag0(Y,2)]

#Remove missing obs.
Y=Y[3:end]
X=X[3:end,:]
T=size(X,1)

#Step 1: set priors and starting values 
# priors for B
#B0=[0;0;0]
sigma0=eye(3)

# priors for sigma2
#T0=1
#D0=0.1

# starting values
#B=B0
#sigma2=1

#reps=5000
#burn=4000

# Step 2: sample B conditional on sigma N(M*, Y*)
isigma0=inv(sigma0)

# inverse value does not exist 
temp=isigma0+(1/sigma2)*(X'*X)
inv(temp)













