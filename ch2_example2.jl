#= This line adds functions to take
   an AR(2) model for US inflation
=#

using DataFrames
include("lag0.jl")
#load inflation datatable
df=readtable(".\\ch2_data\\dataUS.csv", header=false)

# Way of converting DataArray to Array
x1=[df[i,1] for i in [1:size(df,1)]]
x2=[df[i,2] for i in [1:size(df,1)]]
x3=[df[i,3] for i in [1:size(df,1)]]
xx=[x1 x2 x3]
T=size(xx,1)
X=[ones(T,1) lag0(xx,1) lag0(xx,2)]

Y=df
L=2 #number of lags in VAR

#Remove missing obs.
Y=Y[3:end]
X=X[3:end,:]
T=size(X,1)

# compute standard deviation of each series residual via an ols
# regression to be used in setting the prior



# #Step 1: set priors and starting values
# # priors for B
# #B0=[0;0;0]
# sigma0=eye(3)
# # priors for sigma2
# #T0=1
# #D0=0.1
#
# # starting values
# #B=B0
# #sigma2=1
#
# #reps=5000
# #burn=4000
#
# # Step 2: sample B conditional on sigma N(M*, Y*)
# isigma0=inv(sigma0)
#
# # inverse value does not exist
# temp=isigma0+(1/sigma2)*(X'*X)
# inv(temp)
