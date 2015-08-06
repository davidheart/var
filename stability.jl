

#=
beta=[0.1;0.2;0.4;0.5;0.6;0.7]
n=2
l=1
=#

function stability(beta,n,l)
  # coef   (n*l+1)xn matrix with the coef from the VAR
  # l      number of lags
  # n      number of endog variables
  # FF     matrix with all coef
  # S      dummy var: if equal one->stability
  coef=reshape(beta,n*l+1,n)
  #coef
  FF=zeros(n*l, n*l)
  FF[n+1:n*l,1:n*(l-1)]=eye(n*(l-1),n*(l-1))
  temp=reshape(beta,n*l+1,n)
  temp=temp[2:n*l+1,1:n]'
  FF[1:n,1:n*l]=temp
  ee=abs(eigvals(FF))
  ee=max(ee[1], ee[2])
  S=ee>1
end
