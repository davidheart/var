
function lag0(x,p)

#=
For an input vector or matrix x=[a_{1}; a_{2},...,a_{n}) where 
a_{i} are row vectors, it turns the vector or matrix:

log0(x,10)=[0; ... ; 0; a_{1}; a_{2},...,a_{n}] 

In other words, it lags the variable p periods and places zeros in the 
rows corresponding to the first p periods.

Verion 0.1 in June 2015.
=#

#Compute the number of rows and columns of input x and save in
#variables R and C respectively. 

     R::Int32=size(x,1)
     C::Int32=size(x,2)

# Take the first R-p rows of matrix x
     x1=x[1:(R-p),:]
     return out=[zeros(p,C); x1]

end


## Testing 
#x=[1 2 3; 3 4 5; 6 7 8]
#print(lag0(x,2))



