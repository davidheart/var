
function lead0(x,p)

#Compute the number of rows and columns of input x and save in
#variables R and C respectively. 

     R::Int32=size(x,1)
     C::Int32=size(x,2)

# Take the first R-p rows of matrix x
     x1=x[p+1:end, :]
     return out=[x1;zeros(p,C)]

end


## Testing 
#x=[1 2 3; 3 4 5; 6 7 8]
#print(lead0(x,2))



