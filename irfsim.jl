

function irfsim(b, n, l, v, s, t)
	#=
          n=number of variables
	  l=lag length
	  v=A0 matrix
	  s=shock vector
	  t=horizon
	=#

  e=zeros[t+l,n]
  e[l+1,:]=s

  y=zeros[t+1,n]

  for k in [l+1:t]
	  x=zeros( # revisit later to this point

	  for i in [1:l]
		  for j in [1:n]
			  x=[x y[k-i,j]]
		  end
	   end 
	   y[k,:]=([x 1]*b)+(e[k,:]*v)
   end
   y=y[l+1:rows(y)-l,:]

end



