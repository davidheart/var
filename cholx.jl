

#= In Matlab, [R,p]=chol(A) with two output arguments, never
   produces error message. If A is positive definite, then p is 0 
   and R is the same as above. But if A is not positive definite,  
   then p is a positive integer. 
 
   function out=cholx(x)
	   [R,p]=cholx();
	   if p==0
		   out=R;
	   else
		   out=real(sqrtm(x))';
     end 
  
=#

function cholx(x)
	try chol(x)
	catch ex
		println(ex)
	end
end

#A=[1 2 3; 4 5 6; 7 8 9]
#cholx(A)


