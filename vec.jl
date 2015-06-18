

" this is a function for VEC"

#=
orginal MATLAB code 

function x=vec(y)
	x=[];
	for i=1:cols(y)
		x=[x;y(:,i)];
	end
=#

# TEST
#x=[1,2,3,4];
#x=[1 2 3 4; 5 6 7 8]
#x=[]

function f(x)
	y=Array{Float64,1} 
	for i in [1:size(x, 2)]
		y=[y;x[:,i]]
	end
        return y;
end


# version 0.1, June 2015



