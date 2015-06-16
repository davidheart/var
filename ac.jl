


function ac(data)
        y=data
	x=lag0(data,1)
	x=[x ones(size(data,1),1)]
	y=y[2:end]
	x=x[2:end,:]
	b=x\y
	return b[1]
end 


# test 
#x=[1;2;3;4;5;6;7]
#ac(x)

