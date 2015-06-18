

function getqr(a)

	(q,r)=qr(a)
	for i in [1:size(a,1)]
		if r[i,i]<0
			q[:,i]=-q[:,i]
		end
	end
	return q
end


# Test
#a=[1 2 3; 4 5 6; 7 8 9]
#getqr(a)


