

function getqr(a)

	(q,r)=qr(a)
	for i in [i:size(a,1)]
		if r[i,i]<0
			q[:,i]=-q[:,i]
		end
	end
	return q
end


