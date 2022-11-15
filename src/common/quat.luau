return {
	__mul = function(t,a,b)
		t[1]=a[1]*b[1]-a[2]*b[2]-a[3]*b[3]-a[4]*b[4]
		t[2]=a[2]*b[1]+a[1]*b[2]-a[4]*b[3]+a[3]*b[4]
		t[3]=a[3]*b[1]+a[4]*b[2]+a[1]*b[3]-a[2]*b[4]
		t[4]=a[4]*b[1]-a[3]*b[2]+a[2]*b[3]+a[1]*b[4]
	end,
	__pow = function(t,a,b)
		t[1]=a[1]*b[1]
		t[2]=a[2]*b[1]+a[1]*b[2]
		t[3]=a[3]*b[1]+a[1]*b[3]
		t[4]=a[4]*b[1]-a[3]*b[2]+a[2]*b[3]+a[1]*b[4]	
	end,
	__concat = function(t,a,b)
		t[1]=a[4]*b[1]-a[3]*b[2]+a[2]*b[3]+a[1]*b[4]
		t[2]=-a[4]*b[2]-a[2]*b[4]
		t[3]=-a[4]*b[3]-a[3]*b[4]
		t[4]=a[4]*b[4]
	end,
	__len = function(t) t[1],t[2],t[3],t[4]=t[4],t[3],-t[2],t[1] end,
	__unm = function(t) t[1],t[2],t[3],t[4]=-t[1],-t[2],-t[3],-t[4] end,
	rv = function(t) t[1],t[2],t[3],t[4]=t[1],t[2],t[3],-t[4] end,
	iv = function(t) t[1],t[2],t[3],t[4]=t[1],-t[2],-t[3],t[4] end,
	cj = function(t) t[1],t[2],t[3],t[4]=t[1],-t[2],-t[3],-t[4] end, -- aligns with the quaternion conjugate
	lc = function(t,a,b)
		t[1]=a[1]*b[1]-a[2]*b[2]-a[3]*b[3]-a[4]*b[4]
		t[2]=a[1]*b[2]+a[3]*b[4]
		t[3]=a[1]*b[3]-a[2]*b[4]
		t[4]=a[1]*b[4]
	end,
	-- less considerable inner products below
	fd = function(t,a,b)
		t[1]=a[1]*b[1]-a[2]*b[2]-a[3]*b[3]-a[4]*b[4]
		t[2]=a[2]*b[1]+a[1]*b[2]-a[4]*b[3]+a[3]*b[4]
		t[3]=a[3]*b[1]+a[4]*b[2]+a[1]*b[3]-a[2]*b[4]
		t[4]=a[4]*b[1]+a[1]*b[4]
	end,
	ip = function(t,a,b)
		t[1]=-a[2]*b[2]-a[3]*b[3]-a[4]*b[4]
		t[2]=-a[4]*b[3]+a[3]*b[4]
		t[3]=a[4]*b[2]-a[2]*b[4]
	end,
	sp = function(t,a,b)
		t[1]=a[1]*b[1]-a[2]*b[2]-a[3]*b[3]-a[4]*b[4]
	end
}
