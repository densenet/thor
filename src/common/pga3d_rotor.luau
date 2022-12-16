
-- pga3d_rotor := kgen "3,,1;lgate&rgate=0,2,4"
-- contains natural log, which reverses exp from pga3d_bivector

return {
	__mul = function(t,a,b)
t[1]=a[1]*b[1]-a[9]*b[9]-a[10]*b[10]-a[11]*b[11]
t[6]=a[6]*b[1]+a[1]*b[6]+a[9]*b[7]+a[10]*b[8]-a[7]*b[9]-a[8]*b[10]-a[16]*b[11]-a[11]*b[16]
t[7]=a[7]*b[1]-a[9]*b[6]+a[1]*b[7]+a[11]*b[8]+a[6]*b[9]+a[16]*b[10]-a[8]*b[11]+a[10]*b[16]
t[8]=a[8]*b[1]-a[10]*b[6]-a[11]*b[7]+a[1]*b[8]-a[16]*b[9]+a[6]*b[10]+a[7]*b[11]-a[9]*b[16]
t[9]=a[9]*b[1]+a[1]*b[9]+a[11]*b[10]-a[10]*b[11]
t[10]=a[10]*b[1]-a[11]*b[9]+a[1]*b[10]+a[9]*b[11]
t[11]=a[11]*b[1]+a[10]*b[9]-a[9]*b[10]+a[1]*b[11]
t[16]=a[16]*b[1]+a[11]*b[6]-a[10]*b[7]+a[9]*b[8]+a[8]*b[9]-a[7]*b[10]+a[6]*b[11]+a[1]*b[16]
	end,
	__pow = function(t,a,b)
t[1]=a[1]*b[1]
t[6]=a[6]*b[1]+a[1]*b[6]
t[7]=a[7]*b[1]+a[1]*b[7]
t[8]=a[8]*b[1]+a[1]*b[8]
t[9]=a[9]*b[1]+a[1]*b[9]
t[10]=a[10]*b[1]+a[1]*b[10]
t[11]=a[11]*b[1]+a[1]*b[11]
t[16]=a[16]*b[1]+a[11]*b[6]-a[10]*b[7]+a[9]*b[8]+a[8]*b[9]-a[7]*b[10]+a[6]*b[11]+a[1]*b[16]
	end,
	__concat = function(t,a,b)
t[1]=a[16]*b[1]+a[11]*b[6]-a[10]*b[7]+a[9]*b[8]+a[8]*b[9]-a[7]*b[10]+a[6]*b[11]+a[1]*b[16]
t[6]=a[16]*b[6]+a[6]*b[16]
t[7]=a[16]*b[7]+a[7]*b[16]
t[8]=a[16]*b[8]+a[8]*b[16]
t[9]=a[16]*b[9]+a[9]*b[16]
t[10]=a[16]*b[10]+a[10]*b[16]
t[11]=a[16]*b[11]+a[11]*b[16]
t[16]=a[16]*b[16]
	end,
	__len = function(t) t[1],t[6],t[7],t[8],t[9],t[10],t[11],t[16]=t[16],t[11],-t[10],t[9],t[8],-t[7],t[6],t[1] end,
	__unm = function(t) t[1],t[6],t[7],t[8],t[9],t[10],t[11],t[16]=-t[1],-t[6],-t[7],-t[8],-t[9],-t[10],-t[11],-t[16] end,
	rv = function(t) t[1],t[6],t[7],t[8],t[9],t[10],t[11],t[16]=t[1],-t[6],-t[7],-t[8],-t[9],-t[10],-t[11],t[16] end,
	iv = function(t) t[1],t[6],t[7],t[8],t[9],t[10],t[11],t[16]=t[1],t[6],t[7],t[8],t[9],t[10],t[11],t[16] end,
	cj = function(t) t[1],t[6],t[7],t[8],t[9],t[10],t[11],t[16]=t[1],-t[6],-t[7],-t[8],-t[9],-t[10],-t[11],t[16] end,
	lc = function(t,a,b)
t[1]=a[1]*b[1]-a[9]*b[9]-a[10]*b[10]-a[11]*b[11]
t[6]=a[1]*b[6]-a[11]*b[16]
t[7]=a[1]*b[7]+a[10]*b[16]
t[8]=a[1]*b[8]-a[9]*b[16]
t[9]=a[1]*b[9]
t[10]=a[1]*b[10]
t[11]=a[1]*b[11]
t[16]=a[1]*b[16]
	end,
	fd = function(t,a,b)
t[1]=a[1]*b[1]-a[9]*b[9]-a[10]*b[10]-a[11]*b[11]
t[6]=a[6]*b[1]+a[1]*b[6]-a[16]*b[11]-a[11]*b[16]
t[7]=a[7]*b[1]+a[1]*b[7]+a[16]*b[10]+a[10]*b[16]
t[8]=a[8]*b[1]+a[1]*b[8]-a[16]*b[9]-a[9]*b[16]
t[9]=a[9]*b[1]+a[1]*b[9]
t[10]=a[10]*b[1]+a[1]*b[10]
t[11]=a[11]*b[1]+a[1]*b[11]
t[16]=a[16]*b[1]+a[1]*b[16]
	end,
	ip = function(t,a,b)
t[1]=-a[9]*b[9]-a[10]*b[10]-a[11]*b[11]
t[6]=-a[16]*b[11]-a[11]*b[16]
t[7]=a[16]*b[10]+a[10]*b[16]
t[8]=-a[16]*b[9]-a[9]*b[16]
	end,
	sp = function(t,a,b)
t[1]=a[1]*b[1]-a[9]*b[9]-a[10]*b[10]-a[11]*b[11]
	end,
	log = function(t,u)
		-- result is a bivector
		if not u then
			u = table.create(16, 0)
		end
		
		if t[1] == 0 then
			table.move(t, 6, 8, 6, u)
		else
			local ALPHA, B12, B13, B23, BETA = t[1], t[9], t[10], t[11], t[16]
			local a = 1/(1-ALPHA*ALPHA)
			local b = math.acos(ALPHA)*math.sqrt(a)
			local c = a*BETA*(1-ALPHA*b)
			u[6],u[7],u[8],u[9],u[10],u[11]=c*B23+b*t[6],c*B13+b*t[7],c*B12+b*t[8],b*B12,b*B13,b*B23
		end
		
		return u
	end
}
