
-- kgen ":3,,1;lgate&rgate=2"
-- also: direct exp from the paper "PGA Dynamics"

return {
	__mul = function(t,a,b)
t[1]=-a[9]*b[9]-a[10]*b[10]-a[11]*b[11]
t[6]=a[9]*b[7]+a[10]*b[8]-a[7]*b[9]-a[8]*b[10]
t[7]=-a[9]*b[6]+a[11]*b[8]+a[6]*b[9]-a[8]*b[11]
t[8]=-a[10]*b[6]-a[11]*b[7]+a[6]*b[10]+a[7]*b[11]
t[9]=a[11]*b[10]-a[10]*b[11]
t[10]=-a[11]*b[9]+a[9]*b[11]
t[11]=a[10]*b[9]-a[9]*b[10]
t[16]=a[11]*b[6]-a[10]*b[7]+a[9]*b[8]+a[8]*b[9]-a[7]*b[10]+a[6]*b[11]
	end,
	__pow = function(t,a,b)
t[16]=a[11]*b[6]-a[10]*b[7]+a[9]*b[8]+a[8]*b[9]-a[7]*b[10]+a[6]*b[11]
	end,
	__concat = function(t,a,b)
t[1]=a[11]*b[6]-a[10]*b[7]+a[9]*b[8]+a[8]*b[9]-a[7]*b[10]+a[6]*b[11]
	end,
	__len = function(t) t[6],t[7],t[8],t[9],t[10],t[11]=t[11],-t[10],t[9],t[8],-t[7],t[6] end,
	__unm = function(t) t[6],t[7],t[8],t[9],t[10],t[11]=-t[6],-t[7],-t[8],-t[9],-t[10],-t[11] end,
	rv = function(t) t[6],t[7],t[8],t[9],t[10],t[11]=-t[6],-t[7],-t[8],-t[9],-t[10],-t[11] end,
	iv = function(t) t[6],t[7],t[8],t[9],t[10],t[11]=t[6],t[7],t[8],t[9],t[10],t[11] end,
	cj = function(t) t[6],t[7],t[8],t[9],t[10],t[11]=-t[6],-t[7],-t[8],-t[9],-t[10],-t[11] end,
	lc = function(t,a,b)
t[1]=-a[9]*b[9]-a[10]*b[10]-a[11]*b[11]
	end,
	fd = function(t,a,b)
t[1]=-a[9]*b[9]-a[10]*b[10]-a[11]*b[11]
	end,
	ip = function(t,a,b)
t[1]=-a[9]*b[9]-a[10]*b[10]-a[11]*b[11]
	end,
	sp = function(t,a,b)
t[1]=-a[9]*b[9]-a[10]*b[10]-a[11]*b[11]
	end,
  exp = function(t, u)
    if not u then
      u = table.create(16, 0)
    end
    
    local B01, B02, B03, B12, B13, B23 = t[6], t[7], t[8], t[9], t[10], t[11]
    
    local l = B12*B13*B23
    if l == 0 then
      u[1], u[6], u[7], u[8], u[9], u[10], u[11], u[16] = 1, B01, B02, B03, 0, 0, 0, 0
    else
      local m, a = B01*B23+B02*B13+B03*B12, math.sqrt(l)
      local c, s = math.cos(a), math.sin(a)/a
      local v = m/l*(c-s)
      u[1], u[6], u[7], u[8], u[9], u[10], u[11], u[16] = c, s*B01 + v*B23, s*B02 + v*B13, s*B03 + v*B12, s*B12, s*B13, s*B23, m*s
    end
    
    return u
  end,
}