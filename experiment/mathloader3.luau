-- @densenet

-- mathloader
-- a string loader for stringified geometric algebra kernels
-- the input code needs to be conformal to the model below:
-- rules: <subexp> = <"a["<number>"]*b["<number>"]">
--				<exp> = <subexp> <<"+"|"-"><subexp>>*
-- given the rules above, the source can be expressed:
--[=[
<"t["<index>"]="<exp>>* | "\n"*
]=]


local n, x, l, r = bit32.bnot, bit32.bxor, bit32.lshift, bit32.rshift

local repl, band = bit32.replace, bit32.band

type u32 = number
type bitcode = { u32 }

local function kcompile(idxwidth: u32, src: string): bitcode
		-- 0: exit
	-- 1: t[.]= (assign)
	-- 2: +a[.]*b[.] (plusmulaccess)
	-- 3: -a[.]*b[.] (submulaccess)
	
	-- instructions are 2 bits
	-- indices are idxwidth bits, must be a multiple of two.
	
	if idxwidth % 2 ~= 0 then
		error("idxwidth must be an even number")
	end
	
	local program = {0}
	
	do
		local SIZE_IN = #src
		local ptr = 1
		local ins_ptr = 1
		local ins_depth = 0
		
		local function push(width: u32, chunk: u32)
			local fragment = 0
			repeat
				fragment = math.min(32 - ins_depth, width)
				
				program[ins_ptr] = repl(program[ins_ptr], band(chunk, n(l(-1, fragment))), ins_depth, fragment)
				
				chunk = r(chunk, fragment)
				
				ins_depth += fragment
				
				if ins_depth == 32 then
					ins_depth = 0
					ins_ptr += 1
					program[ins_ptr] = 0
				end
				
				width -= fragment
				
				if width < 0 then
					error("width underflow")
				end
				if ins_depth > 32 then
					error("ins_depth overflow")
				end
				if fragment <= 0 then
					error("fragment is not a positive number")
				end
			until width == 0
		end
		
		local byte = -1
		local function parse_index()
			local num = 0
			repeat
				byte = string.byte(src, ptr) - 48
				if byte > 9 or byte < 0 then
					error(string.format("nan: %s; ptr: %d", string.char(byte + 48), ptr))
				end
				num = num * 10 + byte
				ptr += 1
			until string.byte(src, ptr) == 93
			push(idxwidth, num - 1)
		end
		
		-- instruction loop	
		while ptr < SIZE_IN do
			byte = string.byte(src, ptr)
			if byte == 61 then
				if string.byte(src, ptr + 1) == 45 then
					byte = 45
					ptr += 1
				else
					byte = 43
				end
			end
			if byte == 10 then
				ptr += 1
			elseif byte == 116 then
				push(2, 1)
				ptr += 2
				parse_index()
				ptr += 1
			elseif byte == 43 then
				push(2, 2)
				ptr += 3
				parse_index()
				ptr += 4
				parse_index()
				ptr += 1
			elseif byte == 45 then
				push(2, 3)
				ptr += 3
				parse_index()
				ptr += 4
				parse_index()
				ptr += 1
			else
				error(string.format("syntax error: %x (%s)", byte, string.char(byte)))
			end
		end
	end
	
	return program
end

local function kload(idxwidth: u32, raw: bitcode): ({ number }, { number }, { number }) -> ()
	local index0: { u32 } = {}
	local index1: { u32 } = {}
	local index2: { u32 } = {}
	local index3: { u32 } = {}
	local index4: { u32 } = {}
	
	-- load
	do
		local RAWSIZE = #raw
		local ins_ptr = 1
		local ins_depth = 0
		
		local function pop_index(): u32
			local readw = idxwidth
	
			local idx = 0
			repeat
				local fragment = math.min(32 - ins_depth, readw)
		
				idx = repl(idx, bit32.extract(raw[ins_ptr], ins_depth, fragment), idxwidth - readw, fragment)
		
				readw -= fragment
		
				ins_depth += fragment
		
				if ins_depth == 32 then
					ins_ptr += 1
					ins_depth = 0
				end
		
			until readw == 0
	
			return idx + 1
		end
		
		while ins_ptr <= RAWSIZE do
			local ins = bit32.extract(raw[ins_ptr], ins_depth, 2)
			
			ins_depth += 2
			if ins_depth == 32 then
				ins_depth = 0
				ins_ptr += 1
			end
			
			local N = #index0

			if ins == 1 then
				index0[N + 1] = pop_index()
				index1[N + 1] = -1
				index2[N + 1] = -1
				index3[N + 1] = -1
				index4[N + 1] = -1
			elseif ins == 2 then
				if index0[N] == 0 then
					index0[N] = -2
					index3[N] = pop_index()
					index4[N] = pop_index()
				elseif index0[N] == -1 then
					index0[N] = -4
					index3[N] = pop_index()
					index4[N] = pop_index()
				else
					index0[N + 1] = 0
					index1[N + 1] = pop_index()
					index2[N + 1] = pop_index()
					index3[N + 1] = -1
					index4[N + 1] = -1
				end
			elseif ins == 3 then
				if index0[N] == 0 then
					index0[N] = -3
					index3[N] = pop_index()
					index4[N] = pop_index()
				elseif index0[N] == -1 then
					index0[N] = -5
					index3[N] = pop_index()
					index4[N] = pop_index()
				else
					index0[N + 1] = -1
					index1[N + 1] = pop_index()
					index2[N + 1] = pop_index()				
					index3[N + 1] = -1
					index4[N + 1] = -1
				end
			else
				break
			end
		end
	end
	
	-- procedure
	local DATASIZE = #index0
	return function(t,a,b)
		local t_idx = 0
		local acc = 0
		
		local i0, i1, i2, i3, i4 = 0, 0, 0, 0, 0
		
		for i = 1, DATASIZE do
			i0, i1, i2 = index0[i], index1[i], index2[i]
			
			if i0 == -2 then
				i3, i4 = index3[i], index4[i]
				acc=acc+a[i1]*b[i2]+a[i3]*b[i4]
			elseif i0 == -3 then
				i3, i4 = index3[i], index4[i]
				acc=acc+a[i1]*b[i2]-a[i3]*b[i4]
			elseif i0 == -4 then
				i3, i4 = index3[i], index4[i]
				acc=acc-a[i1]*b[i2]+a[i3]*b[i4]
			elseif i0 == -5 then
				i3, i4 = index3[i], index4[i]
				acc=acc-a[i1]*b[i2]-a[i3]*b[i4]
			elseif i0 == 0 then
				acc=acc+a[i1]*b[i2]
			elseif i0 == -1 then
				acc=acc-a[i1]*b[i2]
			else
				t[t_idx] = acc
				acc, t_idx = 0, i0
			end
		end
		if t_idx then
			t[t_idx] = acc
		end
	end
end

local function mload(src: string)
	-- instructions:
	-- t[.]=C -- 1
	-- t[.]=C*t[.] -- 2
	-- t[.]=t[.] -- 3
	-- t[.]=-t[.] -- 4
	local index0: { u32 } = {}
	local index1: { u32 } = {}
	local index2: { u32 } = {}
	local index3: { u32 } = {}
	
	do
		local ptr = 1
		local SRCSIZE = src:len()
		local N
		while ptr <= SRCSIZE do
			N = #index0
			local a, b, idx, neg = string.find(src, "t%[(%d+)%]=(%-?)", ptr)
			if a ~= nil then
				ptr = b+1
				local coef
				a, b, coef = string.find(src, "^(%d+)[*\n]", ptr)
				if coef ~= nil then
					ptr = b+1
				end
				local rhs_idx
				a, b, rhs_idx = string.find(src, "t%[(%d+)%]", ptr)
				if a ~= nil and string.sub(src, ptr-1, ptr-1) ~= "\n" then
					if coef~=nil then
						index0[N+1] = 2
						index1[N+1] = tonumber(idx) :: number
						index2[N+1] = (if neg == "" then 1 else -1) * tonumber(coef) :: number
						index3[N+1] = tonumber(rhs_idx) :: number
					else
						index0[N+1]=if neg=="" then 3 else 4
						index1[N+1]=tonumber(idx)::number
						index2[N+1]=math.huge
						index3[N+1]=tonumber(rhs_idx)::number
					end
					ptr = b+1
				else
					if coef==nil then
						coef=0
					end
					index0[N+1] = 1
					index1[N+1] = tonumber(idx) :: number
					index2[N+1] = (if neg == "" then 1 else -1) * tonumber(coef) :: number
					index3[N+1] = -1
				end
			else
				break
			end
		end
	end
	
	print(table.concat(index0, ","))
	print(table.concat(index1, ","))
	print(table.concat(index2, ","))
	print(table.concat(index3, ","))
	
	local PROGSIZE = #index0
	return function(t)
		local _t = table.create(#t)
		local i0, i1, i2, i3
		for i = 1, PROGSIZE do
			i0, i1, i2, i3 = index0[i], index1[i], index2[i], index3[i]
			if i0 == 3 then
				_t[i1] = t[i3]
			elseif i0 == 4 then
				_t[i1] = -t[i3]
			elseif i0 == 2 then
				_t[i1] = i2*t[i3]
			else
				_t[i1] = i2
			end
		end
		
		return _t
	end
end
--[===[

local pga3d_dual = mload[[
t[1]=t[16]
t[2]=t[15]
t[3]=-t[14]
t[4]=t[13]
t[5]=-t[12]
t[6]=t[11]
t[7]=-t[10]
t[8]=t[9]
t[9]=2*t[8]
t[10]=-t[7]
t[11]=t[6]
t[12]=t[5]
t[13]=-t[4]
t[14]=t[3]
t[15]=-t[2]
t[16]=t[1]
]]

local t = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16}

--print(table.concat(pga3d_dual(t), ","))

local stamp = os.clock()

for i = 1,1e5 do
	pga3d_dual(t)
end

print(os.clock()-stamp)



local pga3d_bitcode = compile(4, [[
t[1]=a[1]*b[1]+a[3]*b[3]+a[4]*b[4]+a[5]*b[5]-a[9]*b[9]-a[10]*b[10]-a[11]*b[11]-a[15]*b[15]
t[2]=a[2]*b[1]+a[1]*b[2]+a[6]*b[3]+a[7]*b[4]+a[8]*b[5]-a[3]*b[6]-a[4]*b[7]-a[5]*b[8]-a[12]*b[9]-a[13]*b[10]-a[14]*b[11]-a[9]*b[12]-a[10]*b[13]-a[11]*b[14]-a[16]*b[15]+a[15]*b[16]
t[3]=a[3]*b[1]+a[1]*b[3]+a[9]*b[4]+a[10]*b[5]-a[4]*b[9]-a[5]*b[10]-a[15]*b[11]-a[11]*b[15]
t[4]=a[4]*b[1]-a[9]*b[3]+a[1]*b[4]+a[11]*b[5]+a[3]*b[9]+a[15]*b[10]-a[5]*b[11]+a[10]*b[15]
t[5]=a[5]*b[1]-a[10]*b[3]-a[11]*b[4]+a[1]*b[5]-a[15]*b[9]+a[3]*b[10]+a[4]*b[11]-a[9]*b[15]
t[6]=a[6]*b[1]-a[3]*b[2]+a[2]*b[3]+a[12]*b[4]+a[13]*b[5]+a[1]*b[6]+a[9]*b[7]+a[10]*b[8]-a[7]*b[9]-a[8]*b[10]-a[16]*b[11]+a[4]*b[12]+a[5]*b[13]+a[15]*b[14]-a[14]*b[15]-a[11]*b[16]
t[7]=a[7]*b[1]-a[4]*b[2]-a[12]*b[3]+a[2]*b[4]+a[14]*b[5]-a[9]*b[6]+a[1]*b[7]+a[11]*b[8]+a[6]*b[9]+a[16]*b[10]-a[8]*b[11]-a[3]*b[12]-a[15]*b[13]+a[5]*b[14]+a[13]*b[15]+a[10]*b[16]
t[8]=a[8]*b[1]-a[5]*b[2]-a[13]*b[3]-a[14]*b[4]+a[2]*b[5]-a[10]*b[6]-a[11]*b[7]+a[1]*b[8]-a[16]*b[9]+a[6]*b[10]+a[7]*b[11]+a[15]*b[12]-a[3]*b[13]-a[4]*b[14]-a[12]*b[15]-a[9]*b[16]
t[9]=a[9]*b[1]-a[4]*b[3]+a[3]*b[4]+a[15]*b[5]+a[1]*b[9]+a[11]*b[10]-a[10]*b[11]+a[5]*b[15]
t[10]=a[10]*b[1]-a[5]*b[3]-a[15]*b[4]+a[3]*b[5]-a[11]*b[9]+a[1]*b[10]+a[9]*b[11]-a[4]*b[15]
t[11]=a[11]*b[1]+a[15]*b[3]-a[5]*b[4]+a[4]*b[5]+a[10]*b[9]-a[9]*b[10]+a[1]*b[11]+a[3]*b[15]
t[12]=a[12]*b[1]+a[9]*b[2]-a[7]*b[3]+a[6]*b[4]+a[16]*b[5]+a[4]*b[6]-a[3]*b[7]-a[15]*b[8]+a[2]*b[9]+a[14]*b[10]-a[13]*b[11]+a[1]*b[12]+a[11]*b[13]-a[10]*b[14]+a[8]*b[15]-a[5]*b[16]
t[13]=a[13]*b[1]+a[10]*b[2]-a[8]*b[3]-a[16]*b[4]+a[6]*b[5]+a[5]*b[6]+a[15]*b[7]-a[3]*b[8]-a[14]*b[9]+a[2]*b[10]+a[12]*b[11]-a[11]*b[12]+a[1]*b[13]+a[9]*b[14]-a[7]*b[15]+a[4]*b[16]
t[14]=a[14]*b[1]+a[11]*b[2]+a[16]*b[3]-a[8]*b[4]+a[7]*b[5]-a[15]*b[6]+a[5]*b[7]-a[4]*b[8]+a[13]*b[9]-a[12]*b[10]+a[2]*b[11]+a[10]*b[12]-a[9]*b[13]+a[1]*b[14]+a[6]*b[15]-a[3]*b[16]
t[15]=a[15]*b[1]+a[11]*b[3]-a[10]*b[4]+a[9]*b[5]+a[5]*b[9]-a[4]*b[10]+a[3]*b[11]+a[1]*b[15]
t[16]=a[16]*b[1]-a[15]*b[2]+a[14]*b[3]-a[13]*b[4]+a[12]*b[5]+a[11]*b[6]-a[10]*b[7]+a[9]*b[8]+a[8]*b[9]-a[7]*b[10]+a[6]*b[11]-a[5]*b[12]+a[4]*b[13]-a[3]*b[14]+a[2]*b[15]+a[1]*b[16]
]])

for i, v in pga3d_bitcode do
	pga3d_bitcode[i] = string.format("0x%x", v)
end
print(table.concat(pga3d_bitcode, ","))


--[=[
0x388a0081,0x6788d123,0x617bbaae,0x68961080,0x8f52d1e3,
0x9ce2f74d,0xf27b8eb7,0x9febbfda,0x38882028,0xe5383d26,
0xe37abae,0x4a8c228c,0xa939ea0a,0xa70491e9,0x8ed023ac,
0x7a3a3a4a,0x8612c165,0x509323b8,0xe1b799a2,0x2b3abf97,
0xafb7deb1,0xc4f0699f,0x34d8c62b,0x59ea6096,0xcba7e7e8,
0xecb52cee,0x4c1e77e6,0x63dcb31,0x709ab59d,0xa9a95a3f,
0xfd3f0bbe,0x8a1f8fa,0x93a3288f,0x2a9e6a80,0x24c26979,
0xa2b428fb,0x9e3ea290,0x34cba0aa,0xa638990e,0x2eb78aa0,
0x3589b188,0xd8b5393e,0x39da067e,0x9f2ab0ab,0xb1f4f9ed,
0xf27c660c,0xe952458f,0x468ddcb6,0xc0aebaba,0x7cee6f62,
0xbe1a836d,0x5ed1a37c,0xe3273d92,0x3b9a869b,0x2f96d0b2,
0xcaa0eb9f,0xf8492239,0xf782a2a4,0xcb61ec3e,0x75a92e3c,
0x6e1e789a,0xeb4e969,0xf0b86d2f,0x0

]=]

local kpga3d = loadkernel(4, { 
0x388a0081,0x6788d123,0x617bbaae,0x68961080,0x8f52d1e3,
0x9ce2f74d,0xf27b8eb7,0x9febbfda,0x38882028,0xe5383d26,
0xe37abae,0x4a8c228c,0xa939ea0a,0xa70491e9,0x8ed023ac,
0x7a3a3a4a,0x8612c165,0x509323b8,0xe1b799a2,0x2b3abf97,
0xafb7deb1,0xc4f0699f,0x34d8c62b,0x59ea6096,0xcba7e7e8,
0xecb52cee,0x4c1e77e6,0x63dcb31,0x709ab59d,0xa9a95a3f,
0xfd3f0bbe,0x8a1f8fa,0x93a3288f,0x2a9e6a80,0x24c26979,
0xa2b428fb,0x9e3ea290,0x34cba0aa,0xa638990e,0x2eb78aa0,
0x3589b188,0xd8b5393e,0x39da067e,0x9f2ab0ab,0xb1f4f9ed,
0xf27c660c,0xe952458f,0x468ddcb6,0xc0aebaba,0x7cee6f62,
0xbe1a836d,0x5ed1a37c,0xe3273d92,0x3b9a869b,0x2f96d0b2,
0xcaa0eb9f,0xf8492239,0xf782a2a4,0xcb61ec3e,0x75a92e3c,
0x6e1e789a,0xeb4e969,0xf0b86d2f,0x0})


local __TIMING = false

if __TIMING then
	local t = table.create(16, 0)
	local a = table.create(16, 0)
	local b = table.create(16, 0)

	a[2],b[7] = 1, -1
	
	local stamp = os.clock()
	
	for i = 1,1e4 do
		kpga3d(t, a, b)
	end
	
	print(string.format("--%.3fms", 1e3*(os.clock()-stamp)))
	
	for i = 1, 16 do
		t[i] = string.format("%.2f", t[i])
	end
	
	print("--", table.concat(t, ","))	
end


local kpga3d_l = function(t,a,b)
t[1]=a[1]*b[1]+a[3]*b[3]+a[4]*b[4]+a[5]*b[5]-a[9]*b[9]-a[10]*b[10]-a[11]*b[11]-a[15]*b[15]
t[2]=a[2]*b[1]+a[1]*b[2]+a[6]*b[3]+a[7]*b[4]+a[8]*b[5]-a[3]*b[6]-a[4]*b[7]-a[5]*b[8]-a[12]*b[9]-a[13]*b[10]-a[14]*b[11]-a[9]*b[12]-a[10]*b[13]-a[11]*b[14]-a[16]*b[15]+a[15]*b[16]
t[3]=a[3]*b[1]+a[1]*b[3]+a[9]*b[4]+a[10]*b[5]-a[4]*b[9]-a[5]*b[10]-a[15]*b[11]-a[11]*b[15]
t[4]=a[4]*b[1]-a[9]*b[3]+a[1]*b[4]+a[11]*b[5]+a[3]*b[9]+a[15]*b[10]-a[5]*b[11]+a[10]*b[15]
t[5]=a[5]*b[1]-a[10]*b[3]-a[11]*b[4]+a[1]*b[5]-a[15]*b[9]+a[3]*b[10]+a[4]*b[11]-a[9]*b[15]
t[6]=a[6]*b[1]-a[3]*b[2]+a[2]*b[3]+a[12]*b[4]+a[13]*b[5]+a[1]*b[6]+a[9]*b[7]+a[10]*b[8]-a[7]*b[9]-a[8]*b[10]-a[16]*b[11]+a[4]*b[12]+a[5]*b[13]+a[15]*b[14]-a[14]*b[15]-a[11]*b[16]
t[7]=a[7]*b[1]-a[4]*b[2]-a[12]*b[3]+a[2]*b[4]+a[14]*b[5]-a[9]*b[6]+a[1]*b[7]+a[11]*b[8]+a[6]*b[9]+a[16]*b[10]-a[8]*b[11]-a[3]*b[12]-a[15]*b[13]+a[5]*b[14]+a[13]*b[15]+a[10]*b[16]
t[8]=a[8]*b[1]-a[5]*b[2]-a[13]*b[3]-a[14]*b[4]+a[2]*b[5]-a[10]*b[6]-a[11]*b[7]+a[1]*b[8]-a[16]*b[9]+a[6]*b[10]+a[7]*b[11]+a[15]*b[12]-a[3]*b[13]-a[4]*b[14]-a[12]*b[15]-a[9]*b[16]
t[9]=a[9]*b[1]-a[4]*b[3]+a[3]*b[4]+a[15]*b[5]+a[1]*b[9]+a[11]*b[10]-a[10]*b[11]+a[5]*b[15]
t[10]=a[10]*b[1]-a[5]*b[3]-a[15]*b[4]+a[3]*b[5]-a[11]*b[9]+a[1]*b[10]+a[9]*b[11]-a[4]*b[15]
t[11]=a[11]*b[1]+a[15]*b[3]-a[5]*b[4]+a[4]*b[5]+a[10]*b[9]-a[9]*b[10]+a[1]*b[11]+a[3]*b[15]
t[12]=a[12]*b[1]+a[9]*b[2]-a[7]*b[3]+a[6]*b[4]+a[16]*b[5]+a[4]*b[6]-a[3]*b[7]-a[15]*b[8]+a[2]*b[9]+a[14]*b[10]-a[13]*b[11]+a[1]*b[12]+a[11]*b[13]-a[10]*b[14]+a[8]*b[15]-a[5]*b[16]
t[13]=a[13]*b[1]+a[10]*b[2]-a[8]*b[3]-a[16]*b[4]+a[6]*b[5]+a[5]*b[6]+a[15]*b[7]-a[3]*b[8]-a[14]*b[9]+a[2]*b[10]+a[12]*b[11]-a[11]*b[12]+a[1]*b[13]+a[9]*b[14]-a[7]*b[15]+a[4]*b[16]
t[14]=a[14]*b[1]+a[11]*b[2]+a[16]*b[3]-a[8]*b[4]+a[7]*b[5]-a[15]*b[6]+a[5]*b[7]-a[4]*b[8]+a[13]*b[9]-a[12]*b[10]+a[2]*b[11]+a[10]*b[12]-a[9]*b[13]+a[1]*b[14]+a[6]*b[15]-a[3]*b[16]
t[15]=a[15]*b[1]+a[11]*b[3]-a[10]*b[4]+a[9]*b[5]+a[5]*b[9]-a[4]*b[10]+a[3]*b[11]+a[1]*b[15]
t[16]=a[16]*b[1]-a[15]*b[2]+a[14]*b[3]-a[13]*b[4]+a[12]*b[5]+a[11]*b[6]-a[10]*b[7]+a[9]*b[8]+a[8]*b[9]-a[7]*b[10]+a[6]*b[11]-a[5]*b[12]+a[4]*b[13]-a[3]*b[14]+a[2]*b[15]+a[1]*b[16]
end

local __TIMING2 = false

if __TIMING2 then
	local t = table.create(16, 0)
	local a = table.create(16, 0)
	local b = table.create(16, 0)
	
	for i = 2, 5 do
		a[i], b[i] = math.random(), math.random()
	end
	
	local stamp = os.clock()
	
	for i = 1,1e5 do
		kpga3d_l(t, a, b)
	end
	
	print(string.format("--%.3fms", 1e3*(os.clock()-stamp)))
end

local __TIMING_DENSE = false

if __TIMING_DENSE then
	local t = table.create(16, 0)
	local a = table.create(16, 0)
	local b = table.create(16, 0)
	
	for i = 1, 16 do
		a[i], b[i] = math.random(1000) - 500, math.random(1000) - 500
	end
	
	local stamp = os.clock()
	
	for i = 1,1e4 do
		kpga3d(t, a, b)
	end
	
	print(string.format("--%.3fms", 1e3*(os.clock()-stamp)))
	
	for i = 1, 16 do
		t[i] = string.format("%.2f", t[i])
	end
	
	print("--", table.concat(t, ","))	
end

local pga3d_regprod_bitcode = compile(4, [[
t[1]=a[16]*b[1]-a[15]*b[2]+a[14]*b[3]-a[13]*b[4]+a[12]*b[5]+a[11]*b[6]-a[10]*b[7]+a[9]*b[8]+a[8]*b[9]-a[7]*b[10]+a[6]*b[11]-a[5]*b[12]+a[4]*b[13]-a[3]*b[14]+a[2]*b[15]+a[1]*b[16]
t[2]=-a[16]*b[2]-a[14]*b[6]+a[13]*b[7]-a[12]*b[8]-a[8]*b[12]+a[7]*b[13]-a[6]*b[14]-a[2]*b[16]
t[3]=-a[16]*b[3]-a[15]*b[6]+a[13]*b[9]-a[12]*b[10]-a[10]*b[12]+a[9]*b[13]-a[6]*b[15]-a[3]*b[16]
t[4]=-a[16]*b[4]-a[15]*b[7]+a[14]*b[9]-a[12]*b[11]-a[11]*b[12]+a[9]*b[14]-a[7]*b[15]-a[4]*b[16]
t[5]=-a[16]*b[5]-a[15]*b[8]+a[14]*b[10]-a[13]*b[11]-a[11]*b[13]+a[10]*b[14]-a[8]*b[15]-a[5]*b[16]
t[6]=a[16]*b[6]-a[13]*b[12]+a[12]*b[13]+a[6]*b[16]
t[7]=a[16]*b[7]-a[14]*b[12]+a[12]*b[14]+a[7]*b[16]
t[8]=a[16]*b[8]-a[14]*b[13]+a[13]*b[14]+a[8]*b[16]
t[9]=a[16]*b[9]-a[15]*b[12]+a[12]*b[15]+a[9]*b[16]
t[10]=a[16]*b[10]-a[15]*b[13]+a[13]*b[15]+a[10]*b[16]
t[11]=a[16]*b[11]-a[15]*b[14]+a[14]*b[15]+a[11]*b[16]
t[12]=-a[16]*b[12]-a[12]*b[16]
t[13]=-a[16]*b[13]-a[13]*b[16]
t[14]=-a[16]*b[14]-a[14]*b[16]
t[15]=-a[16]*b[15]-a[15]*b[16]
t[16]=a[16]*b[16]
]])


for i, v in pga3d_regprod_bitcode do
	pga3d_regprod_bitcode[i] = string.format("0x%x", v)
end
print(table.concat(pga3d_regprod_bitcode, ","))

--[=[
0xd87b0f81,0x6a4b8f32,0x879e269d,0xad3a5a5b,0x2e1b4bc3,0x5dc7f17c,0xadf7bdb2,0x9f1f57c6,0x8c97b2fc,0xf22b9e6f,0xff37cbe5,0xabe366ec,0xf9bd8aeb,0xfb4fd1f3,0xcaeb39d9,0x7d3e7f66,0x2ebcd7e5,0x76f99f5b,0x77dadbaf,0xb72cddfe,0xfb8fa1f7,0xe97e2eba,0x9bb2cee7,0xb7bafa9f,0xffb7eaed,0x3cff1fbe,0xfdf7fd7f,0xf7fbeff9,0x3fe

]=]

]===]


local __RELEASE = false
if __RELEASE then
	return {
		kcompile = kcompile,
		kload = kload,
		mload = mload
	}
end
