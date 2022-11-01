
local n, x, l, r = bit32.bnot, bit32.bxor, bit32.lshift, bit32.rshift

local repl, band = bit32.replace, bit32.band

local function compile(src: string)
	-- 0: exit
	-- 1: t[.]=a[.]*b[.] (assignmulaccess)
	-- 2: +a[.]*b[.] (plusmulaccess)
	-- 3: -a[.]*b[.] (submulaccess)
	
	-- instructions are 2 bits
	-- indices are idxwidth bits
	
	local program = {0}
	local data = {}
	
	do
		local SIZE_IN = #src
		local ptr = 1
		local ins_ptr = 1
		local ins_depth = 0
		
		local function push(width: number, chunk: number)
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
		
		local byte	
		local function parse_index()
			local num = 0
			repeat
				byte = string.byte(src, ptr) - 48
				if byte > 9 or byte < 0 then
					error(string.format("nan: %s; ptr: %d", string.char(byte + 48), ptr))
				end
				--print("digit: ", byte)
				num = num * 10 + byte
				ptr += 1
			until string.byte(src, ptr) == 93

			data[#data + 1] = num
		end
		
		-- instruction loop	
		while ptr < SIZE_IN do
			byte = string.byte(src, ptr)
			if byte == 116 then
				push(2, 1)
				ptr += 2
				parse_index()
				ptr += 4
				parse_index()
				ptr += 4
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
			elseif byte == 10 then
				ptr += 1
			else
				error(string.format("syntax error: %x (%s)", byte, string.char(byte)))
			end
		end
	end
	
	return program, data
end

local function loadkernel(program: { number }, data: { number }): ({ number }, { number }, { number }) -> ()
	local PROGSIZE32 = #program
	return function(t,a,b)
		local ins_ptr = 1
		local ins_depth = 0
		local data_ptr = 1
		
		local ins = 0
		
		local t_idx = 0
		
		local acc = 0
		
		while ins_ptr <= PROGSIZE32 do
			ins = bit32.extract(program[ins_ptr], ins_depth, 2)
			
			ins_depth += 2
			if ins_depth == 32 then
				ins_depth = 0
				ins_ptr += 1
			end
				
			if ins == 1 then
				if t_idx > 0 then
					t[t_idx] = acc
					acc = 0
				end
				
				--print("T")
				
				t_idx = data[data_ptr]
				data_ptr += 1
				local a_idx = data[data_ptr]
				data_ptr += 1
				local b_idx = data[data_ptr]
				data_ptr += 1
				
				if a[a_idx] == 0 or b[b_idx] == 0 then
					continue
				end
				
				--t_idx = pop_index()
				
				--print(t_idx, a_idx, b_idx)
				
				acc+=a[a_idx]*b[b_idx]
			elseif ins == 2 then
				--print("+")
		
				local a_idx = data[data_ptr]
				data_ptr += 1
				local b_idx = data[data_ptr]
				data_ptr += 1
				--print(t_idx, a_idx, b_idx)
				
				if a[a_idx] == 0 or b[b_idx] == 0 then
					continue
				end
				
				acc+=a[a_idx]*b[b_idx]
			elseif ins == 3 then
				--print("-")
				
				local a_idx = data[data_ptr]
				data_ptr += 1
				local b_idx = data[data_ptr]
				data_ptr += 1
				
				if a[a_idx] == 0 or b[b_idx] == 0 then
					continue
				end				
			
				--print(t_idx, a_idx, b_idx)
				
				acc-=a[a_idx]*b[b_idx]
			else
				--print(string.format("eof? ins: %d; ins_ptr: %d; ins_depth: %d", ins, ins_ptr, ins_depth))
				break
			end
		end
		
		if t_idx then
			t[t_idx] = acc
		end
	end
end

__PRINT_PGA3D_BC = false
if __PRINT_PGA3D_BC then
	local pga3d_program, pga3d_data = compile[=[
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
]=]
	
	local out = {}
	for i = 1, #pga3d_program do
		out[i] = string.format("0x%x", pga3d_program[i])
	end
	
	print("-- program: ", table.concat(out, ","))
	print("-- data: ", table.concat(pga3d_data, ","))
end

__ARMCHAIR = true

if __ARMCHAIR then
	local kpga3d = loadkernel(
		{ 0xfea9ffa9,0xffa9bfff,0xebbdbaad,0xfabfaaad,0xabfaaebd,0xffabbefd,0xebbdbaad,0xfab9aeb9,0xeaf9eeba,
			0xeee9baeb,0xaeb9ebae,0xaeeebaed,0x0 },
		{ 1,1,1,3,3,4,4,5,5,9,9,10,10,11,11,15,15,2,2,1,1,2,6,3,7,4,8,5,3,6,4,7,5,8,12,9,13,10,14,11,9,12,10,
			13,11,14,16,15,15,16,3,3,1,1,3,9,4,10,5,4,9,5,10,15,11,11,15,4,4,1,9,3,1,4,11,5,3,9,15,10,5,11,10,
			15,5,5,1,10,3,11,4,1,5,15,9,3,10,4,11,9,15,6,6,1,3,2,2,3,12,4,13,5,1,6,9,7,10,8,7,9,8,10,16,11,4,12
			,5,13,15,14,14,15,11,16,7,7,1,4,2,12,3,2,4,14,5,9,6,1,7,11,8,6,9,16,10,8,11,3,12,15,13,5,14,13,15,10
			,16,8,8,1,5,2,13,3,14,4,2,5,10,6,11,7,1,8,16,9,6,10,7,11,15,12,3,13,4,14,12,15,9,16,9,9,1,4,3,3,4,15,
			5,1,9,11,10,10,11,5,15,10,10,1,5,3,15,4,3,5,11,9,1,10,9,11,4,15,11,11,1,15,3,5,4,4,5,10,9,9,10,1,11,
			3,15,12,12,1,9,2,7,3,6,4,16,5,4,6,3,7,15,8,2,9,14,10,13,11,1,12,11,13,10,14,8,15,5,16,13,13,1,10,2,8,
			3,16,4,6,5,5,6,15,7,3,8,14,9,2,10,12,11,11,12,1,13,9,14,7,15,4,16,14,14,1,11,2,16,3,8,4,7,5,15,6,5,7,
			4,8,13,9,12,10,2,11,10,12,9,13,1,14,6,15,3,16,15,15,1,11,3,10,4,9,5,5,9,4,10,3,11,1,15,16,16,1,15,2,14,
			3,13,4,12,5,11,6,10,7,9,8,8,9,7,10,6,11,5,12,4,13,3,14,2,15,1,16 })
	
		local kpga3d_l = (loadstring([[
return function(t,a,b)
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
]]) :: any)() :: (number, number, number) -> ()
	
	local a = table.create(16, 0)
	local b = table.create(16, 0)
	local t = table.create(16)
	
	for i = 2, 5 do
		a[i] = math.random()
		b[i] = math.random()
	end
	
	local stamp = os.clock()
	
	for i = 1, 1e4 do
		kpga3d_l(t, a, b)
	end
	
	print(string.format("-- kpga3d_l: %.2fms", 1e3* (os.clock() - stamp)))
	
	local stamp = os.clock()
	
	for i = 1, 1e4 do
		kpga3d(t, a, b)
	end
	
	print(string.format("-- kpga3d: %.2fms", 1e3* (os.clock() - stamp)))
end

-- program: 	0xfea9ffa9,0xffa9bfff,0xebbdbaad,0xfabfaaad,0xabfaaebd,0xffabbefd,0xebbdbaad,0xfab9aeb9,0xeaf9eeba,0xeee9baeb,0xaeb9ebae,0xaeeebaed,0x0
-- data: 	1,1,1,3,3,4,4,5,5,9,9,10,10,11,11,15,15,2,2,1,1,2,6,3,7,4,8,5,3,6,4,7,5,8,12,9,13,10,14,11,9,12,10,13,11,14,16,15,15,16,3,3,1,1,3,9,4,10,5,4,9,5,10,15,11,11,15,4,4,1,9,3,1,4,11,5,3,9,15,10,5,11,10,15,5,5,1,10,3,11,4,1,5,15,9,3,10,4,11,9,15,6,6,1,3,2,2,3,12,4,13,5,1,6,9,7,10,8,7,9,8,10,16,11,4,12,5,13,15,14,14,15,11,16,7,7,1,4,2,12,3,2,4,14,5,9,6,1,7,11,8,6,9,16,10,8,11,3,12,15,13,5,14,13,15,10,16,8,8,1,5,2,13,3,14,4,2,5,10,6,11,7,1,8,16,9,6,10,7,11,15,12,3,13,4,14,12,15,9,16,9,9,1,4,3,3,4,15,5,1,9,11,10,10,11,5,15,10,10,1,5,3,15,4,3,5,11,9,1,10,9,11,4,15,11,11,1,15,3,5,4,4,5,10,9,9,10,1,11,3,15,12,12,1,9,2,7,3,6,4,16,5,4,6,3,7,15,8,2,9,14,10,13,11,1,12,11,13,10,14,8,15,5,16,13,13,1,10,2,8,3,16,4,6,5,5,6,15,7,3,8,14,9,2,10,12,11,11,12,1,13,9,14,7,15,4,16,14,14,1,11,2,16,3,8,4,7,5,15,6,5,7,4,8,13,9,12,10,2,11,10,12,9,13,1,14,6,15,3,16,15,15,1,11,3,10,4,9,5,5,9,4,10,3,11,1,15,16,16,1,15,2,14,3,13,4,12,5,11,6,10,7,9,8,8,9,7,10,6,11,5,12,4,13,3,14,2,15,1,16

-- without zero guard
-- kpga3d_l: 21.57ms
-- kpga3d: 125.85ms

-- with
-- kpga3d_l: 21.55ms
-- kpga3d: 113.59ms
