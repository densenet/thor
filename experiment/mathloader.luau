-- seeing how concisely and efficiently i can carry out stringified math kernels in luau
-- the input code needs to be conformal to the model below:
-- rules: <subexp> = <"a["<number>"]*b["<number>"]">
--				<exp> = <subexp> <<"+"|"-"><subexp>>*
--				<index> = <num>
-- and also needs to conform to the semantic rule:
--				<index> = <#<index> + 1>
-- where # translates to a counting operator on amount of index tokens
-- given the rules above, the source can be expressed:
--[=[
<"t["<index>"]="<expansion>>*
]=]

-- local bit2str = require("src\\util\\bit2str")

local n, x, l, r = bit32.bnot, bit32.bxor, bit32.lshift, bit32.rshift

local repl, band = bit32.replace, bit32.band

local function compile(idxwidth: number, src: string)
	-- 0: exit
	-- 1: t[.]=a[.]*b[.] (assignmulaccess)
	-- 2: +a[.]*b[.] (plusmulaccess)
	-- 3: -a[.]*b[.] (submulaccess)
	
	-- instructions are 2 bits
	-- indices are idxwidth bits
	
	if idxwidth % 2 ~= 0 then
		error("idxwidth must be an even number")
	end
	
	local program = {0}
	
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
			push(idxwidth, num - 1)
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
	
	return program
end

local function loadkernel(idxwidth: number, program: { number }): ({ number }, { number }, { number }) -> ()
	local PROGSIZE32 = #program
	return function(t,a,b)
		local ins_ptr = 1
		local ins_depth = 0
		
		--[=[
		local function pop_index(): number
			local readw = idxwidth
		
			local idx = 0
			repeat
				local fragment = math.min(32 - ins_depth, readw)
			
				idx = repl(idx, bit32.extract(program[ins_ptr], ins_depth, fragment), idxwidth - readw, fragment)
			
				readw -= fragment
			
				ins_depth += fragment
			
				if ins_depth == 32 then
					ins_ptr += 1
					ins_depth = 0
				end
			
			until readw == 0
		
			return idx + 1
		end
	]=]
		
		local ins = 0
		
		local t_idx = 0
		
		local acc = 0
		
		while ins_ptr <= PROGSIZE32 do
			ins = bit32.extract(program[ins_ptr], ins_depth, 2)
			local a_idx, b_idx = 0, 0
			
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
				
				--t_idx = pop_index()
				
				do
					local readw = idxwidth
		
					t_idx = 0
					repeat
						local fragment = math.min(32 - ins_depth, readw)
			
						t_idx = repl(t_idx, bit32.extract(program[ins_ptr], ins_depth, fragment), idxwidth - readw, fragment)
			
						readw -= fragment
			
						ins_depth += fragment
			
						if ins_depth == 32 then
							ins_ptr += 1
							ins_depth = 0
						end
					until readw == 0
		
					t_idx += 1
				end
				
				do
					local readw = idxwidth
		
					repeat
						local fragment = math.min(32 - ins_depth, readw)
			
						a_idx = repl(a_idx, bit32.extract(program[ins_ptr], ins_depth, fragment), idxwidth - readw, fragment)
			
						readw -= fragment
			
						ins_depth += fragment
			
						if ins_depth == 32 then
							ins_ptr += 1
							ins_depth = 0
						end
					until readw == 0
		
					a_idx += 1
				end
				
				if a[a_idx] == 0 then
					continue
				end
				
				do
					local readw = idxwidth
		
					repeat
						local fragment = math.min(32 - ins_depth, readw)
			
						b_idx = repl(b_idx, bit32.extract(program[ins_ptr], ins_depth, fragment), idxwidth - readw, fragment)
			
						readw -= fragment
			
						ins_depth += fragment
			
						if ins_depth == 32 then
							ins_ptr += 1
							ins_depth = 0
						end
					until readw == 0
		
					b_idx += 1
				end
				
				if b[b_idx] == 0 then
					continue
				end
	
				--print(t_idx, a_idx, b_idx)
				
				acc+=a[a_idx]*b[b_idx]
			elseif ins == 2 then
				--print("+")
				
				do
					local readw = idxwidth
		
					repeat
						local fragment = math.min(32 - ins_depth, readw)
			
						a_idx = repl(a_idx, bit32.extract(program[ins_ptr], ins_depth, fragment), idxwidth - readw, fragment)
			
						readw -= fragment
			
						ins_depth += fragment
			
						if ins_depth == 32 then
							ins_ptr += 1
							ins_depth = 0
						end
					until readw == 0
		
					a_idx += 1
				end
				
				if a[a_idx] == 0 then
					continue
				end
				
				do
					local readw = idxwidth
		
					repeat
						local fragment = math.min(32 - ins_depth, readw)
			
						b_idx = repl(b_idx, bit32.extract(program[ins_ptr], ins_depth, fragment), idxwidth - readw, fragment)
			
						readw -= fragment
			
						ins_depth += fragment
			
						if ins_depth == 32 then
							ins_ptr += 1
							ins_depth = 0
						end
					until readw == 0
		
					b_idx += 1
				end
				
				if b[b_idx] == 0 then
					continue
				end
				
				--print(t_idx, a_idx, b_idx)
				
				acc+=a[a_idx]*b[b_idx]
			elseif ins == 3 then
				--print("-")
				
				
				do
					local readw = idxwidth
		
					repeat
						local fragment = math.min(32 - ins_depth, readw)
			
						a_idx = repl(a_idx, bit32.extract(program[ins_ptr], ins_depth, fragment), idxwidth - readw, fragment)
			
						readw -= fragment
			
						ins_depth += fragment
			
						if ins_depth == 32 then
							ins_ptr += 1
							ins_depth = 0
						end
					until readw == 0
		
					a_idx += 1
				end
				
				if a[a_idx] == 0 then
					continue
				end				
				
				do
					local readw = idxwidth
		
					repeat
						local fragment = math.min(32 - ins_depth, readw)
			
						b_idx = repl(b_idx, bit32.extract(program[ins_ptr], ins_depth, fragment), idxwidth - readw, fragment)
			
						readw -= fragment
			
						ins_depth += fragment
			
						if ins_depth == 32 then
							ins_ptr += 1
							ins_depth = 0
						end
					until readw == 0
		
					b_idx += 1
				end				
				
				if b[b_idx] == 0 then
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

--[=[
local kquat = loadkernel(2, { 0x57fad701, 0xe29eaf48, 0x66cf7789, 0xca })

local t = {}
kquat(t, {0, 0, 0, 1}, {0, 0, 0, 1})
print(string.format("%d, %d, %d, %d", t[1], t[2], t[3], t[4]))
-- -1, 0, 0, 0

]=]



-- armchair benchmark w.r.t loadstring

local __ARMCHAIR = true

if __ARMCHAIR then
	
	local kpga3d = loadkernel(4, { 
	0xce228001,
	0x99e23448,
	0x45eeeab,
	0x36896108,
	0xd8f52d1e,
	0x79ce2f74,
	0xaf27b8eb,
	0x89febbfd,
	0x98e22080,
	0xbb94e0f4,
	0x8c0cdeae,
	0xa4a8c22,
	0xe9a939ea,
	0xeb29c111,
	0x92a3b408,
	0x155e8e8e,
	0x3b88612c,
	0x9a250932,
	0xf97e1b79,
	0xeb12b3ab,
	0x99fafb7d,
	0x18af13c1,
	0x8258d363,
	0x9fa167a9,
	0xb3bb2e9f,
	0xdf9bb2d4,
	0xcb314c1d,
	0xb59d063d,
	0x5a3f709a,
	0xbbea9a9,
	0xf8fafd3f,
	0xca23c221,
	0x9aa024e8,
	0x265e4aa7,
	0x428fb24c,
	0xea290a2b,
	0x2e82a9e3,
	0xe26438d3,
	0xde2a8298,
	0x89b1882e,
	0xb5393e35,
	0xda067ed8,
	0x2ab0ab39,
	0xf4f9ed9f,
	0x9f198331,
	0x549163fc,
	0xa3772dba,
	0x2baeae91,
	0x3b9bd8b0,
	0xe1a8375f,
	0xed1a37cb,
	0x3273d925,
	0xb9a869be,
	0xf96d0b23,
	0x2a83b9f2,
	0xe12488e7,
	0xde0a8a93,
	0xcb61ec3f,
	0x75a92e3c,
	0x6e1e789a,
	0xeb4e969,
	0xf0b86d2f,
	       0})	
	
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
	local t = table.create(16, 0)
	
	for i = 2, 5 do
		a[i] = math.random()
		b[i] = math.random()
	end
	
	--[=[
	local stamp = os.clock()
	
	for i = 1, 1e4 do
		kpga3d_l(t, a, b)
	end
	
	print(string.format("-- kpga3d_l: %.2fms", 1e3* (os.clock() - stamp)))
	]=]
	
	local stamp = os.clock()
	
	for i = 1, 1e4 do
		kpga3d(t, a, b)
	end
	
	for i = 1, 16 do
		t[i] = string.format("%.2f", t[i])
	end
	
	print("--", table.concat(t, ", "))
	
	print(string.format("-- kpga3d: %.2fms", 1e3* (os.clock() - stamp)))
end

-- the armchair benchmark assigns random values to the grade 1 components of a mv pair & computes product

-- pop_stack() not inlined:
-- kpga3d_l: 21.71ms
-- kpga3d: 419.52ms

-- pop_stack() inlined:
-- kpga3d_l: 21.77ms
-- kpga3d: 322.11ms

-- seems to run faster on Roblox, whatever that means
-- on a Lua-based Lua interpreter,~500ms
-- i'm not happy that we're closer to the expectedly slow interpreter than loadstring
-- but it's not easy to compete

-- i will make a new loadkernel that stores a secondary list of indices instead of encoding it with the instructions

local __PRINT_PGA3D_BC = false

if __PRINT_PGA3D_BC then
	local pga3d_bc = compile(4, [[
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
	
	for i = 1, #pga3d_bc do
		pga3d_bc[i] = string.format("0x%x", pga3d_bc[i])
	end

	print(table.concat(pga3d_bc, ",\n"))
end

--[=[
0xce228001,
0x99e23448,
0x45eeeab,
0x36896108,
0xd8f52d1e,
0x79ce2f74,
0xaf27b8eb,
0x89febbfd,
0x98e22080,
0xbb94e0f4,
0x8c0cdeae,
0xa4a8c22,
0xe9a939ea,
0xeb29c111,
0x92a3b408,
0x155e8e8e,
0x3b88612c,
0x9a250932,
0xf97e1b79,
0xeb12b3ab,
0x99fafb7d,
0x18af13c1,
0x8258d363,
0x9fa167a9,
0xb3bb2e9f,
0xdf9bb2d4,
0xcb314c1d,
0xb59d063d,
0x5a3f709a,
0xbbea9a9,
0xf8fafd3f,
0xca23c221,
0x9aa024e8,
0x265e4aa7,
0x428fb24c,
0xea290a2b,
0x2e82a9e3,
0xe26438d3,
0xde2a8298,
0x89b1882e,
0xb5393e35,
0xda067ed8,
0x2ab0ab39,
0xf4f9ed9f,
0x9f198331,
0x549163fc,
0xa3772dba,
0x2baeae91,
0x3b9bd8b0,
0xe1a8375f,
0xed1a37cb,
0x3273d925,
0xb9a869be,
0xf96d0b23,
0x2a83b9f2,
0xe12488e7,
0xde0a8a93,
0xcb61ec3f,
0x75a92e3c,
0x6e1e789a,
0xeb4e969,
0xf0b86d2f,
0x0
]=]