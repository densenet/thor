local b, l = bit32.btest, bit32.lshift

local LOOKUP8 = (function()
	local ret = {}
	
	for i = 0, 255 do
		local chars = table.create(8)
		
		for j = 0, 7 do
			table.insert(chars, 1, if b(l(1, j), i) then "1" else "0")
		end
		
		ret[i] = chars
	end

	return ret
end)()

local function bit2str(n: number): string
	local s = table.create(32)
	
	do
		local q = 0
		local ptr = 0
		repeat
			table.move(s, 1, ptr, 9)

			q = bit32.extract(n, 0, 8)
			n = bit32.arshift(n, 8)
			s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8]
				= table.unpack(LOOKUP8[q])
			ptr += 8
		until n == 0
	end
	
	return table.concat(s)
end

return bit2str
