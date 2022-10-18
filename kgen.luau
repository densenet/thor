local bi = {
	[0b0000] = 1, -- one
  [0b0001] = 2, -- w
	[0b0010] = 3, -- x
	[0b0100] = 4, -- y
	[0b1000] = 5, -- z
  [0b0011] = 6, -- wx
  [0b0101] = 7, -- wy
  [0b1001] = 8, -- wz
	[0b0110] = 9, -- xy
	[0b1010] = 10, -- xz
	[0b1100] = 11, -- yz
  [0b0111] = 12, -- wxy
  [0b1011] = 13, -- wxz
  [0b1101] = 14, -- wyz
  [0b1110] = 15, -- xyz
	[0b1111] = 16 -- wxyz
}

local bl = {
  0b0000, -- scalar component
  0b0001, 0b0010, 0b0100, 0b1000, --  vector components
  0b0011, 0b0101, 0b1001, 0b0110, 0b1010, 0b1100, -- bi-vector components
  0b0111, 0b1011, 0b1101, 0b1110, -- tri-vector / pseudovector components
  0b1111 -- pseudoscalar component
}

--[=[
local bt = {
	"_", "w", "x", "y", "z", "wx", "wy", "wz", "xy", "xz", "yz", "wxy", "wxz", "wyz", "xyz", "wxyz"
}
]=]

local b = bit32.btest

local e = (function()
  local ret = {}
  for i = 0, 255 do
    local value = 0
    if b(i, 0x2) then
      if b(i, 0x10) then
        value += 1
      end
    end
    if b(i, 0x4) then
      if b(i, 0x10) then
        value += 1
      end
      if b(i, 0x20) then
        value += 1
      end
    end
    if b(i, 0x8) then
      if b(i, 0x10) then
        value += 1
      end
      if b(i, 0x20) then
        value += 1
      end
      if b(i, 0x40) then
        value += 1
      end
    end
    if bit32.band(value, 1) == 1 then
      ret[i] = -1
    else
      ret[i] = 1
    end
  end

  -- to get index of E with two basis nibbles:
  -- bxor(lshift(right_basis, 4), left_basis)

  return ret
end)()

local bitct = (function()
  local ret = {}

  for i = 0, 255 do
    ret[i] = i%2 + (i%4 - i%2)/2  + (i%8 - i%4)/4 + (i%16 - i%8)/8 + (i%32 - i%16)/16 + (i%64 - i%32)/32 + (i%128 - i%64)/64 + (i%256 - i%128)/128
  end

  return ret
end)()

local gbit = {1, 2, 4, 8, 16}
local gbd = {1, 2, 6, 12, 16, 17}
local gl = {1, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 16}

local x = bit32.bxor

-- template for multiplication
--[=[
local i_1
for i = 1, 16 do
  i_1 = bl[i]
  t_prod[i] += t1[i] * t2[1]
    + e[x(i_1, 0x10)] * (t1[bi[x(i_1, 0x1)]] * t2[2])  + e[x(i_1, 0x20)] * (t1[bi[x(i_1, 0x2)]] * t2[3])  + e[x(i_1, 0x40)] * (t1[bi[x(i_1, 0x4)]] * t2[4])  + e[x(i_1, 0x80)] * (t1[bi[x(i_1, 0x8)]] * t2[5])
    + e[x(i_1, 0x30)] * (t1[bi[x(i_1, 0x3)]] * t2[6])  + e[x(i_1, 0x50)] * (t1[bi[x(i_1, 0x5)]] * t2[7])  + e[x(i_1, 0x90)] * (t1[bi[x(i_1, 0x9)]] * t2[8])  + e[x(i_1, 0x60)] * (t1[bi[x(i_1, 0x6)]] * t2[9])  + e[x(i_1, 0xA0)] * (t1[bi[x(i_1, 0xA)]] * t2[10]) + e[x(i_1, 0xC0)] * (t1[bi[x(i_1, 0xC)]] * t2[11])
    + e[x(i_1, 0x70)] * (t1[bi[x(i_1, 0x7)]] * t2[12]) + e[x(i_1, 0xB0)] * (t1[bi[x(i_1, 0xB)]] * t2[13]) + e[x(i_1, 0xD0)] * (t1[bi[x(i_1, 0xD)]] * t2[14]) + e[x(i_1, 0xE0)] * (t1[bi[x(i_1, 0xE)]] * t2[15])
    + e[x(i_1, 0xF0)] * (t1[bi[x(i_1, 0xF)]] * t2[16])
end
]=]

local function get_signature(signature: number, wbasis0: number, wbasis1: number): number
  local retsign = 1

  local common_basis = bit32.band(wbasis0, wbasis1)
  for pos = 1, 4 do
    if b(common_basis, gbit[pos]) then
      local sigdata = bit32.extract(signature, (pos - 1) * 2, 2)
      if sigdata == 1 then
        retsign = -retsign
      elseif sigdata == 2 then
        retsign = 0
        break
      end
    end
  end

  return retsign
end

export type Code<T> = T

-- generate a geometric product kernel
local GATEOFF = 0x20
local f = string.format
local function kgen_geo_prod(gate: number, signature: number): Code<string>
  local buf1 = {}
  local buf2 = {
    T = {},
  }

  for _i = 1, #gbd - 1 do
    --[=[
    if not b(gate, gbit[_i]) then
      continue
    end]=]

    -- ith product component
    for i = gbd[_i], gbd[_i + 1] - 1 do
      -- tg_i: ith basis word
      local tg_i = bl[i]
      local term_ct = 0
      for _j = 1, #gbd - 1 do
        -- enforce rhs gate
        if not b(gate, gbit[_j] * GATEOFF) then
          continue
        end

        -- jth term in ith product component
        for j = gbd[_j], gbd[_j + 1] - 1 do
          -- bw_j: jth basis word (rhs component)
          local bw_j = bl[j]

          local bw_i = x(tg_i, bw_j)
          local sign = e[x(tg_i, bw_j * 0x10)]
          local i_swp = bi[bw_i]

          -- enforce lhs gate
          if not b(gate, gl[i_swp]) then
            continue
          end

          -- enforce signature
          sign = sign * get_signature(signature, bw_i, bw_j)
          if sign == 0 then
            continue
          end
          --[=[
          local bw_co = bit32.band(bw_i, bw_j)
          local nulled = false
          for pos = 1, 4 do
            if b(bw_co, gbit[pos]) then
              local sigdata = bit32.extract(signature, (pos - 1) * 2, 2)
              if sigdata == 1 then
                sign = -sign
              elseif sigdata == 2 then
                nulled = true
                break
              end
            end
          end
          if nulled then
            continue
          end]=]
          --[=[
          if i_swp == 2 and j == 2 then
            print(string.format("%x,%x", bw_i, bw_j))
          end]=]

          term_ct += 1

          buf2.T[term_ct] = f("%sa[%d]*b[%d]", if sign == 1 then if term_ct > 1 then "+" else "" else "-", i_swp, j)
        end
      end

      if term_ct > 0 then
        buf1[#buf1 + 1] = f("t[%d]=%s", i, table.concat(buf2.T))
      end
      table.clear(buf2.T)
    end
  end

  return f("function(t,a,b)\n%s\nend", table.concat(buf1, "\n"))
end

local function _stub(): Code<string>
  return ""
end

-- gate: 10 bit number representing two juxtaposed grade bases.
--       the 5 lsbs is a left multiplier, and 5 msbs right multiplier.
--       on-bit refers to nonzero components in the respective grade.
--       bit position increases with grade dimension, up to the 5th bit.
--       0b0001000010 : kernel to multiply 4D multivectors of grade 1 (pure 4D vectors)
-- signature: determines how a given 4D basis vector multiplies against itself wrt geometric product.
--            two bits per 4D vector basis = 8 bits.
--            00 => e_n*e_n=1
--            01 => e_n*e_n=-1
--            10 => e_n*e_n=0

export type Kernel<T> = {
  __mul: Code<T>,
  __pow: Code<T>,
  __concat: Code<T>
}

print(kgen_geo_prod(0b11111_11111, 0))

return function(gate: number, signature: number): Kernel<string>
  local kernel = {}

  kernel.__mul = kgen_geo_prod(gate, signature)

  kernel.__pow = _stub()

  kernel.__concat = _stub()

  return kernel :: Kernel<string>
end
