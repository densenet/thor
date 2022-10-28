local util = script.Parent.util

local BASIS = require(util.basis)
local EPS = require(util.eps)

local bi = BASIS.IDX
local bl = BASIS.LIST

--[=[
local bt = {
	"_", "w", "x", "y", "z", "wx", "wy", "wz", "xy", "xz", "yz", "wxy", "wxz", "wyz", "xyz", "wxyz"
}
]=]

local b = bit32.btest

local e = EPS.MUL

local gbit = {1, 2, 4, 8, 16}
local gbd = {1, 2, 6, 12, 16, 17}
local gl = {1, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 16}

local x = bit32.bxor

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
local function kgen_geo_prod(_gate: number?, _signature: number?): Code<string>
  local gate: number = if _gate then _gate else 0x3ff
  local signature: number = if _signature then _signature else 0
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
          local i_swp = bi[bw_i]
          local sign = e[x(bw_i, bw_j * 0x10)]

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
--       0b00010_00010 : kernel to multiply 4D multivectors of grade 1 (pure 4D vectors)
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
--print(kgen_geo_prod(0b00001_00001, 0))

export type Signature = {
  Unit: number,
  Anti: number,
  Null: number,
}

local PGA3D = {
  Unit = 3,
  Anti = 0,
  Null = 1
}

return function(_gate: number?, _signature: Signature?): Kernel<Signature>
  local gate = if _gate then _gate else 0
  local signature = if _signature then _signature else PGA3D
  local kernel = {}

  kernel.__mul = kgen_geo_prod(gate, signature)

  kernel.__pow = _stub()

  kernel.__concat = _stub()

  return kernel
end
