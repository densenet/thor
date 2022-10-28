local util = script.Parent.util

local MAX_D = 8

local BINC, inc_word = require(util.binc)(MAX_D), require(util.inc_word)

local b, x, l = bit32.btest, bit32.bxor, bit32.lshift

--local gbd = {1, 2, 6, 12, 16, 17}
--local gl = {1, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 16}

-- idx: the 1-indexed position in an array-based multivector in which a multivector component exists
type idx = number
-- grade_idx: the 1-indexed position of the nth grade word in a list of all grade words. node n starts from zero. this is generally n+1
type grade_idx = number
-- word: the bit representation of a multivector component basis.
--       bit words have a intuitive & concise 1-1 correspondence: on bit refers to presence of a wedged basis vector
type word = number
-- grade_word: similar to word, except represents all grades in the multivector which do not have all components zero.
--             a multivector in 4D which only consists of a sum of vectors has a grade word 0b01000.
type grade_word = number

type GradeMaps = {
  GRADE_IDX_TO_FIRST_IDX: { [grade_idx]: idx },
  IDX_TO_GRADE_WORD: { [idx]: grade_word },
  WORD_TO_IDX: { [word]: idx },
  IDX_TO_WORD: { [idx]: word }
}

local get_sigmaps
do
  local function generate_sigmap_tables(D)
    local gi2i = table.create(D + 2)
    local fr = table.freeze
    
    -- generate the map that takes grade index to first grade-associated MV component index
    -- if the grade is N, it will return the component index of e_0^e_1^...*e_(N-1)
    do
      local t = gi2i
      t[1] = 1
      
      for k = 1, D + 1 do
        t[k + 1] = BINC[D][k] + t[k]
      end
      
      fr(t)
    end
    
    local i2gw = {}
    do
      local t = i2gw
      
      for k = 1, D do
        for i = BINC[D][k], BINC[D][k + 1] do
          t[i] = l(2, k - 1)
        end
      end
      
      fr(t)
    end
    
    local w2i, i2w = {}, {}
    do
      local t, u = w2i, i2w
      
      for k = 1, D do
        local w = l(1, k - 1) - 1
        t[w], u[k] = k, w
        for i = 1, BINC[D][k] - 1 do
          w = inc_word(w, D)
          t[w], u[k + i] = k + i, w
        end
      end
      
      fr(t); fr(u)
    end
    
    return gi2i, i2gw, w2i, i2w
  end 
  
  local function f(signature: Signature): GradeMaps
    local gi2i, i2gw, w2i, i2w = generate_sigmap_tables(signature.Null + signature.Unit + signature.Anti)
    
    return {
      GRADE_IDX_TO_FIRST_IDX = gi2i,
      IDX_TO_GRADE_WORD = i2gw,
      WORD_TO_IDX = w2i,
      IDX_TO_WORD = i2w
    }
  end
  
  get_sigmaps = f
end

local function get_prod_term_sign(signature: Signature, bw_i: number, bw_j: number): number
  local retsign = 1

  local D0 = signature.Null
  local D1 = signature.Unit + D0
  local D2 = signature.Anti + D1

  -- null signature
  local common_basis = bit32.band(bw_i, bw_j)
  for pos = 0, D0 - 1 do
    if b(common_basis, l(1, pos)) then
      return 0
    end
  end
  -- sign produced by reducing (bw_i,bw_j) to product basis
  -- aka epsilon
  do
    local acc = false
    local negate = false
    for i = D2 - 1, 1, -1 do
      acc = if b(bw_i, l(1, i)) then not acc else acc
      if b(bw_j, l(1, i - 1)) and acc then
        negate = not negate
      end
    end
    retsign = if negate then -retsign else retsign
  end
  -- anti signature
  for pos = D1, D2 - 1 do
    if b(common_basis, l(1, pos)) then
      retsign = -retsign
    end
  end

  return retsign
end

export type Code<T> = T

-- generate a geometric product kernel
local GATEOFF = 0x20
local f = string.format
local function kgen_geo_prod(signature: Signature, gate: number): Code<string>
  local buf1 = {}
  local buf2 = {
    T = {},
  }

  local SIGMAP = get_sigmaps(signature)

  local gbd, gl = SIGMAP.GRADE_IDX_TO_FIRST_IDX, SIGMAP.IDX_TO_GRADE_WORD
  local bl, bi = SIGMAP.IDX_TO_WORD, SIGMAP.WORD_TO_IDX

  for _i = 1, #gbd - 1 do
    -- ith product component
    for i = gbd[_i], gbd[_i + 1] - 1 do
      -- tg_i: ith basis word
      local tg_i = bl[i]
      local term_ct = 0
      for _j = 1, #gbd - 1 do
        -- enforce rhs gate
        if not b(gate, l(1, _j - 1) * GATEOFF) then
          continue
        end

        -- jth term in ith product component
        for j = gbd[_j], gbd[_j + 1] - 1 do
          -- bw_j: jth basis word (rhs component)
          local bw_j = bl[j]

          local bw_i = x(tg_i, bw_j)
          local i_swp = bi[bw_i]

          -- enforce lhs gate
          if not b(gate, gl[i_swp]) then
            continue
          end

          -- enforce signature
          local sign = get_prod_term_sign(signature, bw_i, bw_j)
          if sign == 0 then
            continue
          end

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
--       on-bit refers to zero components in the respective grade.
--       bit position increases with grade dimension, up to the 5th bit.
--       0b11101_11101 : kernel to multiply 4D multivectors of grade 1 (pure 4D vectors)

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

local STANDARD_SIG = {
  pga3d = { Unit = 3, Anti = 0, Null = 1 }
}

return function(_signature: (Signature | string), _gate: number): Kernel<string>
  local gate = if _gate then _gate else 0
  
  local signature: Signature = if typeof(_signature) == "string" then STANDARD_SIG[_signature] else _signature

  local DU, DA, DN = signature.Unit, signature.Anti, signature.Null
  if DU + DA + DN > MAX_D then
    error(f("a signature with total dimension greater than 8 is not supported; got (%d, %d, %d)", DU, DA, DN))
  end
  if DU % 1 ~= 0 or DA % 1 ~= 0 or DN % 1 ~= 0 then
    error(f("a signature can't be composed of fractional numbers; got (%.2f, %.2f, %.2f)", DU, DA, DN))
  end
  if DU < 0 or DA < 0 or DN < 0 then
    error(f("a signature must be composed of non negative numbers; got (%d, %d, %d)", DU, DA, DN))
  end
  local kernel = {}

  kernel.__mul = kgen_geo_prod(signature, gate)

  kernel.__pow = _stub()

  kernel.__concat = _stub()

  return kernel
end
