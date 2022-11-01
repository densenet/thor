local MAX_D = 16

local BINC, inc_word
if script then
  local util = script.Parent.util  
  
  BINC, inc_word = require(util.binc), require(util.inc_word)
else
  BINC, inc_word = require("src\\util\\binc"), require("src\\util\\inc_word")
end
BINC = BINC(MAX_D)

local b, x, l, n = bit32.btest, bit32.bxor, bit32.lshift, bit32.bnot

local rz = bit32.countrz

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

type dim = number

type GradeMaps = {
  GRADE_IDX_TO_FIRST_IDX: { [grade_idx]: idx },
  IDX_TO_GRADE_WORD: { [idx]: grade_word },
  WORD_TO_IDX: { [word]: idx },
  IDX_TO_WORD: { [idx]: word }
}

local get_sigmaps
do
  local function generate_sigmap_tables(D: dim): any
    local gi2i = table.create(D + 2)
    local fr = table.freeze
    
    -- generate the map that takes grade index to first grade-associated MV component index
    -- if the grade is N, it will return the component index of e_0^e_1^...*e_(N-1)
    do
      local t = gi2i
      t[1] = 1
      
      for k = 1, D + 1 do
        t[k + 1] = BINC[D + 1][k] + t[k]
      end
      
      fr(t)
    end
    
    -- generate the map that takes MV component index to associated grade word
    local i2gw = {}
    do
      local t = i2gw
      
      for k = 1, D + 1 do
        for i = gi2i[k], gi2i[k + 1] - 1 do
          t[i] = l(1, k - 1)
        end
      end
      
      fr(t)
    end
    
    -- generate the map that takes MV component basis word to respective MV component index, and the reverse map
    local w2i, i2w = {}, {}
    do
      local t, u = w2i, i2w
      
      local i = 1
      local w = 0
      repeat
        t[w], u[i] = i, w
        w = inc_word(w, D)
        i += 1
      until w == 0
      
      fr(t); fr(u)
    end
    
    return gi2i, i2gw, w2i, i2w
  end 
  
  local function f(sig: Signature): GradeMaps
    local gi2i, i2gw, w2i, i2w = generate_sigmap_tables(sig.Null + sig.Unit + sig.Anti)
    
    return {
      GRADE_IDX_TO_FIRST_IDX = gi2i,
      IDX_TO_GRADE_WORD = i2gw,
      WORD_TO_IDX = w2i,
      IDX_TO_WORD = i2w
    }
  end
  
  get_sigmaps = f
end

local function get_prod_term_sign(sig: Signature, bw_i: number, bw_j: number): number
  local retsign = 1

  local D0 = sig.Null
  local D1 = sig.Unit + D0
  local D2 = sig.Anti + D1
  
  -- cobw gives us information about which basis vectors are shared in the MV component basis word
  local cobw = bit32.band(bw_i, bw_j)
  
  -- null signature short-circuit:
  -- if any two null signature elements are multiplied, we get wedge product case: it eliminates the product term
  for pos = 0, D0 - 1 do
    if b(cobw, l(1, pos)) then
      return 0
    end
  end
  -- sign produced by reducing (bw_i,bw_j) to product basis
  -- follows fundamental rules of multiplying two MV component bases: eliminate the shared bases, flip on the non-shared bases
  -- inlined version of ./util/eps.luau
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
  -- anti signature bookkeeping: just have to flip signs if two basis vectors annihilate
  for pos = D1, D2 - 1 do
    if b(cobw, l(1, pos)) then
      retsign = -retsign
    end
  end

  return retsign
end

export type Code<T> = T

local function expand(D: dim, gate: grade_word, cb: (word, word) -> (), cb2: (word) -> ())
  local bw_A, bw_B = 0, 0
  
  local grade_A: grade_idx, grade_B: grade_idx = 1, 1
  repeat
    bw_B = 0
    repeat
      cb(bw_A, bw_B)
      
      bw_B = inc_word(bw_B, D)
      
      if rz(n(bw_B)) == grade_B then
        while true do
          grade_B += 1
          if not b(gate, l(1, D + grade_B)) then
            break
          end
          bw_B = n(l(n(bw_B), 1))
        end
      end
    until bw_B == 0 or grade_B > D + 1
    cb2(bw_A)
   
    bw_A = inc_word(bw_A, D)
  until bw_A == 0 or grade_A > D + 1
end

-- generate a geometric product kernel
local f = string.format
local function kgen_geo_prod(KGEN_ENV): Code<string>
  local D: dim = KGEN_ENV.D
  local SIGMAP: GradeMaps = KGEN_ENV.SIGMAP
  local sig: Signature = KGEN_ENV.Signature
  local gate: grade_word = KGEN_ENV.Gate
  
  local _GRADE2I, IDX2GW = SIGMAP.GRADE_IDX_TO_FIRST_IDX, SIGMAP.IDX_TO_GRADE_WORD
  local _IDX2W, WORD2I = SIGMAP.IDX_TO_WORD, SIGMAP.WORD_TO_IDX
  
  local buf1: { Code<string> } = table.create(#IDX2GW)
  do
    local buf2: { Code<string> } = table.create(#IDX2GW)
    
    local term_ct = 0
   
    expand(D, gate, function(bw_i: word, bw_j: word)
      local j: idx = WORD2I[bw_j]
    
      -- initially, bw_i is actually the product's output basis word
      -- reverse calculate which bw_i would need to be lhs for a given rhs bw_j to get bw_i as product basis word
      bw_i = x(bw_i, bw_j) :: word
      -- "_swp" encodes the trick we did above. it's the multivector component index of calculated bw_i
      local i_swp: idx = WORD2I[bw_i]
      

      local i_gw: grade_word = IDX2GW[i_swp]

      if b(gate, i_gw) then
        return -- case: gate says left multiplicand's current grade should be all zero
      end
     
      local sign: number = get_prod_term_sign(sig, bw_i, bw_j)
      if sign == 0 then
        return -- case: some common multivector component had a null basis vector
      end
      
      term_ct += 1
      
      -- term entry
      buf2[term_ct] = f("%sa[%d]*b[%d]", if sign == 1 then if term_ct > 1 then "+" else "" else "-", i_swp, j)
    end, function(tg_i: word)
      if term_ct > 0 then
        local i: idx = WORD2I[tg_i]
        
        -- product term, composed of all relevant term entries
        buf1[#buf1 + 1] = f("t[%d]=%s", i, table.concat(buf2))
        
        table.clear(buf2)
        
        term_ct = 0
      end
    end)
  end

  return f("function(t,a,b)\n%s\nend", table.concat(buf1, "\n"))
end

local function kgen_wedge_prod(KGEN_ENV): Code<string>
  local D: dim = KGEN_ENV.D
  local SIGMAP: GradeMaps = KGEN_ENV.SIGMAP
  local sig: Signature = KGEN_ENV.Signature
  local gate: grade_word = KGEN_ENV.Gate
  
  local buf1: { Code<string> } = {}
  do
    local buf2: { Code<string> } = {}
    
    
  end
  
  
  return f("function(t,a,b)\n%s\nend", table.concat(buf1, "\n"))
end

local function _stub(): Code<string>
  return ""
end

export type Kernel<T> = {
  __mul: Code<T>,
  __pow: Code<T>,
  __concat: Code<T>,
  __mod: Code<T>,
  __unm: Code<T>
}

export type Signature = {
  Unit: number,
  Anti: number,
  Null: number,
}

local STANDARD_SIG = {
  c = { Unit = 0, Anti = 1, Null = 0 },
  hyp = { Unit = 1, Anti = 0, Null = 0 },
  diff = { Unit = 0, Anti = 0, Null = 1 },
  quat = { Unit = 0, Anti = 2, Null = 0 },
  pga3d = { Unit = 3, Anti = 0, Null = 1 },
  cga3d = { Unit = 4, Anti = 1, Null = 0 },
  octo = { Unit = 0, Anti = 3, Null = 0 },
  dquat = { Unit = 0, Anti = 2, Null = 1 },
}

-- signature: all data relevant to exactly what the behavior of the algebra should be.
--            decisions about the outcome of multiplying basis elements determine the algebra.
--            it's important to study 

-- gate: 10 bit number representing two juxtaposed grade bases.
--       the 5 lsbs is a left multiplier, and 5 msbs right multiplier.
--       on-bit refers to zero components in the respective grade.
--       bit position increases with grade dimension, up to the 5th bit.

local function genkernel(_signature: (Signature | string), _gate: grade_word?): Kernel<string>
  local gate = if _gate then _gate else 0
  
  local sig: Signature = if typeof(_signature) == "string" then STANDARD_SIG[_signature] else _signature

  local DU, DA, DN = sig.Unit, sig.Anti, sig.Null
  local D = DU+DA+DN
  if D > MAX_D then
    error(f("a sig with total dimension greater than 8 is not supported; got (%d, %d, %d)", DU, DA, DN))
  end
  if DU % 1 ~= 0 or DA % 1 ~= 0 or DN % 1 ~= 0 then
    error(f("a sig can't be composed of fractional numbers; got (%.2f, %.2f, %.2f)", DU, DA, DN))
  end
  if DU < 0 or DA < 0 or DN < 0 then
    error(f("a sig must be composed of non negative numbers; got (%d, %d, %d)", DU, DA, DN))
  end
  
  local kernel = {}
  
  do
    local SIGMAP: GradeMaps = get_sigmaps(sig)
    
    local KGEN_ENV = {
      D = D,
      SIGMAP = SIGMAP,
      Signature = sig,
      Gate = gate
    }
  
    kernel.__mul = kgen_geo_prod(KGEN_ENV) -- GEOMETRIC PRODUCT
    kernel.__pow = kgen_wedge_prod(KGEN_ENV) -- WEDGE PRODUCT / MEET
    kernel.__concat = _stub() -- REGRESSIVE PRODUCT / JOIN
    kernel.__len = _stub() -- MV DUAL
    kernel.__unm = _stub() -- UNARY MINUS
    kernel.__mod = _stub() -- SANDWICH PRODUCT (A % B == B * A * B.r)
  end

  return kernel :: Kernel<string>
end

--[=[
do
  local k = genkernel("quat", 1)
  
  print(k.__mul)
end
]=]

return genkernel
