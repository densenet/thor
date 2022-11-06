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
-- positive integer representing signature dimensions
type dim = number
-- sign: -1 | 0 | 1
type sign = number

export type GradeMaps = {
  GRADE_IDX_TO_FIRST_IDX: { [grade_idx]: idx },
  IDX_TO_GRADE_WORD: { [idx]: grade_word },
  WORD_TO_IDX: { [word]: idx },
  IDX_TO_WORD: { [idx]: word }
}

local MAX_D = 16

local BINC, inc_word
if script then
  local util = script.Parent.util  
  
  BINC, inc_word = require(util.binc), require(util.inc_word)
else
  BINC, inc_word = require("src\\util\\binc"), require("src\\util\\inc_word")
end
BINC = BINC(MAX_D)

--local bit2str = require("src\\util\\bit2str")

local b, x, l, n, r = bit32.btest, bit32.bxor, bit32.lshift, bit32.bnot, bit32.rshift
local rz = bit32.countrz

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
  
  local function f(DU: dim, DA: dim, DN: dim): GradeMaps
    local gi2i, i2gw, w2i, i2w = generate_sigmap_tables(DU+DA+DN)
    
    return {
      GRADE_IDX_TO_FIRST_IDX = gi2i,
      IDX_TO_GRADE_WORD = i2gw,
      WORD_TO_IDX = w2i,
      IDX_TO_WORD = i2w
    }
  end
  
  get_sigmaps = f
end

  -- sign produced by reducing (bw_i,bw_j) to product basis
  -- follows fundamental rules of multiplying two MV component bases: eliminate the shared bases, flip on the non-shared bases
local function eps(D: dim, bw_i: word, bw_j: word): sign
  local acc = false
  local negate = false
  for i = D - 1, 1, -1 do
    acc = if b(bw_i, l(1, i)) then not acc else acc
    if b(bw_j, l(1, i - 1)) and acc then
      negate = not negate
    end
  end
  return if negate then -1 else 1
end

local function get_prod_term_sign(DU: dim, DA: dim, DN: dim, bw_i: word, bw_j: word): sign
  local retsign = 1

  local D0, D1, D2 = DN, DU + DN, DU + DN + DA
  
  -- cobw gives us information about which basis vectors are shared in the MV component basis word
  local cobw = bit32.band(bw_i, bw_j)
  
  -- null signature short-circuit:
  -- if any two null signature elements are multiplied, we get wedge product case: it eliminates the product term
  for pos = 0, D0 - 1 do
    if b(cobw, l(1, pos)) then
      return 0
    end
  end

  retsign *= eps(D2, bw_i, bw_j)
    
  -- anti signature bookkeeping: just have to flip signs if two basis vectors annihilate
  for pos = D1, D2 - 1 do
    if b(cobw, l(1, pos)) then
      retsign = -retsign
    end
  end

  return retsign :: sign
end

local f = string.format
return function(opt, args)
  local DU, DA, DN = table.unpack(string.split(opt, ","))
  
  DU, DA, DN = tonumber(DU) or 0, tonumber(DA) or 0, tonumber(DN) or 0
  
  if DU % 1 ~= 0 or DA % 1 ~= 0 or DN % 1 ~= 0 then
    error(f("a sig can't be composed of fractional numbers; got (%.2f, %.2f, %.2f)", DU, DA, DN))
  end
  if DU < 0 or DA < 0 or DN < 0 then
    error(f("a sig must be composed of non negative numbers; got (%d, %d, %d)", DU, DA, DN))
  end
  
  local D = DU+DA+DN
  
  local gate = 0
  
  if args.gate then
    local gb = string.split(args.gate, ",")
    
    for i, v in gb do
      v = tonumber(v)
      gate = x(gate, l(1, v - 1))
    end
  end
  
  local sig = {
    D = D,
  }
  
  local SIGMAP = get_sigmaps(DU, DA, DN)
  
  local WORD2I, IDX2GW = SIGMAP.WORD_TO_IDX, SIGMAP.IDX_TO_GRADE_WORD
  
  local function prod_initialize(bw_i, bw_j)
    local j: idx = WORD2I[bw_j]
  
    bw_i = x(bw_i, bw_j) :: word
    local i_swp: idx = WORD2I[bw_i]
    local i_gw: grade_word = IDX2GW[i_swp]

    if b(gate, i_gw) then
      return -- case: gate says left multiplicand's current grade should be all zero
    end
    
    return bw_i, bw_j, i_swp, j
  end
  
  local function geometric_product(state, _bw_i: word, _bw_j: word)
    local bw_i, bw_j, i, j = prod_initialize(_bw_i, _bw_j)
    
    if not bw_i then
      return
    end
    
    local s: sign = get_prod_term_sign(DU, DA, DN, bw_i, bw_j)
    if s == 0 then
      return
    end
    
    state.term_ct += 1
       
    state.term_buf[state.term_ct] = f("%sa[%d]*b[%d]", if s == 1 then if state.term_ct > 1 then "+" else "" else "-", i, j)   
  end
  
  local function wedge(state, _bw_i: word, _bw_j: word)
    local bw_i, bw_j, i, j = prod_initialize(_bw_i, _bw_j)
    
    if not bw_i then
      return
    end
      
    local s: sign = get_prod_term_sign(0, 0, D, bw_i, bw_j)
    if s == 0 then
      return
    end
    
    state.term_ct += 1
    
    state.term_buf[state.term_ct] = f("%sa[%d]*b[%d]", if s == 1 then if state.term_ct > 1 then "+" else "" else "-", i, j)
  end
  
  local BW_IM = n(l(-1, D))
  local function regressive(state, tg_i: word, _bw_j: word)
    local tg_ir, bw_jr = x(BW_IM, tg_i), x(BW_IM, _bw_j)
     
    local s_tg, s_j = eps(D, tg_ir, tg_i), eps(D, _bw_j, bw_jr)

    local bw_ir = prod_initialize(tg_ir, bw_jr)
    
    local bw_i = x(BW_IM, bw_ir)
    
    local s_i = eps(D, bw_i, bw_ir)
    
    if not bw_i then
      return
    end
    
    local s: sign = get_prod_term_sign(0, 0, D, bw_ir, bw_jr)
    if s == 0 then
      return
    end
    
    s *= s_tg * s_j * s_i
    
    local i, j = WORD2I[bw_i], WORD2I[_bw_j]
    
    state.term_ct += 1
    
    state.term_buf[state.term_ct] = f("%sa[%d]*b[%d]", if s == 1 then if state.term_ct > 1 then "+" else "" else "-", i, j)
  end
  
  local function prod_finalize(state, tg_i, k_buf)
    local i: idx = WORD2I[tg_i]
      
    -- product term, composed of all relevant term entries
    k_buf[#k_buf + 1] = f("t[%d]=%s", i, table.concat(state.term_buf))
  end
  
  local function _stub() end
  
  local HANDLER_MAP = {
    product = geometric_product,
    wedge = wedge,
    regressive = regressive
  }
  
  function sig:handler(op, state)
    state.grade = 1
    if HANDLER_MAP[op] then
      return function(...)
        HANDLER_MAP[op](state, ...)
      end, function(...)
        prod_finalize(state, ...)
      end
    end
    return _stub
  end
  
  function sig.inc_word(state, bw)
    bw = inc_word(bw, D)
    
    local i = WORD2I[bw]
    local gw = IDX2GW[i]
    
    if rz(n(bw)) == rz(gw) then
      while true do
        gw = l(gw, 1)
        if not b(gate, gw) then
          return bw
        end
        bw = n(l(n(bw), 1))
      end
    else
      return bw
    end
    return 0
  end
  
  return sig
end
