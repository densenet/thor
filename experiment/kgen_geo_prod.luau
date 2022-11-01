local b, x, l = bit32.btest, bit32.bxor, bit32.lshift

local MAX_D = 8

local BINC, inc_word = require("src\\util\\binc")(MAX_D), require("src\\util\\inc_word")

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

local function get_prod_term_sign(signature: Signature, bw_i: number, bw_j: number): number
  local retsign = 1

  local D0 = signature.Null
  local D1 = signature.Unit + D0
  local D2 = signature.Anti + D1
  
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

type idx = number
type grade_idx = number
type word = number
type grade_word = number

type dim = number

type GradeMaps = {
  GRADE_IDX_TO_FIRST_IDX: { [grade_idx]: idx },
  IDX_TO_GRADE_WORD: { [idx]: grade_word },
  WORD_TO_IDX: { [word]: idx },
  IDX_TO_WORD: { [idx]: word }
}

type Signature = {
  Unit: number,
  Anti: number,
  Null: number,
}

local function expand(D: dim, cb: (word, word) -> unknown, cb2: (word) -> unknown)
  local bw_A, bw_B = 0, 0
  
  repeat
    bw_B = 0
    repeat
      cb(bw_A, bw_B)
      
      bw_B = inc_word(bw_B, D)
    until bw_B == 0
    cb2(bw_A)
    
    bw_A = inc_word(bw_A, D)
  until bw_A == 0
end

local f = string.format
local function kgen_geo_prod(KGEN_ENV): string
  local D: dim = KGEN_ENV.D
  local SIGMAP: GradeMaps = KGEN_ENV.SIGMAP
  local sig: Signature = KGEN_ENV.Signature
  local gate: grade_word = KGEN_ENV.Gate
  
  local _GRADE2I, IDX2GW = SIGMAP.GRADE_IDX_TO_FIRST_IDX, SIGMAP.IDX_TO_GRADE_WORD
  local _IDX2W, WORD2I = SIGMAP.IDX_TO_WORD, SIGMAP.WORD_TO_IDX
  
  local buf1 = table.create(#IDX2GW)
  do
    local buf2 = table.create(#IDX2GW)
    
    local term_ct = 0
   
    expand(D, function(bw_i: word, bw_j: word)
      local j: idx = WORD2I[bw_j]
    
      -- initially, bw_i is actually the product's output basis word
      -- reverse calculate which bw_i would need to be lhs for a given rhs bw_j to get bw_i as product basis word
      bw_i = x(bw_i, bw_j) :: word
      -- "_swp" encodes the trick we did above. it's the multivector component index of calculated bw_i
      local i_swp: idx = WORD2I[bw_i]
      
      local i_gw: grade_word, j_gw: grade_word = IDX2GW[i_swp], IDX2GW[j]
      if b(gate, i_gw) or b(gate, j_gw * 0x20) then
        return -- case: gate says one multiplicand's current grade should be all zero
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

local pga3d: Signature = { Unit = 3, Anti = 0, Null = 1 }

local KGEN_ENV_PGA3D = {
	D = 4,
	SIGMAP = get_sigmaps(pga3d),
	Signature = pga3d,
	Gate = 0
}

local cga: Signature = { Unit = 4, Anti = 1, Null = 0 }

local KGEN_ENV_CGA3D = {
	D = 5,
	SIGMAP = get_sigmaps(cga),
	Signature = cga,
	Gate = 0
}

local wedge4: Signature = { Unit = 0, Anti = 0, Null = 4 }

local KGEN_ENV_WEDGE4 = {
	D = 4,
	SIGMAP = get_sigmaps(wedge4),
	Signature = wedge4,
	Gate = 0
}

local quat: Signature = { Unit = 0, Anti = 2, Null = 0 }

local KGEN_ENV_QUAT = {
	D = 2,
	SIGMAP = get_sigmaps(quat),
	Signature = quat,
	Gate = 0
}

print(kgen_geo_prod(KGEN_ENV_CGA3D))

--[=[

sig: pga3d
function(t,a,b)
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

sig: cga3d
function(t,a,b)
t[1]=a[1]*b[1]+a[2]*b[2]+a[3]*b[3]+a[4]*b[4]+a[5]*b[5]-a[6]*b[6]-a[7]*b[7]-a[8]*b[8]-a[9]*b[9]+a[10]*b[10]-a[11]*b[11]-a[12]*b[12]+a[13]*b[13]-a[14]*b[14]+a[15]*b[15]+a[16]*b[16]-a[17]*b[17]-a[18]*b[18]+a[19]*b[19]-a[20]*b[20]+a[21]*b[21]+a[22]*b[22]-a[23]*b[23]+a[24]*b[24]+a[25]*b[25]+a[26]*b[26]+a[27]*b[27]-a[28]*b[28]-a[29]*b[29]-a[30]*b[30]-a[31]*b[31]-a[32]*b[32]
t[2]=a[2]*b[1]+a[1]*b[2]+a[7]*b[3]+a[8]*b[4]+a[9]*b[5]-a[10]*b[6]-a[3]*b[7]-a[4]*b[8]-a[5]*b[9]+a[6]*b[10]-a[17]*b[11]-a[18]*b[12]+a[19]*b[13]-a[20]*b[14]+a[21]*b[15]+a[22]*b[16]-a[11]*b[17]-a[12]*b[18]+a[13]*b[19]-a[14]*b[20]+a[15]*b[21]+a[16]*b[22]-a[27]*b[23]+a[28]*b[24]+a[29]*b[25]+a[30]*b[26]+a[23]*b[27]-a[24]*b[28]-a[25]*b[29]-a[26]*b[30]-a[32]*b[31]-a[31]*b[32]
t[3]=a[3]*b[1]-a[7]*b[2]+a[1]*b[3]+a[11]*b[4]+a[12]*b[5]-a[13]*b[6]+a[2]*b[7]+a[17]*b[8]+a[18]*b[9]-a[19]*b[10]-a[4]*b[11]-a[5]*b[12]+a[6]*b[13]-a[23]*b[14]+a[24]*b[15]+a[25]*b[16]+a[8]*b[17]+a[9]*b[18]-a[10]*b[19]+a[27]*b[20]-a[28]*b[21]-a[29]*b[22]-a[14]*b[23]+a[15]*b[24]+a[16]*b[25]+a[31]*b[26]-a[20]*b[27]+a[21]*b[28]+a[22]*b[29]+a[32]*b[30]-a[26]*b[31]+a[30]*b[32]
t[4]=a[4]*b[1]-a[8]*b[2]-a[11]*b[3]+a[1]*b[4]+a[14]*b[5]-a[15]*b[6]-a[17]*b[7]+a[2]*b[8]+a[20]*b[9]-a[21]*b[10]+a[3]*b[11]+a[23]*b[12]-a[24]*b[13]-a[5]*b[14]+a[6]*b[15]+a[26]*b[16]-a[7]*b[17]-a[27]*b[18]+a[28]*b[19]+a[9]*b[20]-a[10]*b[21]-a[30]*b[22]+a[12]*b[23]-a[13]*b[24]-a[31]*b[25]+a[16]*b[26]+a[18]*b[27]-a[19]*b[28]-a[32]*b[29]+a[22]*b[30]+a[25]*b[31]-a[29]*b[32]
t[5]=a[5]*b[1]-a[9]*b[2]-a[12]*b[3]-a[14]*b[4]+a[1]*b[5]-a[16]*b[6]-a[18]*b[7]-a[20]*b[8]+a[2]*b[9]-a[22]*b[10]-a[23]*b[11]+a[3]*b[12]-a[25]*b[13]+a[4]*b[14]-a[26]*b[15]+a[6]*b[16]+a[27]*b[17]-a[7]*b[18]+a[29]*b[19]-a[8]*b[20]+a[30]*b[21]-a[10]*b[22]-a[11]*b[23]+a[31]*b[24]-a[13]*b[25]-a[15]*b[26]-a[17]*b[27]+a[32]*b[28]-a[19]*b[29]-a[21]*b[30]-a[24]*b[31]+a[28]*b[32]
t[6]=a[6]*b[1]-a[10]*b[2]-a[13]*b[3]-a[15]*b[4]-a[16]*b[5]+a[1]*b[6]-a[19]*b[7]-a[21]*b[8]-a[22]*b[9]+a[2]*b[10]-a[24]*b[11]-a[25]*b[12]+a[3]*b[13]-a[26]*b[14]+a[4]*b[15]+a[5]*b[16]+a[28]*b[17]+a[29]*b[18]-a[7]*b[19]+a[30]*b[20]-a[8]*b[21]-a[9]*b[22]+a[31]*b[23]-a[11]*b[24]-a[12]*b[25]-a[14]*b[26]+a[32]*b[27]-a[17]*b[28]-a[18]*b[29]-a[20]*b[30]-a[23]*b[31]+a[27]*b[32]
t[7]=a[7]*b[1]-a[3]*b[2]+a[2]*b[3]+a[17]*b[4]+a[18]*b[5]-a[19]*b[6]+a[1]*b[7]+a[11]*b[8]+a[12]*b[9]-a[13]*b[10]-a[8]*b[11]-a[9]*b[12]+a[10]*b[13]-a[27]*b[14]+a[28]*b[15]+a[29]*b[16]+a[4]*b[17]+a[5]*b[18]-a[6]*b[19]+a[23]*b[20]-a[24]*b[21]-a[25]*b[22]-a[20]*b[23]+a[21]*b[24]+a[22]*b[25]+a[32]*b[26]-a[14]*b[27]+a[15]*b[28]+a[16]*b[29]+a[31]*b[30]-a[30]*b[31]+a[26]*b[32]
t[8]=a[8]*b[1]-a[4]*b[2]-a[17]*b[3]+a[2]*b[4]+a[20]*b[5]-a[21]*b[6]-a[11]*b[7]+a[1]*b[8]+a[14]*b[9]-a[15]*b[10]+a[7]*b[11]+a[27]*b[12]-a[28]*b[13]-a[9]*b[14]+a[10]*b[15]+a[30]*b[16]-a[3]*b[17]-a[23]*b[18]+a[24]*b[19]+a[5]*b[20]-a[6]*b[21]-a[26]*b[22]+a[18]*b[23]-a[19]*b[24]-a[32]*b[25]+a[22]*b[26]+a[12]*b[27]-a[13]*b[28]-a[31]*b[29]+a[16]*b[30]+a[29]*b[31]-a[25]*b[32]
t[9]=a[9]*b[1]-a[5]*b[2]-a[18]*b[3]-a[20]*b[4]+a[2]*b[5]-a[22]*b[6]-a[12]*b[7]-a[14]*b[8]+a[1]*b[9]-a[16]*b[10]-a[27]*b[11]+a[7]*b[12]-a[29]*b[13]+a[8]*b[14]-a[30]*b[15]+a[10]*b[16]+a[23]*b[17]-a[3]*b[18]+a[25]*b[19]-a[4]*b[20]+a[26]*b[21]-a[6]*b[22]-a[17]*b[23]+a[32]*b[24]-a[19]*b[25]-a[21]*b[26]-a[11]*b[27]+a[31]*b[28]-a[13]*b[29]-a[15]*b[30]-a[28]*b[31]+a[24]*b[32]
t[10]=a[10]*b[1]-a[6]*b[2]-a[19]*b[3]-a[21]*b[4]-a[22]*b[5]+a[2]*b[6]-a[13]*b[7]-a[15]*b[8]-a[16]*b[9]+a[1]*b[10]-a[28]*b[11]-a[29]*b[12]+a[7]*b[13]-a[30]*b[14]+a[8]*b[15]+a[9]*b[16]+a[24]*b[17]+a[25]*b[18]-a[3]*b[19]+a[26]*b[20]-a[4]*b[21]-a[5]*b[22]+a[32]*b[23]-a[17]*b[24]-a[18]*b[25]-a[20]*b[26]+a[31]*b[27]-a[11]*b[28]-a[12]*b[29]-a[14]*b[30]-a[27]*b[31]+a[23]*b[32]
t[11]=a[11]*b[1]+a[17]*b[2]-a[4]*b[3]+a[3]*b[4]+a[23]*b[5]-a[24]*b[6]+a[8]*b[7]-a[7]*b[8]-a[27]*b[9]+a[28]*b[10]+a[1]*b[11]+a[14]*b[12]-a[15]*b[13]-a[12]*b[14]+a[13]*b[15]+a[31]*b[16]+a[2]*b[17]+a[20]*b[18]-a[21]*b[19]-a[18]*b[20]+a[19]*b[21]+a[32]*b[22]+a[5]*b[23]-a[6]*b[24]-a[26]*b[25]+a[25]*b[26]-a[9]*b[27]+a[10]*b[28]+a[30]*b[29]-a[29]*b[30]+a[16]*b[31]+a[22]*b[32]
t[12]=a[12]*b[1]+a[18]*b[2]-a[5]*b[3]-a[23]*b[4]+a[3]*b[5]-a[25]*b[6]+a[9]*b[7]+a[27]*b[8]-a[7]*b[9]+a[29]*b[10]-a[14]*b[11]+a[1]*b[12]-a[16]*b[13]+a[11]*b[14]-a[31]*b[15]+a[13]*b[16]-a[20]*b[17]+a[2]*b[18]-a[22]*b[19]+a[17]*b[20]-a[32]*b[21]+a[19]*b[22]-a[4]*b[23]+a[26]*b[24]-a[6]*b[25]-a[24]*b[26]+a[8]*b[27]-a[30]*b[28]+a[10]*b[29]+a[28]*b[30]-a[15]*b[31]-a[21]*b[32]
t[13]=a[13]*b[1]+a[19]*b[2]-a[6]*b[3]-a[24]*b[4]-a[25]*b[5]+a[3]*b[6]+a[10]*b[7]+a[28]*b[8]+a[29]*b[9]-a[7]*b[10]-a[15]*b[11]-a[16]*b[12]+a[1]*b[13]-a[31]*b[14]+a[11]*b[15]+a[12]*b[16]-a[21]*b[17]-a[22]*b[18]+a[2]*b[19]-a[32]*b[20]+a[17]*b[21]+a[18]*b[22]+a[26]*b[23]-a[4]*b[24]-a[5]*b[25]-a[23]*b[26]-a[30]*b[27]+a[8]*b[28]+a[9]*b[29]+a[27]*b[30]-a[14]*b[31]-a[20]*b[32]
t[14]=a[14]*b[1]+a[20]*b[2]+a[23]*b[3]-a[5]*b[4]+a[4]*b[5]-a[26]*b[6]-a[27]*b[7]+a[9]*b[8]-a[8]*b[9]+a[30]*b[10]+a[12]*b[11]-a[11]*b[12]+a[31]*b[13]+a[1]*b[14]-a[16]*b[15]+a[15]*b[16]+a[18]*b[17]-a[17]*b[18]+a[32]*b[19]+a[2]*b[20]-a[22]*b[21]+a[21]*b[22]+a[3]*b[23]-a[25]*b[24]+a[24]*b[25]-a[6]*b[26]-a[7]*b[27]+a[29]*b[28]-a[28]*b[29]+a[10]*b[30]+a[13]*b[31]+a[19]*b[32]
t[15]=a[15]*b[1]+a[21]*b[2]+a[24]*b[3]-a[6]*b[4]-a[26]*b[5]+a[4]*b[6]-a[28]*b[7]+a[10]*b[8]+a[30]*b[9]-a[8]*b[10]+a[13]*b[11]+a[31]*b[12]-a[11]*b[13]-a[16]*b[14]+a[1]*b[15]+a[14]*b[16]+a[19]*b[17]+a[32]*b[18]-a[17]*b[19]-a[22]*b[20]+a[2]*b[21]+a[20]*b[22]-a[25]*b[23]+a[3]*b[24]+a[23]*b[25]-a[5]*b[26]+a[29]*b[27]-a[7]*b[28]-a[27]*b[29]+a[9]*b[30]+a[12]*b[31]+a[18]*b[32]
t[16]=a[16]*b[1]+a[22]*b[2]+a[25]*b[3]+a[26]*b[4]-a[6]*b[5]+a[5]*b[6]-a[29]*b[7]-a[30]*b[8]+a[10]*b[9]-a[9]*b[10]-a[31]*b[11]+a[13]*b[12]-a[12]*b[13]+a[15]*b[14]-a[14]*b[15]+a[1]*b[16]-a[32]*b[17]+a[19]*b[18]-a[18]*b[19]+a[21]*b[20]-a[20]*b[21]+a[2]*b[22]+a[24]*b[23]-a[23]*b[24]+a[3]*b[25]+a[4]*b[26]-a[28]*b[27]+a[27]*b[28]-a[7]*b[29]-a[8]*b[30]-a[11]*b[31]-a[17]*b[32]
t[17]=a[17]*b[1]+a[11]*b[2]-a[8]*b[3]+a[7]*b[4]+a[27]*b[5]-a[28]*b[6]+a[4]*b[7]-a[3]*b[8]-a[23]*b[9]+a[24]*b[10]+a[2]*b[11]+a[20]*b[12]-a[21]*b[13]-a[18]*b[14]+a[19]*b[15]+a[32]*b[16]+a[1]*b[17]+a[14]*b[18]-a[15]*b[19]-a[12]*b[20]+a[13]*b[21]+a[31]*b[22]+a[9]*b[23]-a[10]*b[24]-a[30]*b[25]+a[29]*b[26]-a[5]*b[27]+a[6]*b[28]+a[26]*b[29]-a[25]*b[30]+a[22]*b[31]+a[16]*b[32]
t[18]=a[18]*b[1]+a[12]*b[2]-a[9]*b[3]-a[27]*b[4]+a[7]*b[5]-a[29]*b[6]+a[5]*b[7]+a[23]*b[8]-a[3]*b[9]+a[25]*b[10]-a[20]*b[11]+a[2]*b[12]-a[22]*b[13]+a[17]*b[14]-a[32]*b[15]+a[19]*b[16]-a[14]*b[17]+a[1]*b[18]-a[16]*b[19]+a[11]*b[20]-a[31]*b[21]+a[13]*b[22]-a[8]*b[23]+a[30]*b[24]-a[10]*b[25]-a[28]*b[26]+a[4]*b[27]-a[26]*b[28]+a[6]*b[29]+a[24]*b[30]-a[21]*b[31]-a[15]*b[32]
t[19]=a[19]*b[1]+a[13]*b[2]-a[10]*b[3]-a[28]*b[4]-a[29]*b[5]+a[7]*b[6]+a[6]*b[7]+a[24]*b[8]+a[25]*b[9]-a[3]*b[10]-a[21]*b[11]-a[22]*b[12]+a[2]*b[13]-a[32]*b[14]+a[17]*b[15]+a[18]*b[16]-a[15]*b[17]-a[16]*b[18]+a[1]*b[19]-a[31]*b[20]+a[11]*b[21]+a[12]*b[22]+a[30]*b[23]-a[8]*b[24]-a[9]*b[25]-a[27]*b[26]-a[26]*b[27]+a[4]*b[28]+a[5]*b[29]+a[23]*b[30]-a[20]*b[31]-a[14]*b[32]
t[20]=a[20]*b[1]+a[14]*b[2]+a[27]*b[3]-a[9]*b[4]+a[8]*b[5]-a[30]*b[6]-a[23]*b[7]+a[5]*b[8]-a[4]*b[9]+a[26]*b[10]+a[18]*b[11]-a[17]*b[12]+a[32]*b[13]+a[2]*b[14]-a[22]*b[15]+a[21]*b[16]+a[12]*b[17]-a[11]*b[18]+a[31]*b[19]+a[1]*b[20]-a[16]*b[21]+a[15]*b[22]+a[7]*b[23]-a[29]*b[24]+a[28]*b[25]-a[10]*b[26]-a[3]*b[27]+a[25]*b[28]-a[24]*b[29]+a[6]*b[30]+a[19]*b[31]+a[13]*b[32]
t[21]=a[21]*b[1]+a[15]*b[2]+a[28]*b[3]-a[10]*b[4]-a[30]*b[5]+a[8]*b[6]-a[24]*b[7]+a[6]*b[8]+a[26]*b[9]-a[4]*b[10]+a[19]*b[11]+a[32]*b[12]-a[17]*b[13]-a[22]*b[14]+a[2]*b[15]+a[20]*b[16]+a[13]*b[17]+a[31]*b[18]-a[11]*b[19]-a[16]*b[20]+a[1]*b[21]+a[14]*b[22]-a[29]*b[23]+a[7]*b[24]+a[27]*b[25]-a[9]*b[26]+a[25]*b[27]-a[3]*b[28]-a[23]*b[29]+a[5]*b[30]+a[18]*b[31]+a[12]*b[32]
t[22]=a[22]*b[1]+a[16]*b[2]+a[29]*b[3]+a[30]*b[4]-a[10]*b[5]+a[9]*b[6]-a[25]*b[7]-a[26]*b[8]+a[6]*b[9]-a[5]*b[10]-a[32]*b[11]+a[19]*b[12]-a[18]*b[13]+a[21]*b[14]-a[20]*b[15]+a[2]*b[16]-a[31]*b[17]+a[13]*b[18]-a[12]*b[19]+a[15]*b[20]-a[14]*b[21]+a[1]*b[22]+a[28]*b[23]-a[27]*b[24]+a[7]*b[25]+a[8]*b[26]-a[24]*b[27]+a[23]*b[28]-a[3]*b[29]-a[4]*b[30]-a[17]*b[31]-a[11]*b[32]
t[23]=a[23]*b[1]-a[27]*b[2]+a[14]*b[3]-a[12]*b[4]+a[11]*b[5]-a[31]*b[6]+a[20]*b[7]-a[18]*b[8]+a[17]*b[9]-a[32]*b[10]+a[5]*b[11]-a[4]*b[12]+a[26]*b[13]+a[3]*b[14]-a[25]*b[15]+a[24]*b[16]-a[9]*b[17]+a[8]*b[18]-a[30]*b[19]-a[7]*b[20]+a[29]*b[21]-a[28]*b[22]+a[1]*b[23]-a[16]*b[24]+a[15]*b[25]-a[13]*b[26]+a[2]*b[27]-a[22]*b[28]+a[21]*b[29]-a[19]*b[30]+a[6]*b[31]-a[10]*b[32]
t[24]=a[24]*b[1]-a[28]*b[2]+a[15]*b[3]-a[13]*b[4]-a[31]*b[5]+a[11]*b[6]+a[21]*b[7]-a[19]*b[8]-a[32]*b[9]+a[17]*b[10]+a[6]*b[11]+a[26]*b[12]-a[4]*b[13]-a[25]*b[14]+a[3]*b[15]+a[23]*b[16]-a[10]*b[17]-a[30]*b[18]+a[8]*b[19]+a[29]*b[20]-a[7]*b[21]-a[27]*b[22]-a[16]*b[23]+a[1]*b[24]+a[14]*b[25]-a[12]*b[26]-a[22]*b[27]+a[2]*b[28]+a[20]*b[29]-a[18]*b[30]+a[5]*b[31]-a[9]*b[32]
t[25]=a[25]*b[1]-a[29]*b[2]+a[16]*b[3]+a[31]*b[4]-a[13]*b[5]+a[12]*b[6]+a[22]*b[7]+a[32]*b[8]-a[19]*b[9]+a[18]*b[10]-a[26]*b[11]+a[6]*b[12]-a[5]*b[13]+a[24]*b[14]-a[23]*b[15]+a[3]*b[16]+a[30]*b[17]-a[10]*b[18]+a[9]*b[19]-a[28]*b[20]+a[27]*b[21]-a[7]*b[22]+a[15]*b[23]-a[14]*b[24]+a[1]*b[25]+a[11]*b[26]+a[21]*b[27]-a[20]*b[28]+a[2]*b[29]+a[17]*b[30]-a[4]*b[31]+a[8]*b[32]
t[26]=a[26]*b[1]-a[30]*b[2]-a[31]*b[3]+a[16]*b[4]-a[15]*b[5]+a[14]*b[6]-a[32]*b[7]+a[22]*b[8]-a[21]*b[9]+a[20]*b[10]+a[25]*b[11]-a[24]*b[12]+a[23]*b[13]+a[6]*b[14]-a[5]*b[15]+a[4]*b[16]-a[29]*b[17]+a[28]*b[18]-a[27]*b[19]-a[10]*b[20]+a[9]*b[21]-a[8]*b[22]-a[13]*b[23]+a[12]*b[24]-a[11]*b[25]+a[1]*b[26]-a[19]*b[27]+a[18]*b[28]-a[17]*b[29]+a[2]*b[30]+a[3]*b[31]-a[7]*b[32]
t[27]=a[27]*b[1]-a[23]*b[2]+a[20]*b[3]-a[18]*b[4]+a[17]*b[5]-a[32]*b[6]+a[14]*b[7]-a[12]*b[8]+a[11]*b[9]-a[31]*b[10]+a[9]*b[11]-a[8]*b[12]+a[30]*b[13]+a[7]*b[14]-a[29]*b[15]+a[28]*b[16]-a[5]*b[17]+a[4]*b[18]-a[26]*b[19]-a[3]*b[20]+a[25]*b[21]-a[24]*b[22]+a[2]*b[23]-a[22]*b[24]+a[21]*b[25]-a[19]*b[26]+a[1]*b[27]-a[16]*b[28]+a[15]*b[29]-a[13]*b[30]+a[10]*b[31]-a[6]*b[32]
t[28]=a[28]*b[1]-a[24]*b[2]+a[21]*b[3]-a[19]*b[4]-a[32]*b[5]+a[17]*b[6]+a[15]*b[7]-a[13]*b[8]-a[31]*b[9]+a[11]*b[10]+a[10]*b[11]+a[30]*b[12]-a[8]*b[13]-a[29]*b[14]+a[7]*b[15]+a[27]*b[16]-a[6]*b[17]-a[26]*b[18]+a[4]*b[19]+a[25]*b[20]-a[3]*b[21]-a[23]*b[22]-a[22]*b[23]+a[2]*b[24]+a[20]*b[25]-a[18]*b[26]-a[16]*b[27]+a[1]*b[28]+a[14]*b[29]-a[12]*b[30]+a[9]*b[31]-a[5]*b[32]
t[29]=a[29]*b[1]-a[25]*b[2]+a[22]*b[3]+a[32]*b[4]-a[19]*b[5]+a[18]*b[6]+a[16]*b[7]+a[31]*b[8]-a[13]*b[9]+a[12]*b[10]-a[30]*b[11]+a[10]*b[12]-a[9]*b[13]+a[28]*b[14]-a[27]*b[15]+a[7]*b[16]+a[26]*b[17]-a[6]*b[18]+a[5]*b[19]-a[24]*b[20]+a[23]*b[21]-a[3]*b[22]+a[21]*b[23]-a[20]*b[24]+a[2]*b[25]+a[17]*b[26]+a[15]*b[27]-a[14]*b[28]+a[1]*b[29]+a[11]*b[30]-a[8]*b[31]+a[4]*b[32]
t[30]=a[30]*b[1]-a[26]*b[2]-a[32]*b[3]+a[22]*b[4]-a[21]*b[5]+a[20]*b[6]-a[31]*b[7]+a[16]*b[8]-a[15]*b[9]+a[14]*b[10]+a[29]*b[11]-a[28]*b[12]+a[27]*b[13]+a[10]*b[14]-a[9]*b[15]+a[8]*b[16]-a[25]*b[17]+a[24]*b[18]-a[23]*b[19]-a[6]*b[20]+a[5]*b[21]-a[4]*b[22]-a[19]*b[23]+a[18]*b[24]-a[17]*b[25]+a[2]*b[26]-a[13]*b[27]+a[12]*b[28]-a[11]*b[29]+a[1]*b[30]+a[7]*b[31]-a[3]*b[32]
t[31]=a[31]*b[1]+a[32]*b[2]-a[26]*b[3]+a[25]*b[4]-a[24]*b[5]+a[23]*b[6]+a[30]*b[7]-a[29]*b[8]+a[28]*b[9]-a[27]*b[10]+a[16]*b[11]-a[15]*b[12]+a[14]*b[13]+a[13]*b[14]-a[12]*b[15]+a[11]*b[16]+a[22]*b[17]-a[21]*b[18]+a[20]*b[19]+a[19]*b[20]-a[18]*b[21]+a[17]*b[22]-a[6]*b[23]+a[5]*b[24]-a[4]*b[25]+a[3]*b[26]+a[10]*b[27]-a[9]*b[28]+a[8]*b[29]-a[7]*b[30]+a[1]*b[31]+a[2]*b[32]
t[32]=a[32]*b[1]+a[31]*b[2]-a[30]*b[3]+a[29]*b[4]-a[28]*b[5]+a[27]*b[6]+a[26]*b[7]-a[25]*b[8]+a[24]*b[9]-a[23]*b[10]+a[22]*b[11]-a[21]*b[12]+a[20]*b[13]+a[19]*b[14]-a[18]*b[15]+a[17]*b[16]+a[16]*b[17]-a[15]*b[18]+a[14]*b[19]+a[13]*b[20]-a[12]*b[21]+a[11]*b[22]-a[10]*b[23]+a[9]*b[24]-a[8]*b[25]+a[7]*b[26]+a[6]*b[27]-a[5]*b[28]+a[4]*b[29]-a[3]*b[30]+a[2]*b[31]+a[1]*b[32]
end

sig: wedge4
function(t,a,b)
t[1]=a[1]*b[1]
t[2]=a[2]*b[1]+a[1]*b[2]
t[3]=a[3]*b[1]+a[1]*b[3]
t[4]=a[4]*b[1]+a[1]*b[4]
t[5]=a[5]*b[1]+a[1]*b[5]
t[6]=a[6]*b[1]-a[3]*b[2]+a[2]*b[3]+a[1]*b[6]
t[7]=a[7]*b[1]-a[4]*b[2]+a[2]*b[4]+a[1]*b[7]
t[8]=a[8]*b[1]-a[5]*b[2]+a[2]*b[5]+a[1]*b[8]
t[9]=a[9]*b[1]-a[4]*b[3]+a[3]*b[4]+a[1]*b[9]
t[10]=a[10]*b[1]-a[5]*b[3]+a[3]*b[5]+a[1]*b[10]
t[11]=a[11]*b[1]-a[5]*b[4]+a[4]*b[5]+a[1]*b[11]
t[12]=a[12]*b[1]+a[9]*b[2]-a[7]*b[3]+a[6]*b[4]+a[4]*b[6]-a[3]*b[7]+a[2]*b[9]+a[1]*b[12]
t[13]=a[13]*b[1]+a[10]*b[2]-a[8]*b[3]+a[6]*b[5]+a[5]*b[6]-a[3]*b[8]+a[2]*b[10]+a[1]*b[13]
t[14]=a[14]*b[1]+a[11]*b[2]-a[8]*b[4]+a[7]*b[5]+a[5]*b[7]-a[4]*b[8]+a[2]*b[11]+a[1]*b[14]
t[15]=a[15]*b[1]+a[11]*b[3]-a[10]*b[4]+a[9]*b[5]+a[5]*b[9]-a[4]*b[10]+a[3]*b[11]+a[1]*b[15]
t[16]=a[16]*b[1]-a[15]*b[2]+a[14]*b[3]-a[13]*b[4]+a[12]*b[5]+a[11]*b[6]-a[10]*b[7]+a[9]*b[8]+a[8]*b[9]-a[7]*b[10]+a[6]*b[11]-a[5]*b[12]+a[4]*b[13]-a[3]*b[14]+a[2]*b[15]+a[1]*b[16]
end

sig: quat
function(t,a,b)
t[1]=a[1]*b[1]-a[2]*b[2]-a[3]*b[3]-a[4]*b[4] -- 1
t[2]=a[2]*b[1]+a[1]*b[2]-a[4]*b[3]+a[3]*b[4] -- k
t[3]=a[3]*b[1]+a[4]*b[2]+a[1]*b[3]-a[2]*b[4] -- j
t[4]=a[4]*b[1]-a[3]*b[2]+a[2]*b[3]+a[1]*b[4] -- i
end


]=]