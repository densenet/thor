local D = 4

local BINC = require("src\\util\\binc")(D)
local inc_word = require("src\\util\\inc_word")

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
  Unit: dim,
  Anti: dim,
  Null: dim
}
local _b, _x, l, n = bit32.btest, bit32.bxor, bit32.lshift, bit32.bnot

local get_sigmaps
do
  local function generate_sigmap_tables(D: dim): any
    local gi2i = table.create(D + 2)
    local fr = table.freeze
    
    print(D)
    
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

local pga3d = {Unit = 3, Anti = 0, Null = 1}


for k, v in get_sigmaps(pga3d) do
  print(k)
  for k2, v2 in v do
    print("-", k2, v2)
  end    
end



--[=[
for k, v in BINC do
  print(k)
  for k2, v2 in v do
    print("-", k2, v2)
  end
end
]=]


--[=[
4
GRADE_IDX_TO_FIRST_IDX
-	1	1
-	2	2
-	3	6
-	4	12
-	5	16
-	6	17
IDX_TO_GRADE_WORD
-	1	1
-	2	2
-	3	2
-	4	2
-	5	2
-	6	4
-	7	4
-	8	4
-	9	4
-	10	4
-	11	4
-	12	8
-	13	8
-	14	8
-	15	8
-	16	16
WORD_TO_IDX
-	1	2
-	2	3
-	3	6
-	4	4
-	5	7
-	6	9
-	7	12
-	8	5
-	9	8
-	10	10
-	11	13
-	12	11
-	13	14
-	14	15
-	15	16
-	0	1
IDX_TO_WORD
-	1	0
-	2	1
-	3	2
-	4	4
-	5	8
-	6	3
-	7	5
-	8	9
-	9	6
-	10	10
-	11	12
-	12	7
-	13	11
-	14	13
-	15	14
-	16	15


]=]

--[=[
1
-	1	1
-	2	0
-	3	0
-	4	0
-	5	0
2
-	1	1
-	2	1
-	3	0
-	4	0
-	5	0
3
-	1	1
-	2	2
-	3	1
-	4	0
-	5	0
4
-	1	1
-	2	3
-	3	3
-	4	1
-	5	0
5
-	1	1
-	2	4
-	3	6
-	4	4
-	5	1



GRADE_IDX_TO_FIRST_IDX
-	1	1
-	2	2
-	3	5
-	4	8
-	5	9
-	6	9
IDX_TO_GRADE_WORD
-	1	2
-	2	2
-	3	4
WORD_TO_IDX
-	1	2
-	2	3
-	3	3
-	4	4
-	5	4
-	7	4
-	9	5
-	0	1
IDX_TO_WORD
-	1	0
-	2	1
-	3	3
-	4	7
-	5	9

]=]
