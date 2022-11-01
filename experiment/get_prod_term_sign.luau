local bit2str = require("src\\util\\bit2str")

local b, _x, l = bit32.btest, bit32.bxor, bit32.lshift

type Signature = {
  Unit: number,
  Anti: number,
  Null: number,
}

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

local bw1 = math.random(0, 15)
local bw2 = math.random(0, 15)

print(string.format("bw1: %s; bw2: %s", bit2str(bw1), bit2str(bw2)))

local pga3d: Signature = { Unit = 3, Anti = 0, Null = 1 }

print(string.format("get_prod_term_sign(pga3d): %d", get_prod_term_sign(pga3d, bw1, bw2)))

--[=[
bw1: 00000111; bw2: 00000110
get_prod_term_sign(pga3d): -1
bw1: 00000001; bw2: 00001010
get_prod_term_sign(pga3d): 1
bw1: 00000010; bw2: 00001000
get_prod_term_sign(pga3d): 1
bw1: 00000000; bw2: 00001101
get_prod_term_sign(pga3d): 1
bw1: 00001111; bw2: 00000010
get_prod_term_sign(pga3d): 1
bw1: 00001000; bw2: 00000100
get_prod_term_sign(pga3d): -1
bw1: 00001000; bw2: 00001111
get_prod_term_sign(pga3d): -1
bw1: 00000101; bw2: 00000100
get_prod_term_sign(pga3d): 1
bw1: 00000001; bw2: 00000001
get_prod_term_sign(pga3d): 0

]=]
