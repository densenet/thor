local D_TEST = 4

local b, x, l = bit32.btest, bit32.bxor, bit32.lshift

local function inc_word(w: number, D: number)
  -- want to swap most significant (0,1) pair, as long as there's no bit at Dth bit.
  -- if at Dth bit:
  --  * slide Dth bit and all adjacent on-bits up to first (0,1) pair
  --  * swap aforementioned (0,1) pair
  --  if no (0,1) pair exists, slide to almost-LSB & turn on LSB

  if b(w, l(1, D - 1)) then
   
  end
  
 
end

local random_word = math.pow(2, D_TEST)