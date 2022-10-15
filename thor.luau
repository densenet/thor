-- polychromatist

-- see mv-primer.txt for information

-- note that two basis elements juxtaposed together is a wedge product between them
local basis_index = {
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

local basis_list = {
  0b0000, -- scalar component
  0b0001, 0b0010, 0b0100, 0b1000, --  vector components
  0b0011, 0b0101, 0b1001, 0b0110, 0b1010, 0b1100, -- bi-vector components
  0b0111, 0b1011, 0b1101, 0b1110, -- tri-vector / pseudovector components
  0b1111 -- pseudoscalar component
}

local basis_tokens = {
	"_", "w", "x", "y", "z", "wx", "wy", "wz", "xy", "xz", "yz", "wxy", "wxz", "wyz", "xyz", "wxyz"
}

local btest = bit32.btest
local EPSREAL = (function()
  -- we must compute a certain count value from our index set.
  -- this count value directly can determine the sign of a product between the two component bases
  -- if we have a map index (DCBA, HGFE) where each letter represents a bit, then
  -- the map value is as follows:
  --  value: 0
  --  for each on bit X of DCBA
  --    store position in P (i.e. A is 1, B is 2, etc)
  --    for 1 <= i < P
  --      if HGFE[P] is on
  --        value ++;
  -- what the map value represents is, if you can picture both bases juxtaposed to each other,
  -- how many swaps between distinct wedged vector component basis pairs is required until we reduce back down.
  -- for example: (wxz)(xy)
  -- this corresponds to map index (0b1011, 0b0110)
  -- there are 3 on bits for DCBA: bit 1, bit 2, and bit 4
  -- now we iterate all relatively rightward bit positions on the right basis, 3 times
  -- - there are never any bits right of bit 1 (value = 0)
  -- - there are no bits right of bit 2 in this case, bit 1 of right basis is zero (value = 0)
  -- - there are two bits right of bit 4 (value = 2)
  -- so there are at least two swaps that must be done, and indeed,
  -- +wxzxy -> -wxxzy -> -wzy -> +wyz
  -- note that two identical wedged vector component basis pairs annihilate each other when juxtaposed: xx = 1.
  -- since the value is even, the EPSREAL sign value will be positive. if it were odd, it would be negative.

  local ret = {}
  for i = 0, 255 do
    local value = 0
    if btest(i, 0b00000010) then
      if btest(i, 0b00010000) then
        value += 1
      end
    end
    if btest(i, 0b00000100) then
      if btest(i, 0b00010000) then
        value += 1
      end
      if btest(i, 0b00100000) then
        value += 1
      end
    end
    if btest(i, 0b00001000) then
      if btest(i, 0b00010000) then
        value += 1
      end
      if btest(i, 0b00100000) then
        value += 1
      end
      if btest(i, 0b01000000) then
        value += 1
      end
    end
    if bit32.band(value, 1) == 1 then
      ret[i] = -1
    else
      ret[i] = 1
    end
  end

  -- to get index of EPSREAL with two basis nibbles:
  -- bxor(lshift(right_basis, 4), left_basis)
  -- could it be a 2D array instead? yes. but as for now, it's not

  return ret
end)()

-- template for multiplication
--[=[
local i_1
for i = 1, 16 do
  i_1 = basis_list[i]
  t_prod[i] += t1[i] * t2[1]
    + EPSREAL[bxor(i_1, 0x10)] * (t1[basis_index[bxor(i_1, 0x1)]] * t2[2])  + EPSREAL[bxor(i_1, 0x20)] * (t1[basis_index[bxor(i_1, 0x2)]] * t2[3])  + EPSREAL[bxor(i_1, 0x40)] * (t1[basis_index[bxor(i_1, 0x4)]] * t2[4])  + EPSREAL[bxor(i_1, 0x80)] * (t1[basis_index[bxor(i_1, 0x8)]] * t2[5])
    + EPSREAL[bxor(i_1, 0x30)] * (t1[basis_index[bxor(i_1, 0x3)]] * t2[6])  + EPSREAL[bxor(i_1, 0x50)] * (t1[basis_index[bxor(i_1, 0x5)]] * t2[7])  + EPSREAL[bxor(i_1, 0x90)] * (t1[basis_index[bxor(i_1, 0x9)]] * t2[8])  + EPSREAL[bxor(i_1, 0x60)] * (t1[basis_index[bxor(i_1, 0x6)]] * t2[9])  + EPSREAL[bxor(i_1, 0xA0)] * (t1[basis_index[bxor(i_1, 0xA)]] * t2[10]) + EPSREAL[bxor(i_1, 0xC0)] * (t1[basis_index[bxor(i_1, 0xC)]] * t2[11])
    + EPSREAL[bxor(i_1, 0x70)] * (t1[basis_index[bxor(i_1, 0x7)]] * t2[12]) + EPSREAL[bxor(i_1, 0xB0)] * (t1[basis_index[bxor(i_1, 0xB)]] * t2[13]) + EPSREAL[bxor(i_1, 0xD0)] * (t1[basis_index[bxor(i_1, 0xD)]] * t2[14]) + EPSREAL[bxor(i_1, 0xE0)] * (t1[basis_index[bxor(i_1, 0xE)]] * t2[15])
    + EPSREAL[bxor(i_1, 0xF0)] * (t1[basis_index[bxor(i_1, 0xF)]] * t2[16])
end
]=]


local bxor = bit32.bxor
local abs = math.abs

local thor
local Thor = {}
Thor.unary = {}
Thor.unary.unpack = table.unpack
Thor.unary._ = function(t) return t[1] end
Thor.unary.w = function(t) return t[2] end
Thor.unary.x = function(t) return t[3] end
Thor.unary.y = function(t) return t[4] end
Thor.unary.z = function(t) return t[5] end
Thor.unary.wx = function(t) return t[6] end
Thor.unary.wy = function(t) return t[7] end
Thor.unary.wz = function(t) return t[8] end
Thor.unary.xy = function(t) return t[9] end
Thor.unary.xz = function(t) return t[10] end
Thor.unary.yz = function(t) return t[11] end
Thor.unary.wxy = function(t) return t[12] end
Thor.unary.wxz = function(t) return t[13] end
Thor.unary.wyz = function(t) return t[14] end
Thor.unary.xyz = function(t) return t[15] end
Thor.unary.wxyz = function(t) return t[16] end
Thor.unary.v = function(t) return t[2], t[3], t[4], t[5] end
Thor.unary.bv = function(t) return t[6], t[7], t[8], t[9], t[10], t[11] end
Thor.unary.pv = function(t) return t[12], t[13], t[14], t[15] end
function Thor.unary.step(t)
  if t[16] > 1e-10 or t[16] < -1e-10 then
    return 4
  elseif t[15] > 1e-10 or t[15] < -1e-10 or t[14] > 1e-10 or t[14] < -1e-10 or t[13] > 1e-10 or t[13] < -1e-10 or t[12] > 1e-10 or t[12] < -1e-10 then
    return 3
  elseif t[11] > 1e-10 or t[11] < -1e-10 or t[10] > 1e-10 or t[10] < -1e-10 or t[9] > 1e-10 or t[9] < -1e-10 or t[8] > 1e-10 or t[8] < -1e-10 or t[7] > 1e-10 or t[7] < -1e-10 or t[6] > 1e-10 or t[6] < -1e-10 then
    return 2
  elseif t[5] > 1e-10 or t[5] < -1e-10 or t[4] > 1e-10 or t[4] < -1e-10 or t[3] > 1e-10 or t[3] < -1e-10 or t[2] > 1e-10 or t[2] < -1e-10 then
    return 1
  end
  return 0
end

local kernel_static_map = {}

-- the gate is a 5 bit word representing the grades that will be used in the multiplication.
local function getComputeKernel(gate)
  if kernel_static_map[gate] then return kernel_static_map[gate] end

  local K = {}
  local k

  K.unary = table.clone(Thor.unary)
  K.static = {kernel = gate}
  K.lib = {}

  function K.unary.scalar(t)
    local new_t = table.create(16, 0)
    new_t[1] = t[1]
    return setmetatable(t, K)
  end
  function K.unary.vector(t)
    local new_t = table.create(16, 0)
    new_t[2], new_t[3], new_t[4], new_t[5] = t[2], t[3], t[4], t[5]
    return setmetatable(t, K)
  end
  function K.unary.bivector(t)
    local new_t = table.create(16, 0)
    new_t[6], new_t[7], new_t[8], new_t[9], new_t[10], new_t[11] = t[6], t[7], t[8], t[9], t[10], t[11]
    return setmetatable(t, K)
  end
  function K.unary.trivector(t)
    local new_t = table.create(16, 0)
    new_t[12], new_t[13], new_t[14], new_t[15] = t[12], t[13], t[14], t[15]
    return setmetatable(t, K)
  end
  function K.unary.pseudovector(t)
    local new_t = table.create(16, 0)
    new_t[16] = t[16]
    return setmetatable(t, K)
  end

    -- for intuition about multivector-context operations,
    -- see reference: Functions of Multivector Variables, doi: 10.1371/journal.pone.0116943

  -- REVERSION: for each blade A_k represented by a grade projection of k on a multivector M,
  -- A_k = pow(-1, k*(k - 1) / 2) * A_k
  function K.unary.rev(t)
    local new_t = table.clone(t)
    new_t[6], new_t[7], new_t[8], new_t[9], new_t[10], new_t[11] = -new_t[6], -new_t[7], -new_t[8], -new_t[9], -new_t[10], -new_t[11]
    new_t[12], new_t[13], new_t[14], new_t[15] = -new_t[12], -new_t[13], -new_t[14], -new_t[15]
    return setmetatable(t, K)
  end
    -- GRADE INVOLUTION: for each grade A_k,
    -- A_k = pow(-1, k) * A_k
  function K.unary.caret(t)
    local new_t = table.clone(t)
    new_t[2], new_t[3], new_t[4], new_t[5] = -new_t[2], -new_t[3], -new_t[4], -new_t[5]
    new_t[12], new_t[13], new_t[14], new_t[15] = -new_t[12], -new_t[13], -new_t[14], -new_t[15]
    return setmetatable(t, K)
  end
  -- CLIFFORD CONJUGATION
  -- A_k  = pow(-1, k*(k+1) / 2) * A_k
  --      = rev(caret(A_k))
  --      = caret(rev(A_k))
  function K.unary.bar(t)
    local new_t = table.clone(t)
    new_t[2], new_t[3], new_t[4], new_t[5] = -new_t[2], -new_t[3], -new_t[4], -new_t[5]
    new_t[6], new_t[7], new_t[8], new_t[9], new_t[10], new_t[11] = -new_t[6], -new_t[7], -new_t[8], -new_t[9], -new_t[10], -new_t[11]
    return setmetatable(t, K)
  end

  function K.unary.sq(t)
    return t * t
  end

  function K.unary.sqrt(t)
    -- stubbed
  end

  function K.unary.s_mag(t)
    return math.sqrt(t.rev:scalarprod(t))
  end

    -- INVERSE
    -- sources
    -- https://arxiv.org/pdf/1712.05204.pdf (p. 6)
    -- https://core.ac.uk/download/pdf/74374477.pdf (p. 5,6)
  function K.unary.inv(t)
    local t_bar = t.bar
    local b = t.bar * t
    local b_bar = b:minusgate(0b11000)
    
    local s = (b_bar * b)._

    if math.abs(s) < 1e-9 then
      return Thor.lib.zero
    end

    return (1 / s) * b_bar * t_bar
  end

  return setmetatable({}, {
    __index = K.static,
    __call = k
  })
end

thor = getComputeKernel(0x1f)

Thor = {
	lib = {
    --[=[
		inner = function(t1, t2)
			return 0.5 * (t1 * t2 + t2.reversion * t1.reversion)
		end,
		outer = function(t1, t2)
			return 0.5 * (t1 * t2 - t2.reversion * t1.reversion)
		end,
    ]=]
		exterior = function(t1, t2)
      --[=[
			local prod = thor()
			
			local t1r, t2s, term_rs
			for r = 0, 4 do
				t1r = thor(t1:gradeproject(r))
				if t1r == Thor.lib.zero then
					continue
				end
				for s = 0, 4 do
					t2s = thor(t2:gradeproject(s))
					if r + s > 4 or t2s == Thor.lib.zero then
						continue
					end
					
					prod += (t1r * t2s):gradeproject(r + s)
				end
			end
      ]=]
			
			return t1 ^ t2
		end,
		-- the scalar part of the geometric product
		-- note that sometimes "the scalar product" refers to the regular inner product
		-- also note the following identity:
		-- t1:scalarprod(t2) + t1:fatdot(t2) = t1:left(t2) + t1:right(t2)
		scalarprod = function(t1, t2)
			return (t1 * t2).scalar
		end,
		fatdot = function(t1, t2)
			local prod = thor()

			local t1r, t2s, term_rs
			for r = 0, 4 do
				t1r = thor(t1:gradeproject(r))
        if t1r == Thor.lib.zero then
          continue
        end
				for s = 0, 4 do
					t2s = thor(t2:gradeproject(s))
          if t2s == Thor.lib.zero then
            continue
          end

					term_rs = (t1r * t2s):gradeproject(math.abs(s - r))

					if term_rs then
						prod += term_rs
					end

					term_rs = nil
				end
			end

			return prod
		end,
		left = function(t1, t2)
			local prod = thor()

			local t1r, t2s
			for r = 0, 4 do
				t1r = thor(t1:gradeproject(r))
        if t1r == Thor.lib.zero then
          continue
        end
				for s = r, 4 do
					t2s = thor(t2:gradeproject(s))
          if t2s == Thor.lib.zero then
            continue
          end

					prod += (t1r * t2s):gradeproject(s - r)
				end
			end

			return prod
		end,
		right = function(t1, t2)
			local prod = thor()

			local t1r, t2s
      for s = 0, 4 do
        local t2s = thor(t2:gradeproject(s))
        if t2s == Thor.lib.zero then
          continue
        end
        for r = s, 4 do
          t1r = thor(t1:gradeproject(r))
          if t1r == Thor.lib.zero then
            continue
          end

          prod += (t1r * t2s):gradeproject(s - r)
        end
      end

			return prod
		end,
		regressive = function(t1, t2)
      return t1 .. t2
		end,
		project = function(t1, t2)
			return (t1 < t2) * t2.linverse
		end,
		reject = function(t1, t2)
			return (t1 ^ t2) * t2.linverse
		end,
		reflect = function(t1, t2)
			return -t2 * t1 * t2.linverse
		end,
		-- get a multivector's k-th grade component
    -- for multigrade projection use zerogate
		gradeproject = function(t1, k)
			if k == 0 then
				return t1.scalar
			elseif k == 1 then
				return t1.vector
			elseif k == 2 then
				return t1.bivector
      elseif k == 3 then
        return t1.trivector
			elseif k == 4 then
				return t1.pscalar
      end
      return Thor.lib.zero
		end,
    -- the following functions, zerogate and minusgate, are called grade-specific maps
    -- gb is a 5 bit integer of which each bit position corresponds to a multivector grade
    -- if any bit positions are on, then those grades will be set to zero
    zerogate = function(t1, gb)
      if gb == 0 then
        return t1
      end

      local new_t1 = table.create(16, 0)

      if bit32.bxor(gb, 0b00001) == 1 then
        new_t1[1] = t1[1]
      end
      if bit32.bxor(gb, 0b00010) == 1 then
        new_t1[2], new_t1[3], new_t1[4], new_t1[5] = t1[2], t1[3], t1[4], t1[5]
      end
      if bit32.bxor(gb, 0b00100) == 1 then
        new_t1[6], new_t1[7], new_t1[8], new_t1[9], new_t1[10], new_t1[11] = t1[6], t1[7], t1[8], t1[9], t1[10], t1[11]
      end
      if bit32.bxor(gb, 0b01000) == 1 then
        new_t1[12], new_t1[13], new_t1[14], new_t1[15] = t1[12], t1[13], t1[14], t1[15]
      end
      if bit32.bxor(gb, 0b10000) == 1 then
        new_t1[16] = t1[16]
      end

      return setmetatable(new_t1, thor)
    end,
    -- similar to zerogate, except if any bit positions are on, then those grades will be negated
    -- negategate? minusgate
    minusgate = function(t1, gb)
      if gb == 0 then
        return t1
      end

      local new_t1 = table.clone(new_t1)

      local new_t1 = {-t1[1],
        -t1[2], -t1[3], -t1[4], -t1[5],
        -t1[6], -t1[7], -t1[8], -t1[9], -t1[10], -t1[11],
        -t1[12], -t1[13], -t1[14], -t1[15],
        -t1[16]}
      if bit32.bxor(gb, 0b00001) == 1 then
        new_t1[1] = t1[1]
      end
      if bit32.bxor(gb, 0b00010) == 1 then
        new_t1[2], new_t1[3], new_t1[4], new_t1[5] = t1[2], t1[3], t1[4], t1[5]
      end
      if bit32.bxor(gb, 0b00100) == 1 then
        new_t1[6], new_t1[7], new_t1[8], new_t1[9], new_t1[10], new_t1[11] = t1[6], t1[7], t1[8], t1[9], t1[10], t1[11]
      end
      if bit32.bxor(gb, 0b01000) == 1 then
        new_t1[12], new_t1[13], new_t1[14], new_t1[15] = t1[12], t1[13], t1[14], t1[15]
      end
      if bit32.bxor(gb, 0b10000) == 1 then
        new_t1[16] = t1[16]
      end

      return setmetatable(new_t1, thor)
    end,
    exp = function()

    end,
    pow = function()

    end,
    fuzzyEq = function(t1, t2, eps)
      eps = eps or 1e-9

    end,
	},
  -- regular unary minus
	__unm = function(t)
    local new_t = table.create(16, 0)
    new_t[1], new_t[2], new_t[3], new_t[4], new_t[5], new_t[6], new_t[7], new_t[8], new_t[9], new_t[10], new_t[11], new_t[12], new_t[13], new_t[14], new_t[15], new_t[16] = 
      -t[1], -t[2], -t[3], -t[4], -t[5], -t[6], -t[7], -t[8], -t[9], -t[10], -t[11], -t[12], -t[13], -t[14], -t[15], -t[16]

		return setmetatable(new_t, Thor)
	end,
  -- dual
  -- t = t_dual * I
  __len = function(t)
    local new_t = table.clone(t)
    new_t[1], new_t[2], new_t[4], new_t[6], new_t[8], new_t[10], new_t[12], new_t[14] = new_t[16], new_t[15], new_t[13], new_t[11], new_t[9], new_t[7], new_t[5], new_t[3]
    new_t[3], new_t[5], new_t[7], new_t[9], new_t[11], new_t[13], new_t[15], new_t[16] = -new_t[14], -new_t[12], -new_t[10], -new_t[8], -new_t[6], -new_t[4], -new_t[2], -new_t[1]
    return setmetatable(t, Thor)
  end,
  -- componentwise addition
	__add = function(t1, t2)
		if typeof(t1) == "number" then
			return thor(t1 + t2[1], t2[2], t2[3], t2[4], t2[5], t2[6], t2[7], t2[8], t2[9], t2[10], t2[11], t2[12], t2[13], t2[14], t2[15], t2[16])
		elseif typeof(t2) == "number" then
			return thor(t2 + t1[1], t1[2], t1[3], t1[4], t1[5], t1[6], t1[7], t1[8], t1[9], t1[10], t1[11], t1[12], t1[13], t1[14], t1[15], t1[16])
		end
		return thor(t1[1] + t2[1],
      t1[2] + t2[2], t1[3] + t2[3], t1[4] + t2[4], t1[5] + t2[5],
      t1[6] + t2[6], t1[7] + t2[7], t1[8] + t2[8], t1[9] + t2[9], t1[10] + t2[10], t1[11] + t2[11],
      t1[12] + t2[12], t1[13] + t2[13], t1[14] + t2[14], t1[15] + t2[15],
      t1[16] + t2[16])
	end,
  -- componentwise subtraction
	__sub = function(t1, t2)
		if typeof(t1) == "number" then
			return thor(t1 - t2[1], -t2[2], -t2[3], -t2[4], -t2[5], -t2[6], -t2[7], -t2[8], -t2[9], -t2[10], -t2[11], -t2[12], -t2[13], -t2[14], -t2[15], -t2[16])
		elseif typeof(t2) == "number" then
			return thor(t1[1] - t2, t1[2], t1[3], t1[4], t1[5], t1[6], t1[7], t1[8], t1[9], t1[10], t1[11], t1[12], t1[13], t1[14], t1[15], t1[16])
		end
		return thor(t1[1] - t2[1],
      t1[2] - t2[2], t1[3] - t2[3], t1[4] - t2[4], t1[5] - t2[5],
      t1[6] - t2[6], t1[7] - t2[7], t1[8] - t2[8], t1[9] - t2[9], t1[10] - t2[10], t1[11] - t2[11],
      t1[12] - t2[12], t1[13] - t2[13], t1[14] - t2[14], t1[15] - t2[15],
      t1[16] - t2[16])
	end,
  -- geometric product, also known as clifford product
	__mul = function(t1, t2)
		if typeof(t1) == "number" then
			return thor(t1 * t2[1],
        t1 * t2[2], t1 * t2[3], t1 * t2[4], t1 * t2[5],
        t1 * t2[6], t1 * t2[7], t1 * t2[8], t1 * t2[9], t1 * t2[10], t1 * t2[11],
        t1 * t2[12], t1 * t2[13], t1 * t2[14], t1 * t2[15],
        t1 * t2[16])
		elseif typeof(t2) == "number" then
			return thor(t2 * t1[1],
        t2 * t1[2], t2 * t1[3], t2 * t1[4], t2 * t1[5],
        t2 * t1[6], t2 * t1[7], t2 * t1[8], t2 * t1[9], t2 * t1[10], t2 * t1[11],
        t2 * t1[12], t2 * t1[13], t2 * t1[14], t2 * t1[15],
        t2 * t1[16])
		end

		local t = table.create(16, 0)

    local t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t1a, t1b, t1c, t1d, t1e, t1f = table.unpack(t1)
    local t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t2a, t2b, t2c, t2d, t2e, t2f = table.unpack(t2)

    t[1] = t10 * t20
       + t11 * t21 + t12 * t22 + t13 * t23 + t14 * t24
       + t15 * t25 + t16 * t26 + t17 * t27 + t18 * t28 + t19 * t29 + t1a * t2a
       + t1b * t2b + t1c * t2c + t1d * t2d + t1e * t2e
       + t1f * t2f
    t[2] = t11 * t20
       + t10 * t21 + t15 * t22 + t16 * t23 + t17 * t24
       + t12 * t25 + t13 * t26 + t14 * t27 + t1b * t28 + t1c * t29 + t1d * t2a
       + t18 * t2b + t19 * t2c + t1a * t2d + t1f * t2e
       + t1e * t2f
    t[3] = t12 * t20
       - t15 * t21 + t10 * t22 + t18 * t23 + t19 * t24
       - t11 * t25 - t1b * t26 - t1c * t27 + t13 * t28 + t14 * t29 + t1e * t2a
       - t16 * t2b - t17 * t2c - t1f * t2d + t1a * t2e
       - t1d * t2f
    t[4] = t13 * t20
       - t16 * t21 - t18 * t22 + t10 * t23 + t1a * t24
       + t1b * t25 - t11 * t26 - t1d * t27 - t12 * t28 - t1e * t29 + t14 * t2a
       + t15 * t2b + t1f * t2c - t17 * t2d - t19 * t2e
       + t1c * t2f
    t[5] = t14 * t20
       - t17 * t21 - t19 * t22 - t1a * t23 + t10 * t24
       + t1c * t25 + t1d * t26 - t11 * t27 + t1e * t28 - t12 * t29 - t13 * t2a
       - t1f * t2b + t15 * t2c + t16 * t2d + t18 * t2e
       - t1b * t2f
    t[6] = t15 * t20
       - t12 * t21 + t11 * t22 + t1b * t23 + t1c * t24
       - t10 * t25 - t18 * t26 - t19 * t27 + t16 * t28 + t17 * t29 + t1f * t2a
       - t13 * t2b - t14 * t2c - t1e * t2d + t1d * t2e
       - t1a * t2f
    t[7] = t16 * t20
       - t13 * t21 - t1b * t22 + t11 * t23 + t1d * t24
       + t18 * t25 - t10 * t26 - t1a * t27 - t15 * t28 - t1f * t29 + t17 * t2a
       + t12 * t2b + t1e * t2c - t14 * t2d - t1c * t2e
       + t19 * t2f
    t[8] = t17 * t20
       - t14 * t21 - t1c * t22 - t1d * t23 + t11 * t24
       + t19 * t25 + t1a * t26 - t10 * t27 + t1f * t28 - t15 * t29 - t16 * t2a
       - t1e * t2b + t12 * t2c + t13 * t2d + t1b * t2e
       - t18 * t2f
    t[9] = t18 * t20
       + t1b * t21 - t13 * t22 + t12 * t23 + t1e * t24
       - t16 * t25 + t15 * t26 + t1f * t27 - t10 * t28 - t1a * t29 + t19 * t2a
       - t11 * t2b - t1d * t2c + t1c * t2d - t14 * t2e
       - t17 * t2f
    t[10] = t19 * t20
       + t1c * t21 - t14 * t22 - t1e * t23 + t12 * t24
       - t17 * t25 - t1f * t26 + t15 * t27 + t1a * t28 - t10 * t29 - t18 * t2a
       + t1d * t2b - t11 * t2c - t1b * t2d + t13 * t2e
       + t16 * t2f
    t[11] = t1a * t20
       + t1d * t21 + t1e * t22 - t14 * t23 + t13 * t24
       + t1f * t25 - t17 * t26 + t16 * t27 - t19 * t28 + t18 * t29 - t10 * t2a
       - t1c * t2b + t1b * t2c - t11 * t2d - t12 * t2e
       - t15 * t2f
    t[12] = t1b * t20
       + t18 * t21 - t16 * t22 + t15 * t23 + t1f * t24
       - t13 * t25 + t12 * t26 + t1e * t27 - t11 * t28 - t1d * t29 + t1c * t2a
       - t10 * t2b - t1a * t2c + t19 * t2d - t17 * t2e
       - t14 * t2f
    t[13] = t1c * t20
       + t19 * t21 - t17 * t22 - t1f * t23 + t15 * t24
       - t14 * t25 - t1e * t26 + t12 * t27 + t1d * t28 - t11 * t29 - t1b * t2a
       + t1a * t2b - t10 * t2c - t18 * t2d + t16 * t2e
       + t13 * t2f
    t[14] = t1d * t20
       + t1a * t21 + t1f * t22 - t17 * t23 + t16 * t24
       + t1e * t25 - t14 * t26 + t13 * t27 - t1c * t28 + t1b * t29 - t11 * t2a
       - t19 * t2b + t18 * t2c - t10 * t2d - t15 * t2e
       - t12 * t2f
    t[15] = t1e * t20
       - t1f * t21 + t1a * t22 - t19 * t23 + t18 * t24
       - t1d * t25 + t1c * t26 - t1b * t27 - t14 * t28 + t13 * t29 - t12 * t2a
       + t17 * t2b - t16 * t2c + t15 * t2d - t10 * t2e
       + t11 * t2f
    t[16] = t1f * t20
       - t1e * t21 + t1d * t22 - t1c * t23 + t1b * t24
       - t1a * t25 + t19 * t26 - t18 * t27 - t17 * t28 + t16 * t29 - t15 * t2a
       + t14 * t2b - t13 * t2c + t12 * t2d - t11 * t2e
       + t10 * t2f

    setmetatable(t, Thor)

    return t
	end,
  -- product with inverse
	__div = function(t1, t2)
		if typeof(t1) == "number" then
			return t1 * t2.linverse
		elseif typeof(t2) == "number" then
			return t1 * (1 / t2)
		end
		
		return t1 * t2.linverse
	end,
  -- (exterior algebra) wedge product, also called join
	__pow = function(t1, t2)
    local t = table.create(16, 0)

    local t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t1a, t1b, t1c, t1d, t1e, t1f = table.unpack(t1)
    local t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t2a, t2b, t2c, t2d, t2e, t2f = table.unpack(t2)

    t[1] = t10 * t20
    t[2] = t10 * t21 + t11 * t20
    t[3] = t10 * t22 + t12 * t20
    t[4] = t10 * t23 + t13 * t20
    t[5] = t10 * t24 + t14 * t20
    t[6] = t11 * t22 - t12 * t21 + t15 * t20 + t10 * t25
    t[7] = t11 * t23 - t13 * t21 + t16 * t20 + t10 * t26
    t[8] = t11 * t24 - t14 * t21 + t17 * t20 + t10 * t27
    t[9] = t12 * t23 - t13 * t22 + t18 * t20 + t10 * t28
    t[10] = t12 * t24 - t14 * t22 + t19 * t20 + t10 * t29
    t[11] = t13 * t24 - t14 * t23 + t1a * t20 + t10 * t2a
    t[12] = t11 * t28 - t12 * t26 + t13 * t25 + t15 * t23 - t16 * t22 + t18 * t21 + t10 * t2b
    t[13] = t11 * t29 - t12 * t27 + t14 * t25 + t15 * t24 - t17 * t22 + t19 * t21 + t10 * t2c
    t[14] = t11 * t2a - t13 * t27 + t14 * t26 + t16 * t24 - t17 * t23 + t1a * t21 + t10 * t2d
    t[15] = t12 * t2a - t13 * t29 + t14 * t28 + t18 * t24 - t19 * t23 + t1a * t22 + t1e * t20 + t10 * t2e
    t[16] = t10 * t2f
        + t11 * t2e - t12 * t2d + t13 * t2c  - t14 * t2b
        + t15 * t2a - t16 * t29 + t17 * t28  + t18 * t27 - t19 * t26 + t1a * t25
        + t1b * t24 - t1c * t23 + t1d * t22  - t1e * t21
        + t1f * t20

		return setmetatable(t, Thor)
	end,
  -- regressive product, also called meet
	__concat = function(t1, t2)
    -- t1 .. t2 := #(#t1 ^ #t2)

    -- dual permutation:
    -- (f, e, -d, c, -b, a, -9, 8, -7, 6, -5, 4, -3, 2, -1, -0)

    local t = table.create(16, 0)

    local t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t1a, t1b, t1c, t1d, t1e, t1f = table.unpack(t1)
    local t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t2a, t2b, t2c, t2d, t2e, t2f = table.unpack(t2)
    --[=[
    t10, t11,  t12, t13,  t14, t15,  t16, t17,  t18, t19,  t1a, t1b,  t1c, t1d , t1e, t1f =
      t1f, t1e, -t1d, t1c, -t1b, t1a, -t19, t18, -t17, t16, -t15, t14, -t13, t12, -t11, -t10
    t20, t21,  t22, t23,  t24, t25 , t26, t27,  t28, t29,  t2a, t2b,  t2c, t2d,  t2e, t2f =
      t2f, t2e, -t2d, t2c, -t2b, t2a, -t29, t28, -t27, t26, -t25, t24, -t23, t22, -t21, -t20
      ]=]

    t[16] = t1f * t2f
    t[15] = t1f * t2e + t1e * t2f
    t[14] = t1f * t2d + t1d * t2f
    t[13] = t1f * t2c + t1c * t2f
    t[12] = t1f * t2b + t1b * t2f
    t[11] = t1e * -t2d + t1d * t2e + t1a * t2f + t1f * t2a
    t[10] = -t1e * t2c + t1c * t2e + t19 * t2f + t1f * t29
    t[9] = t1e * -t2b + t1b * t2e + t18 * t2f + t1f * t28
    t[8] = t1d * t2c + t1c * t2d + t17 * t2f + t1f * t27
    t[7] = t1d * t2b - t1b * t2d + t16 * t2f + t1f * t26
    t[6] = t1c * t2b - t1b * t2c + t15 * t2f + t1f * t25
    t[5] = t1e * -t27 - t1d * t29 + t1c * t2a + t1a * t2c - t19 * t2d - t17 * t2e + t1f * t24
    t[4] = t1e * t26 - t1d * t28 + t1b * t2a + t1a * t2b - t18 * t2d - t16 * t2e + t1f * t23
    t[3] = t1e * -t25 - t1c * t28 + t1b * t29 + t19 * t2b - t18 * t2c - t15 * t2e + t1f * t22
    t[2] = -t1d * t25 + t1c * t26 - t1b * t27 - t17 * t2b + t16 * t2c - t15 * t2d + t11 * t2f + t1f * t21
    t[1] = t1f * t20
        + t1e * t21 - t1d * t22 + t1c * t23 - t1b * t24
        + t1a * t25 - t19 * t26 + t18 * t27 + t17 * t28 - t16 * t29 + t15 * t2a
        + t14 * t2b - t13 * t2c + t12 * t2d - t11 * t2e
        + t10 * t2f

		return setmetatable(t, Thor)
    --return #(#t1 ^ #t2)
	end,
  -- sandwich product
	__mod = function(t1, t2)

	end,
	__eq = function(t1, t2)
		if typeof(t1) == "number" then
			for i = 2, 16 do
				if t2[i] ~= 0 then
					return false
				end
			end
			return t2[1] == t1
		end
		if typeof(t2) == "number" then
			for i = 2, 16 do
				if t1[i] ~= 0 then
					return false
				end
			end
			return t1[1] == t2
		end
		for i = 1, 16 do
			if t1[i] ~= t2[i] then
				return false
			end
		end
		return true
	end,
	__tostring = function(t)
		local strt = {}
		for i = 1, 16 do
			if t[i] ~= 0 then
				strt[#strt + 1] = string.format("%s%.3f%s", t[i] < 0 and "－" or (#strt > 0 and "＋" or ""), math.abs(t[i]), basis_tokens[i])
			end
		end
		if #strt == 0 then
			return "zero"
		end
		return table.concat(strt)
	end,
}
Thor.unary.Magnitude = Thor.unary.magnitude
Thor.unary.Unit = Thor.unary.unit
Thor.unary.Inverse = Thor.unary.inverse
Thor.unary.Reversion = Thor.unary.reversion
Thor.unary.Bar = Thor.unary.bar
Thor.unary.grade = Thor.unary.step
Thor.lib.join = Thor.lib.regressive
Thor.lib.meet = Thor.lib.wedge
Thor.__index = function(t, k)
	if Thor.unary[k] then
		return Thor.unary[k](t)
	end
	if Thor.lib[k] then
		return Thor.lib[k]
	end
	return Thor[k]
end
thor = function(a, w, x, y, z, wx, wy, wz, xy, xz, yz, wxy, wxz, wyz, xyz, wxyz)
	if getmetatable(a) == Thor then
		return a
	end
	local v = setmetatable({
    a or 0,
    w or 0, x or 0, y or 0, z or 0,
    wx or 0, wy or 0, wz or 0, xy or 0, xz or 0, yz or 0,
    wxy or 0, wxz or 0, wyz or 0, xyz or 0,
    wxyz or 0
  }, Thor)
	
	return v
end

Thor.lib.im = thor(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)
Thor.lib.iminv = thor(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1)
Thor.lib.zero = thor()
Thor.lib.one = thor(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

return thor
