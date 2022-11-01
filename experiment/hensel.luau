-- this experiment regards the video on Error-Free Geometric Algebra

-- from src\util\extended_gcd
-- needed for MEEA: modified extended Euclidean algorithm
-- EEA is useful for calculating y if we have h = xy^-1 mod p
local function extended_gcd(a: number, b: number): (number, number, number)
  local m00, m01, m10, m11 = 1, 0, 0, 1
  local d, m
  while b ~= 0 do
    m = a % b
    d = (a - m) / b
    a, b = b, m
    m00, m01, m10, m11 =
      m01, m00 - d * m01,
      m11, m10 - d * m11
  end

  return a, m00, m10
end


