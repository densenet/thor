-- this experiment regards the video on Error-Free Geometric Algebra
-- https://www.youtube.com/watch?v=f-0HthoQ188

-- from src\util\extended_gcd
-- EEA is useful for calculating y if we have h = xy^-1 mod p
local function eea(a: number, b: number): (number, number, number)
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

local function meea(p: number, h: number): (number, number, number)
  local x0, x1 = p, h
  local y0, y1 = 1, 0
  local z0, z1 = 0, 1
  local i = 1
  
  local N = math.floor(math.sqrt((p - 1) / 2))
  while x1 > N do
    local q = math.floor(x0 / x1)
    x0, x1 = x1, x0 - q*x1
    y0, y1 = y1, y0 - q*y1
    z0, z1 = z1, z0 - q*z1
    i += 1
  end
  
  if i % 2 == 1 then
    return x1, z1, -y1
  else
    return -x1, -z1, y1
  end
end

local P = 95520949

local __SANITY = false

if __SANITY then
  local h
  do
    local x, y = math.random(1, 10000), math.random(1, 10000)

    print(string.format("frac: %d / %d", x, y))

    local _, y_inv, _ = eea(y, P)
    y_inv = y_inv % P

    print(string.format("y^-1: %d", y_inv))

    h = (x * y_inv) % P
  end

  print(string.format("h = %d", h))

  print("recapturing (x,y,z) from (P,h) with meea")

  local x,y,z = meea(P,h)
  print(string.format("(x,y,z)=(%.2f,%.2f,%.2f)", x,y,z))
  
  print("recapturing h with diophantine equation")
  
  print(string.format("h = %d", (x+z*P)/y))
end

local __ADD = true

if __ADD then
  local x1, y1 = math.random(1, 1000), math.random(1, 1000)
  local x2, y2 = math.random(1, 1000), math.random(1, 1000)
  
  print(string.format("frac1: %d / %d", x1, y1))
  print(string.format("frac2: %d / %d", x2, y2))
  print(string.format("direct add: %d / %d", x1 * y2 + x2 * y1, y1 * y2))
  
  local _, y1_inv, _ = eea(y1, P)
  local _, y2_inv, _ = eea(y2, P)
  y1_inv %= P
  y2_inv %= P
  local h1 = (x1 * y1_inv) % P
  local h2 = (x2 * y2_inv) % P
  local hsum = (h1 + h2) % P
  
  print(string.format("(x,y,z)=(%.2f,%.2f,%.2f)", meea(P, hsum)))
end
