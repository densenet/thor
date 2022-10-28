local e = {}
setmetatable(e, e)
local D = 4

local b, l, x = bit32.btest, bit32.lshift, bit32.bxor

e.__index = function(t, cobw)
  --local c = 0
  local acc = false
  local negate = false
  for i = D - 1, 1, -1 do
    acc = if b(cobw, l(1, i)) then not acc else acc
    if b(cobw, l(1, (D - 1) + i)) and acc then
      negate = not negate
    end
    --c = if acc then x(l(1, i), c) else c
  end

  e[cobw] = if negate then -1 else 1

  return e[cobw]
end

return e
