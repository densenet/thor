-- this util should be inlined

local e = {}
setmetatable(e, e)

local b, l = bit32.btest, bit32.lshift

e.__index = function(t, cobw)
  local acc = false
  local negate = false
  local D = t.D
  for i = D - 1, 1, -1 do
    acc = if b(cobw, l(1, i)) then not acc else acc
    if b(cobw, l(1, (D - 1) + i)) and acc then
      negate = not negate
    end
  end

  e[cobw] = if negate then -1 else 1

  return e[cobw]
end

return e
