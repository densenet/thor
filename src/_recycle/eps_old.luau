local b = bit32.btest

local mul = (function()
  local ret = {}
  for i = 0, 255 do
    local value = 0
    if b(i, 0x2) then
      if b(i, 0x10) then
        value += 1
      end
    end
    if b(i, 0x4) then
      if b(i, 0x10) then
        value += 1
      end
      if b(i, 0x20) then
        value += 1
      end
    end
    if b(i, 0x8) then
      if b(i, 0x10) then
        value += 1
      end
      if b(i, 0x20) then
        value += 1
      end
      if b(i, 0x40) then
        value += 1
      end
    end
    if bit32.band(value, 1) == 1 then
      ret[i] = -1
    else
      ret[i] = 1
    end
  end

  -- to get index of E with two basis nibbles:
  -- bxor(lshift(right_basis, 4), left_basis)

  return ret
end)()

return {
  MUL = mul
}
