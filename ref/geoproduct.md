## geometric product

```luau
local x = bit32.bxor
-- e: epsilon - product term signs
local i_1
for i = 1, 16 do
  i_1 = bl[i]
  t_prod[i] = t1[i] * t2[1]
    + e[x(i_1, 0x10)] * (t1[bi[x(i_1, 0x1)]] * t2[2])  + e[x(i_1, 0x20)] * (t1[bi[x(i_1, 0x2)]] * t2[3])  + e[x(i_1, 0x40)] * (t1[bi[x(i_1, 0x4)]] * t2[4])  + e[x(i_1, 0x80)] * (t1[bi[x(i_1, 0x8)]] * t2[5])
    + e[x(i_1, 0x30)] * (t1[bi[x(i_1, 0x3)]] * t2[6])  + e[x(i_1, 0x50)] * (t1[bi[x(i_1, 0x5)]] * t2[7])  + e[x(i_1, 0x90)] * (t1[bi[x(i_1, 0x9)]] * t2[8])  + e[x(i_1, 0x60)] * (t1[bi[x(i_1, 0x6)]] * t2[9])  + e[x(i_1, 0xA0)] * (t1[bi[x(i_1, 0xA)]] * t2[10]) + e[x(i_1, 0xC0)] * (t1[bi[x(i_1, 0xC)]] * t2[11])
    + e[x(i_1, 0x70)] * (t1[bi[x(i_1, 0x7)]] * t2[12]) + e[x(i_1, 0xB0)] * (t1[bi[x(i_1, 0xB)]] * t2[13]) + e[x(i_1, 0xD0)] * (t1[bi[x(i_1, 0xD)]] * t2[14]) + e[x(i_1, 0xE0)] * (t1[bi[x(i_1, 0xE)]] * t2[15])
    + e[x(i_1, 0xF0)] * (t1[bi[x(i_1, 0xF)]] * t2[16])
end
```

## epsilon

```luau
local b = bit32.btest

local e = (function()
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
```
