local bitct = (function()
  local ret = {}

  for i = 0, 255 do
    ret[i] = i%2 + (i%4 - i%2)/2  + (i%8 - i%4)/4 + (i%16 - i%8)/8 + (i%32 - i%16)/16 + (i%64 - i%32)/32 + (i%128 - i%64)/64 + (i%256 - i%128)/128
  end

  return ret
end)()

return bitct
