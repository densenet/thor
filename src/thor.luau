local modroot = script.Parent

local kgen = require(modroot.kgen)

export type Kernel = kgen.Kernel<string>

export type Thor<K = Kernel> = {
  unpack: (K) -> (number,
    number, number, number, number,
    number, number, number, number, number, number,
    number, number, number, number,
    number),
  _: (K) -> number,
  w: (K) -> number,
  x: (K) -> number,
  y: (K) -> number,
  z: (K) -> number,
  wx: (K) -> number,
  wy: (K) -> number,
  wz: (K) -> number,
  xy: (K) -> number,
  xz: (K) -> number,
  yz: (K) -> number, 
  wxy: (K) -> number,
  wxz: (K) -> number,
  wyz: (K) -> number,
  xyz: (K) -> number, 
  wxyz: (K) -> number,
  v: (K) -> (number, number, number, number),
  bv: (K) -> (number, number, number, number, number, number),
  pv: (K) -> (number, number, number, number),
  step: (K) -> number
}

local _Thor = {
  unpack = table.unpack,
  _ = function(t) return t[1] end,
  w = function(t) return t[2] end,
  x = function(t) return t[3] end,
  y = function(t) return t[4] end,
  z = function(t) return t[5] end,
  wx = function(t) return t[6] end,
  wy = function(t) return t[7] end,
  wz = function(t) return t[8] end,
  xy = function(t) return t[9] end,
  xz = function(t) return t[10] end,
  yz = function(t) return t[11] end,
  wxy = function(t) return t[12] end,
  wxz = function(t) return t[13] end,
  wyz = function(t) return t[14] end,
  xyz = function(t) return t[15] end,
  wxyz = function(t) return t[16] end,
  v = function(t) return t[2], t[3], t[4], t[5] end,
  bv = function(t) return t[6], t[7], t[8], t[9], t[10], t[11] end,
  pv = function(t) return t[12], t[13], t[14], t[15] end,
}

local _KCache = {}

local Thor = {}
Thor.__index = _Thor
Thor.__call = function(t)

end

return Thor
