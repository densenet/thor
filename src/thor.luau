local modroot = script.Parent

local kgen = require(modroot.kgen)

export type Kernel = kgen.Kernel<string>

export type Thor<K = Kernel> = {
  unpack: (K) -> ...number,
}

local _Thor = {
  unpack = table.unpack,
}

local Thor = {}
Thor.__index = _Thor
Thor.__call = function(t, f)
  
end

return Thor
