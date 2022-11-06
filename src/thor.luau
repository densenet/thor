local modroot = script.Parent

local kgen = require(modroot.kgen)

export type Kernel = kgen.Kernel<string>

local Thor = {}
Thor.__call = function(t, f)
   
end

return Thor
