local modroot = script.Parent

local kgen = require(modroot.kgen)

export type Kernel = kgen.Kernel<string>

local Thor = {}
Thor.__call = function(t, sig_list, f)
   
end

return Thor
