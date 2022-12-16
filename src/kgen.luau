export type Code<T = string | any> = T
export type Kernel<T> = {
  __mul: Code<T>,
  __pow: Code<T>,
  __concat: Code<T>,
  __unm: Code<T>,
  __len: Code<T>,
  lc: Code<T>,
  rv: Code<T>,
  iv: Code<T>,
  cj: Code<T>,
}

export type Signature = {
  -- # of basis 1-vectors
  D: dim,
  -- handler dynamically decides & implements the cayley table
  handler: Handler,
  special: { Code <string> },
  inc_word: (state: any, w0: word) -> word,
}

type word = number
type dim = number

type HandlerState = { term_ct: number, term_buf: { Code<string> } }

-- specification for function generators which incrementally build cayley tables
type ProductGenerator = (bw_i: word, bw_j: word) -> ()
-- special case of above for automorphisms e.g. dual, reversion, conjugate
type MutateGenerator = (bw: word) -> ()

type Handler = (sig: Signature, op: string, state: HandlerState) -> (ProductGenerator & MutateGenerator, ...any)

local MAX_D = 16

local BINC, inc_word
if script then
  local util = script.Parent.util  
  
  BINC, inc_word = require(util.binc), require(util.inc_word)
else
  BINC, inc_word = require("src\\util\\binc"), require("src\\util\\inc_word")
end
BINC = BINC(MAX_D)

local function expand(state, sig: Signature, cb: (word, word) -> (), cb2: (word) -> ())
  local bw_A, bw_B = 0, 0
  
  local sig_inc_word = sig.inc_word
  
  repeat
    bw_B = 0
    repeat
      cb(bw_A, bw_B)
      
      bw_B = sig_inc_word(state, bw_B)
    until bw_B == 0
    cb2(bw_A)
   
    bw_A = inc_word(bw_A, sig.D)
  until bw_A == 0
end

local function replace(state, sig: Signature, cb: (word) -> ())
  local bw = 0
  
  local sig_inc_word = sig.inc_word
  
  repeat
    cb(bw)
    
    bw = sig_inc_word(state, bw)
  until bw == 0
end

local f = string.format
local function kgen(sig: Signature, op: string): Code<string>
  local state = {
    term_ct = 0,
  }
  
  local term_buf: { Code<string> } = {}
  state.term_buf = term_buf
  
  local gen: ProductGenerator, finalize = sig:handler(op, state)
  
  local k_buf: { Code<string> } = {}
  
  expand(state, sig, gen, function(bw_i)
    if state.term_ct > 0 then
      finalize(bw_i, k_buf)
      
      table.clear(term_buf)
      
      state.term_ct = 0
    end
  end)

  return table.concat(k_buf, "\n")
end

local function mgen(sig: Signature, op: string): Code<string>
  local state = {
    term_ct = 0
  } 
  
  local m_buf: { Code<string> } = {}
  state.term_buf = m_buf
  
  local gen: MutateGenerator, _ = sig:handler(op, state)
  
  replace(state, sig, gen)
  
  if state.term_ct == 0 then
    return ""
  end
  
  local buf1, buf2 = {}, {}
  
  for i = 1, state.term_ct do
    buf1[i],buf2[i] = table.unpack(m_buf[i]:split"=")
  end
  
  return table.concat(buf1, ",") .. "=" .. table.concat(buf2, ",")
end

local function _stub(): Code<string>
  return ""
end

local SIG_ALIAS = {
  -- complex numbers i^2=-1
  c   = ":,1",
  -- hyperbolic numbers i^2=1
  hyp = ":1",
  -- dual numbers i^2=0
  d   = ":,,1",
  -- quaternions
  quat  = ":,2",
  -- 2D plane-based geometric algebra
  pga2d = ":2,,1",
  -- 2D conformal geometric algebra
  cga2d = ":3,1",
  -- 3D plane-based geometric algebra
  pga3d = ":3,,1",
  -- 3D conformal geometric algebra
  cga3d = ":4,1",
  -- octonions
  octo  = ":,3",
  -- dual quaternions
  dquat = ":,2,1",
}

local _schemes = {}
local _dynsig_cache: { [string]: Signature } = {}
local function parse_dynsig(body: string): Signature
  if _dynsig_cache[body] then
    return _dynsig_cache[body]
  end
  
  local args: { string } = body:split(";")
  
  local scheme, opt = table.unpack((table.remove(args, 1) :: string):split(":"))

  if opt== nil then
    error("invalid signature: could not find a 'scheme:args' separator ':'")
  end
  
  for i = #args, 1, -1 do
    local k, v = table.unpack(args[i]:split"=")
    
    if v == nil then
      v = true
    end
    
    for _, k_sub in k:split"&" do
      args[k_sub] = v
    end
  end
  
  if _schemes[scheme] then
    local sig = _schemes[scheme](opt, args)
    
    assert(sig ~= nil, f("dynsig generator '%s' failed to produce a signature", body))
    assert(typeof(sig.D) == "number", f("dynsig generator '%s' produced invalid type '%s' for dimension field D", body, typeof(sig.D)))
    
    return sig
  end

  error(f("invalid signature: scheme '%s' is unknown", scheme))
end

local function defscheme(names: { string }, f: (opt: string, args: { string }) -> Signature): ()
  for i, v in names do
    _schemes[v] = f
  end
end

local EMPTY_STRING = ""
defscheme({EMPTY_STRING, "clifford"}, function(opt, args)
  return require("src\\clifford")(opt, args)
end)

defscheme({"dual"}, function(opt, args)
  error("'dual' scheme not implemented")
end)

local function genkernel(_signature: string): Kernel<string>
  local sig: Signature = parse_dynsig(SIG_ALIAS[_signature] or _signature)
  
  if sig.D > MAX_D then
    error(f("the signature '%s' carries too high dimension %d (must be >= %d)", _signature, sig.D, MAX_D))
  end

  local kernel = {
    __special = {}
  }
  
  do
    kernel.__mul = kgen(sig, "product") -- GEOMETRIC PRODUCT
    kernel.__pow = kgen(sig, "wedge") -- WEDGE PRODUCT / MEET
    kernel.__concat = kgen(sig, "regressive") -- REGRESSIVE PRODUCT / JOIN
    kernel.__len = mgen(sig, "dual") -- MV DUAL
    kernel.__unm = mgen(sig, "unary_minus")
    kernel.rv = mgen(sig, "reversion")
    kernel.iv = mgen(sig, "involution")
    kernel.cj = mgen(sig, "conjugate")
    kernel.lc = kgen(sig, "left_contraction") -- LEFT CONTRACTION
    kernel.fd = kgen(sig, "fat_dot") -- DOT PRODUCT / "FAT DOT"
    kernel.ip = kgen(sig, "inner_product") -- INNER PRODUCT
    kernel.sp = kgen(sig, "scalar_product") -- SCALAR PRODUCT
  end
  
  -- scheme specific extensions purely handled by the signature implementor
  -- may vary depending on parameters like dimension, or unspecified scheme-dependent parameters
  for k, v in sig.special do
    kernel.__special[k] = v
  end

  return kernel :: Kernel<string>
end


local function __template(sig)
  local k = genkernel(sig)
  
  local src = [[
return {
	__mul = function(t,a,b)
%s
	end,
	__pow = function(t,a,b)
%s
	end,
	__concat = function(t,a,b)
%s
	end,
	__len = function(t) %s end,
	__unm = function(t) %s end,
	rv = function(t) %s end,
	iv = function(t) %s end,
	cj = function(t) %s end,
	lc = function(t,a,b)
%s
	end,
	fd = function(t,a,b)
%s
	end,
	ip = function(t,a,b)
%s
	end,
	sp = function(t,a,b)
%s
	end
}]]
  
  return string.format(src, k.__mul, k.__pow, k.__concat,
    k.__len, k.__unm, k.rv, k.iv, k.cj,
    k.lc, k.fd, k.ip, k.sp)
end

print(__template ":3,,1;lgate&rgate=0,2,4")

--[[
return {
	__mul = function(t,a,b)
t[1]=a[1]*b[1]-a[9]*b[9]-a[10]*b[10]-a[11]*b[11]
t[6]=a[6]*b[1]+a[1]*b[6]+a[9]*b[7]+a[10]*b[8]-a[7]*b[9]-a[8]*b[10]-a[16]*b[11]-a[11]*b[16]
t[7]=a[7]*b[1]-a[9]*b[6]+a[1]*b[7]+a[11]*b[8]+a[6]*b[9]+a[16]*b[10]-a[8]*b[11]+a[10]*b[16]
t[8]=a[8]*b[1]-a[10]*b[6]-a[11]*b[7]+a[1]*b[8]-a[16]*b[9]+a[6]*b[10]+a[7]*b[11]-a[9]*b[16]
t[9]=a[9]*b[1]+a[1]*b[9]+a[11]*b[10]-a[10]*b[11]
t[10]=a[10]*b[1]-a[11]*b[9]+a[1]*b[10]+a[9]*b[11]
t[11]=a[11]*b[1]+a[10]*b[9]-a[9]*b[10]+a[1]*b[11]
t[16]=a[16]*b[1]+a[11]*b[6]-a[10]*b[7]+a[9]*b[8]+a[8]*b[9]-a[7]*b[10]+a[6]*b[11]+a[1]*b[16]
	end,
	__pow = function(t,a,b)
t[1]=a[1]*b[1]
t[6]=a[6]*b[1]+a[1]*b[6]
t[7]=a[7]*b[1]+a[1]*b[7]
t[8]=a[8]*b[1]+a[1]*b[8]
t[9]=a[9]*b[1]+a[1]*b[9]
t[10]=a[10]*b[1]+a[1]*b[10]
t[11]=a[11]*b[1]+a[1]*b[11]
t[16]=a[16]*b[1]+a[11]*b[6]-a[10]*b[7]+a[9]*b[8]+a[8]*b[9]-a[7]*b[10]+a[6]*b[11]+a[1]*b[16]
	end,
	__concat = function(t,a,b)
t[1]=a[16]*b[1]+a[11]*b[6]-a[10]*b[7]+a[9]*b[8]+a[8]*b[9]-a[7]*b[10]+a[6]*b[11]+a[1]*b[16]
t[6]=a[16]*b[6]+a[6]*b[16]
t[7]=a[16]*b[7]+a[7]*b[16]
t[8]=a[16]*b[8]+a[8]*b[16]
t[9]=a[16]*b[9]+a[9]*b[16]
t[10]=a[16]*b[10]+a[10]*b[16]
t[11]=a[16]*b[11]+a[11]*b[16]
t[16]=a[16]*b[16]
	end,
	__len = function(t) t[1],t[6],t[7],t[8],t[9],t[10],t[11],t[16]=t[16],t[11],-t[10],t[9],t[8],-t[7],t[6],t[1] end,
	__unm = function(t) t[1],t[6],t[7],t[8],t[9],t[10],t[11],t[16]=-t[1],-t[6],-t[7],-t[8],-t[9],-t[10],-t[11],-t[16] end,
	rv = function(t) t[1],t[6],t[7],t[8],t[9],t[10],t[11],t[16]=t[1],-t[6],-t[7],-t[8],-t[9],-t[10],-t[11],t[16] end,
	iv = function(t) t[1],t[6],t[7],t[8],t[9],t[10],t[11],t[16]=t[1],t[6],t[7],t[8],t[9],t[10],t[11],t[16] end,
	cj = function(t) t[1],t[6],t[7],t[8],t[9],t[10],t[11],t[16]=t[1],-t[6],-t[7],-t[8],-t[9],-t[10],-t[11],t[16] end,
	lc = function(t,a,b)
t[1]=a[1]*b[1]-a[9]*b[9]-a[10]*b[10]-a[11]*b[11]
t[6]=a[1]*b[6]-a[11]*b[16]
t[7]=a[1]*b[7]+a[10]*b[16]
t[8]=a[1]*b[8]-a[9]*b[16]
t[9]=a[1]*b[9]
t[10]=a[1]*b[10]
t[11]=a[1]*b[11]
t[16]=a[1]*b[16]
	end,
	fd = function(t,a,b)
t[1]=a[1]*b[1]-a[9]*b[9]-a[10]*b[10]-a[11]*b[11]
t[6]=a[6]*b[1]+a[1]*b[6]-a[16]*b[11]-a[11]*b[16]
t[7]=a[7]*b[1]+a[1]*b[7]+a[16]*b[10]+a[10]*b[16]
t[8]=a[8]*b[1]+a[1]*b[8]-a[16]*b[9]-a[9]*b[16]
t[9]=a[9]*b[1]+a[1]*b[9]
t[10]=a[10]*b[1]+a[1]*b[10]
t[11]=a[11]*b[1]+a[1]*b[11]
t[16]=a[16]*b[1]+a[1]*b[16]
	end,
	ip = function(t,a,b)
t[1]=-a[9]*b[9]-a[10]*b[10]-a[11]*b[11]
t[6]=-a[16]*b[11]-a[11]*b[16]
t[7]=a[16]*b[10]+a[10]*b[16]
t[8]=-a[16]*b[9]-a[9]*b[16]
	end,
	sp = function(t,a,b)
t[1]=a[1]*b[1]-a[9]*b[9]-a[10]*b[10]-a[11]*b[11]
	end
}

]]

return genkernel

