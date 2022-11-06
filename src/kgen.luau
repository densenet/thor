export type Code<T = string | any> = T
export type Kernel<T> = {
  __mul: Code<T>,
  __pow: Code<T>,
  __concat: Code<T>,
  __mod: Code<T>,
  __unm: Code<T>
}

export type Signature = {
  -- # of basis 1-vectors
  D: dim,
  -- handler dynamically decides & implements the cayley table
  handler: Handler,
  inc_word: (state: any, w0: word) -> word,
  SchemeData: any
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
  
  local gen: MutateGenerator = sig:handler(op, state)
  
  replace(state, sig, gen)
  
  return table.concat(m_buf, "\n")
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
    local k, v = table.unpack(string.split("="))
    
    if v == nil then
      v = true
    end
    
    args[k] = v
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
  error("not implemented")
end)

local function genkernel(_signature: string): Kernel<string>
  local sig: Signature = parse_dynsig(SIG_ALIAS[_signature] or _signature)
  
  if sig.D > MAX_D then
    error(f("the signature '%s' carries too high dimension %d (must be >= %d)", _signature, sig.D, MAX_D))
  end

  local kernel = {}
  
  do
    kernel.__mul = kgen(sig, "product") -- GEOMETRIC PRODUCT
    kernel.__pow = kgen(sig, "wedge") -- WEDGE PRODUCT / MEET
    kernel.__concat = kgen(sig, "regressive") -- REGRESSIVE PRODUCT / JOIN
    kernel.__len = _stub() -- MV DUAL
    kernel.__unm = _stub() -- UNARY MINUS
    kernel.__mod = _stub() -- SANDWICH PRODUCT (A % B == B * A * B.r)
  end

  return kernel :: Kernel<string>
end


do
  local k = genkernel ":,,6"
  
  print(k.__mul)
end

--[=[


]=]

--return genkernel

