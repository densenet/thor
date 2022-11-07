local BASIS = {}

BASIS.IDX = {
	[0b0000] = 1, -- one
  [0b0001] = 2, -- w
	[0b0010] = 3, -- x
	[0b0100] = 4, -- y
	[0b1000] = 5, -- z
  [0b0011] = 6, -- wx
  [0b0101] = 7, -- wy
  [0b1001] = 8, -- wz
	[0b0110] = 9, -- xy
	[0b1010] = 10, -- xz
	[0b1100] = 11, -- yz
  [0b0111] = 12, -- wxy
  [0b1011] = 13, -- wxz
  [0b1101] = 14, -- wyz
  [0b1110] = 15, -- xyz
	[0b1111] = 16 -- wxyz
}

BASIS.LIST = {
  0b0000, -- scalar component
  0b0001, 0b0010, 0b0100, 0b1000, --  vector components
  0b0011, 0b0101, 0b1001, 0b0110, 0b1010, 0b1100, -- bi-vector components
  0b0111, 0b1011, 0b1101, 0b1110, -- tri-vector / pseudovector components
  0b1111 -- pseudoscalar component
}

BASIS.LABEL = {
	"_", "w", "x", "y", "z", "wx", "wy", "wz", "xy", "xz", "yz", "wxy", "wxz", "wyz", "xyz", "wxyz"
}

return BASIS
