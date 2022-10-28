return function(MAX_N: number)
  local t = table.create(MAX_N + 1, 0)
  for i = 1, MAX_N + 1 do
    t[i] = table.create(MAX_N + 1, 0)
  end

  t[1][1] = 1
  for n = 2, MAX_N + 1 do
    t[n][1], t[n][n] = 1, 1
    for k = 2, n do
      t[n][k] = t[n - 1][k - 1] + t[n - 1][k]
    end
  end

  return t
end
