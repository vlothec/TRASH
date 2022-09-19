
movav = function(x, n = 5)
{
  filter(x, rep(1 / n, n), sides = 2)
}
