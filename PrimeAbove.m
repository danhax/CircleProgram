
function y = PrimeAbove(xin)
x = xin;
%x = ceil(xin+1e-14);
while true
  y = primes(x);
  y = y(end);
  if y<=xin
    x = x+1;
  else
    break
  end
end
end

