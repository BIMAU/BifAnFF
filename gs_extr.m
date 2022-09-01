function [x,nrs]=gs_extr(A, b, di, maxit, extrint)
  D = tril(A, di);
  B = D - A;
  x = zeros(size(b));
  dif = 1;
  xold = 0;
  for j = 1 : maxit
    x = D \ (b + B * x);
    x(1)
    difold = dif;
    dif = x - xold;
    r = (dif' * difold) / (difold' * difold)
    if mod(j, extrint) == 0
      'extrapolate'
      x = x + r * dif / (1 - r);
    end
    xold = x;
    nrs(j)=norm(b - A * x)
  end
end
