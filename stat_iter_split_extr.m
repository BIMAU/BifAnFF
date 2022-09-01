function [x,rns]=stat_iter_split_extr(A, b, K, w, maxit, extr_intv)
  B = K - A;
  x = zeros(size(b));
  dif = 1;
  xold = 0;
  nb = norm(K \ b);
  for j = 1 : maxit
    x = (1 - w) * x + K \ (b + B * x) * w;
    x(1)
    difold = dif;
    dif = x - xold;
    r = (dif' * difold) / (difold' * difold)
    if mod(j, extr_intv) == 0
      'extrapolate'
      x = x + r * dif / (1 - r);
    end
    xold = x;
    rns(j) = norm(K \ (b - A * x)) / nb;
    rns(j)
  end
end
