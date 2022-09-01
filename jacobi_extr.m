function x = jacobi_extr(A, b, maxit, extr_intv)
  D = diag(diag(A))
  B = D - A;
  x = zeros(size(b));
  dif = 1;
  x_old = 0;
  for j = 1 : maxit
    x = D \ (b + B * x);
    x(1)
    dif_old = dif;
    dif = x - x_old;
    r = (dif' * dif_old) / (dif_old' * dif_old)
    if mod(j, extr_intv) == 0
      x = (x - r * x_old) / (1 - r);
    end
    x_old = x;
    norm(b - A * x)
  end
end
