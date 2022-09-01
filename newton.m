function  [x, iter] = newton(x,f,df,opts)
  iter  = 0;
  nrmdx = Inf;
  while ((nrmdx > opts.tol) & (iter < opts.maxit))
    iter  = iter + 1;
    dx    = -df(x) \ f(x);
    x     = x + dx;
    nrmdx = norm(dx);
  end
end