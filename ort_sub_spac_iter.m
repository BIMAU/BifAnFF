function  [V, D] = ort_sub_spac_iter(A, m)
n = size(A, 1);
maxit = 50;
V = rand(n, m) %+i*rand(n,m); 
[V, R] = qr(V, 0);
  for i = 1 : maxit
    W = A \ V;
    V' *W;
    [U, D] = eig(V' * W);
  %  W=W*U;
    [V, R] = qr(W, 0);
  end
  1 ./ diag(D)
  eigvals = sort(eig(A));
  eigvals(1:m)
end
  