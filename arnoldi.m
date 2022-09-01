function  [V, D] = arnoldi(A, maxit)
n = size(A, 1);
v = rand(n, 1);
V = v / norm(v);
last = 1; 
H = [];
  for i = 1 : maxit
    v = A \ V(:, last);
    h = V' * v;
    v = v - V * h;
    hd = norm(v); 
    v = v / hd;
    H = [H, h; zeros(1, last - 1), hd];
    V = [V, v]; 
    last = last + 1;
    [W, D] = eig(H(1 : last - 1, :)); 
    %1 ./ diag(D(1 : min(last - 1, 10), 1 : min(last - 1, 10)))'
    %diag(V' * V)'
  endfor
  V = V(:, 1 : size(W, 1)) * W;
  D = diag(1 ./ diag(D));
  diag(D)
  eigvals = sort(eig(A));
  eigvals(1 : min(size(A, 1), maxit))
end
  