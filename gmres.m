function  [V, D, x, rns] = gmres(A, b, K, maxit)
  n = size(A,1);
  v = K \ b;
  bn = norm(v);
  V = v / bn;
  last = 1; 
  H = [];
  Q = [];
  rn = 1;
  rns = rn;
  R = [];
  for i = 1:maxit
    v = K \ (A * V(:, last));
    h = V' * v;
    v = v - V * h;
    hd = norm(v); 
    v = v / hd;
    H = [H, h; zeros(1, last - 1), hd];
    V = [V, v]; last = last + 1;
    [W, D] = eig(H(1:last - 1, :)); 
    1 ./ diag(D(1:min(last - 1, 10), 1:min(last - 1, 10)))'
    r = Q' * h; 
    if i==1, 
      q = h; 
    else 
      q = h - Q * r; 
    end, 
    rd = norm([q;hd]); 
    q = q / rd;
    qd = hd / rd;
    if i==1, 
      Q = [q; qd]; 
    else 
      Q = [Q, q; zeros(1, last - 2), qd]; 
    end
    R = [R, r; zeros(1, last - 2), rd]; 
    rn = [rn - q(1) * q; - q(1) * qd ];
    yhat = bn * (R \ Q(1, :)');
    rns(last) = norm(rn);
    bn * norm(rn)    
  end
  %V = V(:, last - 1) * W;
  D = diag(1 ./ diag(D));
  yhat = bn * (R \ Q(1, :)');
  x = V(:, 1:last - 1) * yhat;
  norm(K \ (b - A * x))
  norm(b - A * x)
end
  