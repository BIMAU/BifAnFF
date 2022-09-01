function main_alts_orbit(m)
  load OrbitJ
  NMp2 = size(Jhu, 1);
  K = Jhu;
  NM = NMp2 - 2;
  M = 20;
  K(1 : NM / M, NM - NM / M + 1 : NM) = 0;
  [P, eigvalsA] = eig(full(Jhu));
  [eigsort, ind] = sort(abs(diag(eigvalsA)), 'ascend');
  eigsort(1 : 10)
  figure
  plot(P(ind(1),  : ))
  %pause
  eigvals = eig(full(K \ Jhu))'
  figure
  plot(eigvals)
  figure
  plot(eigvals / 2 - 1)
  max(abs(eigvals / 2 - 1))
  switch 6
    case 1
      [V, D] = ort_sub_spac_iter(Jhu, m);
    case 2
      Y = jacobi_extr(Jhu, rand(NMp2, 1), 400, 420);
    case 3
      Y = gs_extr(Jhu, rand(NMp2, 1), NM / M, 400, 420);
    case 4
      rf = 0.5;
      [Y, rns] = stat_iter_split_extr(Jhu, rand(NMp2, 1), K, rf, 100, 120);
    case 5
      rf = 1;
      [Y, rns] = stat_iter_split_extr(Jhu, rand(NMp2, 1), K, rf, 100, 20);
    case 6
      [V, D, Y, rns] = gmres(Jhu, rand(NMp2, 1), K, 30);
  end  
  clear gl
  figure
  Hax=axes()
  semilogy(Hax,rns,'k', 'LineWidth',2);
  set(Hax,'FontSize',18)
  xlabel(Hax,"iteration")%,'fontsize', 14)
  ylabel(Hax,"relative residual")%,'fontsize', 14)
  %
  %semilogy(rns,'k','linewidth',2);
  %xlabel("iteration",'fontsize',15)
  %ylabel("relative residual",'fontsize',15)
end

function out =  Jh(Y, fxcur, Yp, zeta, theta, gl, dt, M) 
  NM = length(Y) - 2;
  N = NM / M;
  indx = 1 : N;
  gl.g1 = Y(NM + 1) + j;
  p = Y(NM + 2);

  Js = sparse(N, NM);
  f_g1s = zeros(NM, 1);
  fX = zeros(NM, 1);
  for i = 0 : M - 1
    Js( : , i * N + indx) = gl.J(Y(i * N + indx));
    f_g1s(i * N + indx) = gl.f_g1(Y(i * N + indx));
    fX(i * N + indx) = gl.f(Y(i * N + indx));
  end
  out = sparse(NM + 2, NM + 2);
  for i = 0 : M - 1
    out(i * N + indx, i * N + indx) = speye(N, N) - p * dt * theta * Js( : , i * N + indx); 
    out(i * N + indx, mod(i - 1, M) * N + indx) =  - (speye(N, N)  + p * dt * (1 - theta) * Js( : , mod(i - 1, M) * N + indx)) ;
    out(i * N + indx, NM + 1) =  - p * dt * (theta * f_g1s(i * N + indx)  +  (1 - theta) * f_g1s(mod(i - 1, M) * N + indx));
    out(i * N + indx, NM + 2) =  - dt * (theta * fX(i * N + indx)  +  (1 - theta) * fX(mod(i - 1, M) * N + indx));   
  end
  %out(NM + 1,  : ) = Yp';
  out(NM + 1, NM + 1) = Yp(NM + 1)';
  out(NM + 2, 1 : NM) = fxcur'; 
end
