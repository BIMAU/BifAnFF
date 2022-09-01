function rns = main_alts
  load xsteady
  n = length(x)
  gl =GLmodel(n / 2);
  x
  par
  gl.g1 = par + j;
  A = gl.J(x);
  full(A)
  switch 5
    case 1
      [V, D] = ort_sub_spac_iter(A, m);
    case 2
      y = jacobi_extr(A, rand(size(x)), 50, 55)
    case 3
      [y,rns] = gs_extr(A, rand(size(x)), 0, 500, -1)
    case 4
      [V, D] = arnoldi(A, 50);
    case 5
      K = tril(A);
      [V, D, x, rns] = gmres(A, rand(size(x)), K, 30);
  end
  clear gl
  Hax=axes()
    semilogy(Hax,rns,'k', 'LineWidth',2);
  set(Hax,'FontSize',18)
  %ax.XAxis.FontSize=20;
xlabel(Hax,"iteration")%,'fontsize', 14)
ylabel(Hax,"relative residual")%,'fontsize', 14)
end


