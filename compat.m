function compat(f,J,n)
  epsi=1e-3
  x=rand(n,1); v=rand(n,1);

  if length(x)<100
    figure(6) 
    spy(J(x))
  end
  
  Jxv=J(x)*v; fx=f(x);
  figure(5)
  for i=1:10
    epsi=0.1*epsi
    dif=(f(x+epsi*v)-fx)/epsi-Jxv;
    norm(dif), plot(dif)
    hold on
  end
  hold off
end
