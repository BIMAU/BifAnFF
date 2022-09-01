function [x,par]=PsArcCont(gl,x,par0,ds,Ns,NewtonOpts)
  n=length(x)/2;
  N=2*n;
  yold=[x;par0]
  ycur=yold;
  y=ycur;
  gl.g1=ycur(N+1)+j
  indx=1:N
  yp(indx,1)=zeros(N,1)
  zeta=1.0;
  yp(N+1,1)=ds;
 
  compat( @(y) g(y,ycur,yp,zeta,gl,ds), @(y) Jg(y,yp,zeta,gl,ds), N+1);
  figure(1)
  for i=1:Ns
   % gl.g1=ycur(N+1)+j;
    y=2*ycur-yold;
    yold=ycur;
    
    [y, iter] = newton(y, @(y) g(y,ycur,yp,zeta,gl,ds), @(y) Jg(y,yp,zeta,gl,ds), NewtonOpts);
    ycur=y;
    yp=ycur-yold;
    pars(i)=y(N+1);
    iters(i)=iter;
    lambdas=eigs(gl.J(y(indx)),4,'lr')
    lambdas1(i)=lambdas(1);
    if rem(i,2) == 0 %second argument should be adapted to change plotting interval
      subplot(2,1,1)
      plot(y(1:n))
      hold on
      subplot(2,1,2)
      plot(y(n+1:2*n))
      hold on
      pause(0.1)
    end
  end
  subplot(2,1,1), hold off
  subplot(2,1,2), hold off
  figure(2)
  subplot(2,1,1)
  plot(pars,iters) 
  subplot(2,1,2)
  plot(pars,real(lambdas1))
  hold on
  plot(pars,abs(imag(lambdas1)))
  hold off
  figure(3)
  plot(pars)
  pars
  x=y(indx)
  par=y(N+1)
end
function out= g(y,ycur,yp,zeta,gl,ds) 
  N=length(y)-1;
  indx=1:N;
  gl.g1=y(N+1)+j;
  out=[gl.f(y(indx)); zeta*yp(indx)'*(y(indx)-ycur(indx))+yp(N+1)*(y(N+1)-ycur(N+1))-ds^2];
end 
function out= Jg(y,yp,zeta,gl,ds) 
  N=length(y)-1;
  indx=1:N;
  gl.g1=y(N+1)+j;
  out=[gl.J(y(indx)), gl.f_g1(y(indx)); zeta*yp(indx)', yp(N+1)]; 
end