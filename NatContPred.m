function [x]=NatContPred(gl,x,pars,NewtonOpts)
  n=length(x)/2;
  figure(1)
  xold=x;
  xcur=xold;
  for i=1:length(pars)
    gl.g1=pars(i)+j;
    x=2*xcur-xold;
    xold=xcur;
    [x, iter] = newton(x, @(x) gl.f(x), @(x) gl.J(x), NewtonOpts);
    xcur=x;
    iters(i)=iter;
    lambdas=eigs(gl.J(x),4,'lr')
    lambdas1(i)=lambdas(1);
    if rem(i,2) == 0 %second argument should be adapted to change plotting interval
      subplot(2,1,1)
      plot(x(1:n))
      hold on
      subplot(2,1,2)
      plot(x(n+1:2*n))
      hold on
      pause(0.25)
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
end