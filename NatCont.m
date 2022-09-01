function [x]=NatCont(gl,x,pars,NewtonOpts)
  n=length(x)/2;
  figure(1)
  for i=1:length(pars)
    gl.g1=pars(i)+j;
    [x, iter] = newton(x, @(x) gl.f(x), @(x) gl.J(x), NewtonOpts);
    lambdas=eigs(gl.J(x),4,'lr')
    if rem(i,2) == 0
      subplot(2,1,1)
      plot(x(1:n))
      hold on
      subplot(2,1,2)
      plot(x(n+1:2*n))
      hold on
    end
    pause(0.25)
    iters(i)=iter;
    lambdas1(i)=lambdas(1);
  end
  subplot(2,1,1), hold off
  subplot(2,1,2), hold off
  iters
  figure(2)
  subplot(2,1,1)
  plot(pars,iters) 
  subplot(2,1,2)
  plot(pars,real(lambdas1))
  hold on
  plot(pars,abs(imag(lambdas1)))
  hold off
end
