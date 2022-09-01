function [x] = time_int(gl, x, T, dt, newton_opts)
  % compute number of time steps and adjust time step to find T - Nt * dt  =  0
  Nt = round(T / dt);
  theta = 0.5;
  dt = T / Nt;  
  N = length(x);
  n = N / 2;
  x_cur = x;
  x_old = x;
  fix_part = 0;
  compat_test(@(x)g(x, fix_part, gl, dt, theta), @(x)Jg(x, gl, dt, theta), N);
  figure(1)
  for i = 1 : Nt
    x = 2 * x_cur - x_old;
    x_old = x_cur;
    fix_part = x_cur + dt * (1 - theta) * gl.f(x_cur); 
    [x,  iter]  =  newton(x, @(x)g(x, fix_part, gl, dt, theta), @(x)Jg(x, gl, dt, theta), newton_opts);
    x_cur = x;
    iters(i) = iter;
    xc(i) = x(round(n / 2));
    yc(i) = x(n + round(n / 2));
    if rem(i, 2)  ==  0 %second argument should be adapted to change plotting interval
      subplot(3, 1, 1)
      plot(x(1 : n))
      hold on
      subplot(3, 1, 2)
      plot(x(n + 1 : 2 * n))
      hold on
    end
    if (i > 1)
        subplot(3, 1, 3)
        plot([xc(i - 1), xc(i)], [yc(i - 1), yc(i)])
        hold on
    end
    pause(0.05)
  end
  subplot(3, 1, 1),  hold off
  subplot(3, 1, 2),  hold off
  subplot(3, 1, 3),  hold off
  figure(2)
  plot(iters) 
end
function out =  g(x, fix_part, gl, dt, theta) 
  out = x - dt * theta * gl.f(x) - fix_part; 
end 
function out =  Jg(x, gl, dt, theta)
  N = length(x);
  out = speye(N, N) - dt * theta * gl.J(x); 
end