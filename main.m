function main(n)
  %Make an instance of the problem
  gl=GLmodel(n);

  %Test compatibility of right-hand side and Jacobian matrix for small n
  if (n < 21)
    compat_test(@(x)gl.f(x),@(x)gl.J(x),2*n)
  end

  %Compute a nontrivial solution
  opts.tol=1e-10;
  opts.maxit = 20;
  x=zeros(2*n,1);
  [x, iter] = newton(x, @(x) gl.f(x), @(x) gl.J(x), opts);


  gl.g4=0.1+0.1*j; % sets right boundary condition
  %gl.g2=0.1
  gl.g1=0;

  switch 3
    case 1
      %Perform natural continuation in lin. term
      %Compute initial solution.
      pars=linspace(0,20,41)
      x=zeros(2*n,1);
      [x]=NatCont(gl,x,pars,opts)
    case 2
      %Same but now with a predictor added.
      pars=linspace(0,20,41)
      x=zeros(2*n,1);
      [x]=NatContPred(gl,x,pars,opts)
    case 3
      %Same but now with Pseudo-Arclength method
      x=zeros(2*n,1);
      par0=0;
      ds=0.2;
      Ns=80
      [x,par]=PsArcCont(gl,x,par0,ds,Ns,opts)
      save xsteady.mat x par
      [V,D]=OrtSubSpacIte(gl.J(x),3);
    case 4
      %Compute an orbit
      x=zeros(2*n,1);
      gl.g1=30+j;
      %Since this results in an imaginary part close to pm i
      % the period 2*pi
      T=40, dt =2*pi/20
      [x]=time_int(gl,x, T, dt, opts)
    case 5
      %Continuation of periodic orbit
      %construct an initial orbit from time integration
      x=zeros(2*n,1);
      gl.g1=30+j;
      %Since this results in an imaginary part close to pm i
      % the period 2*pi
      T=12; dt =2*pi/20;
      [x]=time_int(gl,x,T,dt, opts);
      %Now split up the orbit in M parts assuming period 2*pi
      M=20; % we could make other choices.
      p=2*pi;
      T=p/M;
      X=zeros(2*n*M,1);
      dt=T;
      for i=0:M-1
        [x]=time_int(gl,x,T,dt, opts); %only one time step
        X(i*2*n+[1:2*n])=x;
      end
      % perform continuation
      Ns=10;
      par0=30;
      ds=1.8;
      dt=1/M; %redefinition
      [X,par,p]=PsArcContOrbit(gl,X,par0,p,ds,dt,M,Ns,opts);
      save orbit.mat X par p
    %  par0=par
    %  gl.g1=par+j;
    %  ds=0.18
    %  %[X,par,p]=PsArcContOrbit(gl,X,par0,p,ds,dt,M,Ns,opts);
     case 6
      %Continuation of periodic orbit in parameter gamma_2
      %construct an initial orbit from time integration
      x=zeros(2*n,1);
      gl.g1=30+j;
      gl.g2=0.8;
      %Since this results in an imaginary part close to pm i
      % the period 2*pi
      T=12; dt =2*pi/20;
      [x]=time_int(gl,x,T,dt, opts);
      %Now split up the orbit in M parts assuming period 2*pi
      M=20; % we could make other choices.
      p=2*pi;
      T=p/M;
      X=zeros(2*n*M,1);
      dt=T;
      for i=0:M-1
        [x]=time_int(gl,x,T,dt, opts); %only one time step
        X(i*2*n+[1:2*n])=x;
      end
      % perform continuation
      Ns=8;
      par0=1;
      ds=0.2;
      dt=1/M; %redefinition
      [X,par,p]=PsArcContOrbit2(gl,X,par0,p,ds,dt,M,Ns,opts);
      save orbit.mat X par p
     % par0=par
      %gl.g1=par+j;
      %ds=0.18
      %[X,par,p]=PsArcContOrbit(gl,X,par0,p,ds,dt,M,Ns,opts);
  end
  clear gl
end


