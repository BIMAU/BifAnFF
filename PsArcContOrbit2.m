function [X,par,p]=PsArcContOrbit(gl,X,par0,p,ds,dt,M,Ns,NewtonOpts)
  N=length(X)/M;
  NM=N*M;
  Yold=[X;par0;p];
  Ycur=Yold;
  Y=Ycur;
  gl.g2=Ycur(NM+1);
  indx=1:N;
  Yp=zeros(NM,1);
  zeta=0;
  theta=0.5;
  Yp(NM+1,1)=ds;
  Yp(NM+2,1)=0;
  fXcur=ones(NM,1);
  for i=0:M-1
    fXcur(i*N+indx)=gl.f(Y(i*N+indx));
  end
  compat( @(Y) h(Y,Ycur,fXcur,Yp,zeta,theta,gl,ds,dt,M), ...
        @(Y) Jh(Y,fXcur,Yp,zeta,theta,gl,dt,M),NM+2);
  figure(1)
  clf
  for i=1:Ns
   % Y=2*Ycur-Yold;
    Yold=Ycur;
    for l=0:M-1
       fXcur(l*N+indx)=gl.f(Ycur(l*N+indx));
    end
    [Y, iter] = newton(Y, @(Y) h(Y,Ycur,fXcur,Yp,zeta,theta,gl,ds,dt, M), ...
                 @(Y) Jh(Y,fXcur,Yp,zeta,theta,gl,dt,M),NewtonOpts);
    Ycur=Y;
    Yp=Ycur-Yold;
    pars(i)=Y(NM+1);
    ps(i)=Y(NM+2)
    iters(i)=iter;
    Jhu=Jh(Y,fXcur,Yp,zeta,theta,gl,dt,M);
    [A,B]=splitmatrix(Jhu(1:NM,1:NM),M);
    %lambdas=eigs(A,B,4,'lm') % not correct
    lambdas=eig(full(A),full(B));
    lambdas=lambdas.^M;
    [labs,indsort]=sort(abs(lambdas),'descend');
    lambdas(indsort(1:50))'
    labs(1:50)'
    lambdas1(i,1:4)=lambdas(indsort(1:M:4*M))
    nrmA(i)=norm(Y(1:N)+j*Y(N+1:2*N));
    if rem(i,1) == 0 %second argument should be adapted to change plotting interval
      plot([Y(round(N/4):N:NM);Y(round(N/4))],[Y(N/2+round(N/4):N:NM);Y(N/2+round(N/4))])
      hold on
      pause(0.1)
    end
  end
  hold off
  figure(2)
  subplot(2,1,1)
  plot(pars,iters)
  subplot(2,1,2)
  plot(pars,real(lambdas1))
  hold on
  plot(pars,abs(imag(lambdas1)))
  hold off
  figure(3)
  subplot(2,1,1)
  plot(pars)
  subplot(2,1,2)
  plot(pars,ps)
  figure(4)
  plot(pars,nrmA);
  X=Y(1:NM);
  par=Y(NM+1)
  p=Y(NM+2)
  save OrbitJ.mat Jhu
end
function out= h(Y,Ycur,fxcur,Yp,zeta,theta,gl,ds,dt,M)
  NM=length(Y)-2;
  N=NM/M;
  indx=1:N;
  gl.g2=Y(NM+1);
  p=Y(NM+2);
  fX=zeros(NM,1);
  for i=0:M-1
    fX(i*N+indx)=gl.f(Y(i*N+indx));
  end
  dumprev=Y(1:NM)+p*dt*(1-theta)*fX;
  dum=dumprev(NM-N+indx);
  dumprev(N+1:NM)=dumprev(1:NM-N);
  dumprev(indx)=dum;
  out=zeros(NM+2,1);
  out(1:NM)=Y(1:NM)-p*dt*theta*fX-dumprev;
  %out(NM+1)=Yp'*(Y-Ycur)-ds^2;
  out(NM+1)=Yp(NM+1)'*(Y(NM+1)-Ycur(NM+1))-ds^2;
  out(NM+2)=fxcur'*Y(1:NM);
end
function out= Jh(Y,fxcur,Yp,zeta,theta,gl,dt,M)
  NM=length(Y)-2;
  N=NM/M;
  indx=1:N;
  gl.g2=Y(NM+1);
  p=Y(NM+2);

  Js=sparse(N,NM);
  f_g1s=zeros(NM,1);
  fX=zeros(NM,1);
  for i=0:M-1
    Js(:,i*N+indx)=gl.J(Y(i*N+indx));
    f_g1s(i*N+indx)=gl.f_g2(Y(i*N+indx));
    fX(i*N+indx)=gl.f(Y(i*N+indx));
  end
  out=sparse(NM+2,NM+2);
  for i=0:M-1
    out(i*N+indx,i*N+indx)=speye(N,N)-p*dt*theta*Js(:,i*N+indx);
    out(i*N+indx,mod(i-1,M)*N+indx)=-(speye(N,N) +p*dt*(1-theta)*Js(:,mod(i-1,M)*N+indx)) ;
    out(i*N+indx,NM+1)=-p*dt*(theta*f_g1s(i*N+indx) + (1-theta)*f_g1s(mod(i-1,M)*N+indx));
    out(i*N+indx,NM+2)=-dt*(theta*fX(i*N+indx) + (1-theta)*fX(mod(i-1,M)*N+indx));
  end
  %out(NM+1,:)=Yp';
  out(NM+1,NM+1)=Yp(NM+1)';
  out(NM+2,1:NM)=fxcur';
end
function [A,B]=splitmatrix(J,M);
  NM=length(J);
  N=NM/M;
  indx=1:N;
  A=J;
  B=sparse(NM,NM);
    for i=0:M-1
      B(i*N+indx,i*N+indx)=J(i*N+indx,i*N+indx);
      A(i*N+indx,i*N+indx)=sparse(N,N);
    end
end

