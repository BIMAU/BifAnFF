classdef GLmodel < handle
% GLmodel: Ginzburg-Landau PDE, normalized to the interval
% [0,1]
% Complex formulation
% d/dt u = g1*u+g2*(d/dx)^2 u - g3*|u|^2 u
% u(0)=u(1)=0
% Real formulation

  properties
    % linear damping parameter
    g1=-0.1
    % diffusivity coefficient
    g2=1
    % nonlinear damping parameter
    % Assumed real
    g3=1.0
    % Dirichlet boundary condition right
    g4=1+1*i
    % number of internal grid points/ dofs
    N  = 200
    % grid increment size
    dx
    % second derivative, central discretization
    D2
  end

  methods
    function self = GLmodel(N)
      % constructors
      self.N  = N;
      self.dx = 1 / (self.N+1);
      self.computeMatrix();
    end

    function computeMatrix(self)
      % discretization second-order differential
      self.D2 = spdiags(kron([1, -2, 1],ones(self.N,1)), [-1:1], self.N, self.N)/self.dx^2;
    end

    function [out] = fcomplex(self, y)
      % RHS
      assert(numel(y) == self.N);
      out = (self.g1*speye(self.N,self.N) +  self.g2*self.D2)*y -self.g3*y.*y.*conj(y);
      out(self.N)=out(self.N)+ self.g2*self.g4/self.dx^2;
    end

    function [out] = f(self, y)
      assert(numel(y) == 2*self.N);
      dum=self.fcomplex(y(1:self.N) + i * y(self.N+1:2*self.N));
      out=[real(dum);imag(dum)];
    end

    function [out] = J(self, y)
      assert(numel(y) == 2*self.N);
      out=self.Jsub(y(1:self.N) + i * y(self.N+1:2*self.N));
    end

    function [out] = Jsub(self, y)
       % Jacobian
      assert(numel(y) == self.N);
      Jlin=self.g1*speye(self.N,self.N) +  self.g2*self.D2;
      ReJl=real(Jlin); ImJl=imag(Jlin);
      dum=spdiags(y.*conj(y),0,self.N,self.N);
      dum2=[diag(sparse(real(y)));diag(sparse(imag(y)))];
      out=[ReJl-self.g3*dum, -ImJl;ImJl ,ReJl-self.g3*dum]-2*self.g3*dum2*dum2';
    end

    function [out] = f_g1(self, y)
      % Derivative of f wrt g1
      assert(numel(y) == 2*self.N);
      out=y;
    end

    function [out] = f_g2(self, y)
      % Derivative of f wrt g1
      assert(numel(y) == 2*self.N);
      out=kron(speye(2,2),self.D2)*y;
    end

  end
end
