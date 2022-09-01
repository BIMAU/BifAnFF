load orbit
X2D=reshape(X,20,20);
Xup=X2D(1:10,:),Xlow=X2D(11:20,:);
x = linspace(0,1,10);
y = linspace(0,1,21);
[X,Y] = meshgrid(x,y);
figure(10)
Z=sqrt(Xup.*Xup+Xlow.*Xlow)'
  Z=[Z;Z(1,:)];
colormap( flipud(gray(256)) );
contourf(X,Y,Z)
colorbar
xlabel('x','FontSize',20)
ylabel('t','FontSize',20)
p
par
