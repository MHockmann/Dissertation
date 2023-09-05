%% Setup
d=2;
nn=1233;
xx=linspace(0,1,nn);
NN=7*(nn-1)+1;
[XX,YY]=meshgrid(linspace(0,1,NN),linspace(0,1,NN));
X=[0.4, 0.25, 0.8, 0.1];
Y=[0.6, 0.15, 0.7, 0.5];


%% Compute Moments and SVD
n=20;
kk=0:n;
[kx, ky]=meshgrid(kk,kk);
lambda=[0.3 0.4 0.2 0.1];
kx=kx(:);
ky=ky(:);

A=exp(-2*pi*1i*(kx*X+ky*Y));
T=A*diag(lambda)*A';
[U,S,~]=svd(T);

%% 4 Points -> r=4
r=4;
U=U(:,1:r);
s=diag(S(1:r,1:r));
U_reshaped=zeros(n+1,n+1,r);
for j=1:r
    U_reshaped(:,:,j)=reshape(U(:,j),n+1,n+1);
end
p1=sum(abs(NN*NN*ifft2(U_reshaped,NN,NN)).^2,3)/(n+1)^2;
u3=abs(NN*NN*ifft2(U_reshaped(:,:,3),NN,NN)).^2/(n+1)^2;

%% Defining path through 2 points and plot
t=linspace(0.05,0.9,nn);
yt=0.2/0.7*(t-0.1)+0.5;

figure(1)
plot(X,Y,'rx')
hold on
c=gray;
c=flipud(gray);
imagesc(linspace(0,1,NN),linspace(0,1,NN),p1)
colormap(c)
hold on
plot(t,yt,'m-')
hold on
plot(X,Y,'rx')
axis on
axis([0,1,0,1])
axis equal
axis tight
ax=gca;
ax.LineWidth = 1;
xlabel('$x_1$','Interpreter','Latex')
ylabel('$x_2$','Interpreter','Latex')
%title('$p_{1,20}$ for discrete measure','Interpreter','Latex')

%% Interpolate along path
p1_path=interp2(XX,YY,p1,t,yt);
u3_path=interp2(XX,YY,u3,t,yt);
x2supps=max(abs(t'*ones(1,4)-ones(size(t'))*X),abs(yt'*ones(1,4)-ones(size(yt'))*Y));
dist2supp=min(max(abs(t'*ones(1,4)-ones(size(t'))*X),abs(yt'*ones(1,4)-ones(size(yt'))*Y)),[],2);
bound_local=1-(2*d-1)*3^(d-1)/(2^d*d^(2+d/2))*(n+1)^2*dist2supp(dist2supp<sqrt(2)/(n+1)).^2;
global_bound=1/(n+1)^2*max(lambda)/min(lambda)/3*sum(1./x2supps.^2,2);
global_constant=(1-3*pi^2/32/d^(d/2))*ones(length(t(dist2supp>sqrt(2)/(n+1))),1);
nearestxi=zeros(size(nn,1));
for i=1:nn
    nearestxi(i)=find(dist2supp(i)==x2supps(i,:));
end
%lower_bound=Fejer(t-X(nearestxi),n)/(n+1)...
   % .*Fejer(yt-Y(nearestxi),n)/(n+1);
lower_bound=max(Fejer(t'*ones(1,4)-ones(size(t'))*X,n).*...
    Fejer(yt'*ones(1,4)-ones(size(t'))*Y,n),[],2)/(n+1)^2;
%1-2*d^2*pi^2*(n)^2*dist2supp(dist2supp<sqrt(2)/(n+1)).^2;

%% Compute singular values numerically
tsmall=linspace(0.05,0.9,341);
ytsmall=0.2/0.7*(tsmall-0.1)+0.5;
sigma_min=zeros(1,length(tsmall));
for k=1:length(tsmall)
    Atilde=exp(-2*pi*1i*(kx*[X tsmall(k)]+ky*[Y ytsmall(k)]));
    sigma_min(k)=min(svd(Atilde));
end

%% Plot bounds
figure(2)%t(dist2supp>sqrt(2)/(n+1)),global_constant,'b.',...
plot(t,p1_path,'m-',t(dist2supp<sqrt(2)/(n+1)),bound_local,'g.',...
    t,global_bound,'r--',...
    tsmall,1-sigma_min.^2/(n+1)^d,'c-')
axis([0.05,0.9,0,2])
axis square
xlabel('$x_1$','Interpreter','Latex')
legend('Cross section of $p_{1,20}$', ...
    '$p_{1,n}(x)\leq 1-const\cdot n^2\cdot\min_t \|x-t\|^2$',...
    '$p_{1,n}(x)\leq const\cdot n^{-2} \sum_t |x-t|^{-2}$',...
    '$p_{1,n}(x)\leq 1-\sigma_{\min}(\tilde{A}_{n,x})^2/N$',...
    '$|u_3(x)|^2$',...
    'Interpreter','Latex','FontSize',12)
%title('Cross section and pointwise bounds')

%% Relation between singular vectores and Vandermonde matrix
B=(1/(n+1)*A)\U;
B_prod=B*B';
D=B_prod^(-1);
D_shifted=eye(4)-D;
%%
figure(3)
plot(t(2:end-1),diff(p1_path,2),t(2:end-1),diff(lower_bound,2),t,(1-p1_path)./(1-lower_bound'))
legend('difference','third order','second order')