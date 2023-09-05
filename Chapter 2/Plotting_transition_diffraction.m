%% Setup of parameters
n=30;
k=-n:n;
[kx,ky]=meshgrid(k,k);
I=find(kx.^2+ky.^2<=n^2);
kx=kx(I); ky=ky(I);
aufwand=15;
number=floor(n/1.4);
delta=linspace(0.75/n,1.23/n,aufwand);
s=zeros(aufwand,number);
k=0;
for j=1:number
    for i=1:aufwand
        a1=delta(i)*[1; 0]; a2=delta(i)*[0;1];%[0.5; sqrt(3)/2];
        [w1,w2]=meshgrid(0:j,0:j);
        w1=w1(:); w2=w2(:);
        T=w1*a1'+w2*a2';
        tx=T(:,1);
        ty=T(:,2);
        A0=exp(-2*pi*1i*(tx*kx(:)'+ty*ky(:)'));
        A1=(-2*pi*1i*ones((j+1)^2,1)*kx(:)').*A0;
        A2=(-2*pi*1i*ones((j+1)^2,1)*ky(:)').*A0;
        A0=A0'; A1=A1'; A2=A2';
        tmp=reshape([A0(:); A1(:); A2(:)],[],3*(j+1)^2);
        s(i,j)=svds(tmp,1,'smallest');
        k=k+1
    end
end
s=s';
%% Plotting
%%
condition=n./s;
%%
fig1=figure(1);
imagesc(delta*n,(2:number+1).^2,condition)
xlabel('sep\,$Y \cdot n$','Interpreter','Latex','Fontsize',14)
ylabel('Number of sources $|Y|$','Interpreter','Latex','Fontsize',14)
colorbar
set(gca,'ColorScale','log')
xticks([0.4 0.6 0.8 1 1.15 1.22 1.4 1.6])
caxis([5e-1,1e10]);
title('Visualising diffraction limit in $d=2$','Interpreter','Latex','Fontsize',16)
