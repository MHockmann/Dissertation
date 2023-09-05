%close all
aufwand=1001;
n=1;
q=sqrt(2)/n;
x=linspace(-3*q,3*q,aufwand);
%% Phi initialisieren
phi= @(x) double(abs(x)<q/2).*cos(pi*x/q).^2;
phi_sampled=phi(x);
phi_gefaltet= @(x) double(abs(x)<q).*((q-abs(x))/4.*(1+0.5*cos(2*pi*x/q))+3/8*q/(2*pi)*sin(2*pi/q*abs(x)));
phi_gefaltet_sampled=phi_gefaltet(x);
second_derivative= @(x) double(abs(x)<q).*(-4*pi^2/q^2*phi_gefaltet(x)...
    +pi^2/q^2*(q/(2*pi)*sin(2*pi*abs(x)/q)+q-abs(x)));
second_derivative_sampled=second_derivative(x);
%% 2D Fensterfunktion durch Tensorproduktansatz
[X,Y]=meshgrid(x,x);
psi=((2*pi*n)^2-8*pi^2/q^2).*phi_gefaltet(X).*phi_gefaltet(Y)+...
    pi^2/q^2*double(abs(X)<q).*double(abs(Y)<q).*(((q/(2*pi)*sin(2*pi/q*abs(X))+q-abs(X)).*phi_gefaltet(Y)+...
    (q/(2*pi)*sin(2*pi/q*abs(Y))+q-abs(Y)).*phi_gefaltet(X)));
%psi=(2*pi*n)^2*phi_gefaltet(X).*phi_gefaltet(Y)+...
%    second_derivative(X).*phi_gefaltet(Y)+second_derivative(Y).*phi_gefaltet(X);
%% Higher powers in phi
k=4;
phi_k=@(x) double(abs(x)<q/2).*cos(pi*x/q).^k;
phi_sampled_k=phi_k(x);
phi_hat_k=fftshift(ifft(ifftshift(phi_sampled_k)));
phi_k_gefaltet=real(fftshift(fft(ifftshift((phi_hat_k.^2)))));
%% Modifiziertes psi über Tensorprodukt und inverse FT
v=(-aufwand+1)/2:(aufwand-1)/2;
v=v/(x(end)-x(1));
[VX,VY]=meshgrid(v,v);
[tensorprod1,tensorprod2]=meshgrid(phi_hat_k.^2,phi_hat_k.^2);
psi_k=4*pi^2*real(fftshift(fft2(ifftshift((n^2-VX.^2-VY.^2).*tensorprod1.*tensorprod2))));
%% 1D Faltung der Funktionen phi

figure(1)
plot(x,(1-(2*x/q).^2).*(abs(x)<q/2),'m:','Linewidth',3)
hold on
plot(x,cos(pi*x/q).*(abs(x)<q/2),'k.-')
hold on
plot(x,phi_sampled,'r-')
hold on
plot(x,phi_gefaltet_sampled,'b--')
legend('$\varphi_1(x)=1-\left(\frac{2x}{q}\right)^2$','$\varphi_2(x)=\cos\left(\frac{\pi x}{q}\right)$',...
    '$\varphi(x)=\cos^2\left(\frac{\pi x}{q}\right)$',...
    '$\varphi*\varphi(x)$','Interpreter','Latex','Fontsize',15)
xlabel('x','Interpreter','Latex')
% Test ob Faltung richtig ist
test1=integral(phi,-q,q)^2;
test2=integral(phi_gefaltet,-q,q);
error=abs(test1-test2);


%% 
h=figure(101);
%subplot(1,2,1)
scalefactor = @(x) x.^3;
descalefactor = @(x) sign(x).*abs(x).^(1/3);
Z3 = descalefactor(psi);

contourf(x,x,Z3,20,'LineStyle','none');
%axis equal
%axis tight
%imagesc(x,x,Z3)
%hsv = [0.6 1 0.7; .6 .8 0.8; .6 .6 0.9; .6 .4 1; .6 0.3 1; 0.6 0.2 1; 0.6 0.1 1; 0.6 0.05 1];
%rgb=hsv2rgb(hsv);
%colormap(rgb);

ma=max(Z3(:));
mi=min(Z3(:));
t0=linspace(mi,ma,1000)';
o=ones(1000,1);
t(t0<0)=1-(abs(t0(t0<0))/mi).^0.6;
t(t0>0)=1-(abs(t0(t0>0))/ma).^0.8;
r=o;
r(t0<0)=t(t0<0);
b=o;
b(t0>0)=t(t0>0);
c=min([o,r,r],[b,b,o]);
colormap(c);

hold on
rectangle('position',[-q -q 2*q 2*q],'Linestyle','--','Linewidth',1.5)
title('$\psi(x)$','Interpreter','Latex')
axis square
xlim([-2*q 2*q])
ylim([-2*q 2*q])
grid off
xlabel('$x_1$','Interpreter','Latex')
ylabel('$x_2$','Interpreter','Latex')
%colormap(jet)
%caxis([0,2.5])
a3 = colorbar;
ticks = get(a3,'YTick');
ticks = scalefactor(ticks);
set(a3,'YTickLabel',ticks);

% plot hatching region:
% zone(40:60,15:20) = nan; % force errorneous case
[~,h2] = contourf(x,x,Z3,[0.00001 0.00001]); % plots only the 0.00001 contour
set(h2,'linestyle','none','Tag','HatchingRegion');
hold off;                                % if you want to have more control
ax1 = gca;

% % Example 1: Default hatching
 hp = findobj(ax1,'Tag','HatchingRegion');
 hh = hatchfill2(hp,'cross','LineWidth',1,'Fill','off','HatchColor','blue');
 savefig(h,'psi')
 %title('Example 1: hatchfill2(hp,''HatchColor'',''w'',''FaceColor'',''none'')');

%% Fouriertrafo
h7=figure(102);
psi_hat=fftshift(real(ifft2(ifftshift(psi))))*aufwand;
%imagesc(v,v,descalefactor(psi_hat))
de_psi_hat=descalefactor(psi_hat);
contourf(v,v,de_psi_hat,20,'LineStyle','none');
hold on

%caxis([-0.1,0.4])

ma=max(de_psi_hat(:));
mi=min(de_psi_hat(:));
t0=linspace(mi,ma,1000)';
o=ones(1000,1);
t(t0<0)=1-(abs(t0(t0<0))/ma).^0.75;
t(t0>=0)=1-(abs(t0(t0>=0))/ma).^0.75;
r=o;
r(t0<0)=t(t0<0);
b=o;
b(t0>0)=t(t0>0);
c=min([o,r,r],[b,b,o]);
colormap(c);
% 
hold on
rectangle('Position',[-n -n 2*n 2*n],'Curvature',[1,1],'EdgeColor','k','Linestyle','--','Linewidth',1.5); %1,1, gives circle
xlim([-2*n 2*n])
ylim([-2*n 2*n])
xlabel('$v_1$','Interpreter','Latex')
ylabel('$v_2$','Interpreter','Latex')
axis square
title('$\hat{\psi}(v)$','Interpreter','Latex')
% 
a3 = colorbar;
ticks = get(a3,'YTick');
ticks = scalefactor(ticks);
set(a3,'YTickLabel',ticks);
% 
% plot hatching region:
% zone(40:60,15:20) = nan; % force errorneous case
[c2,h2] = contourf(v,v,de_psi_hat...
    ,[0.001 0.001]); % plots only the 0.001 contour
set(h2,'linestyle','none','Tag','HatchingRegion');
hold off;                                % if you want to have more control
ax1 = gca;
ax2 = copyobj(ax1,figure);

% % Example 1: Default hatching
 hp = findobj(ax1,'Tag','HatchingRegion');
 hh = hatchfill2(hp,'cross','LineWidth',1,'Fill','off','HatchColor','blue');
 
savefig(h7,'psi_hat')

% h1 = openfig('psi.fig','reuse'); % open figure
% ax1 = gca; % get handle to axes of figure
% h2 = openfig('psi_hat.fig','reuse');
% ax2 = gca;
% % test1.fig and test2.fig are the names of the figure files which you would % like to copy into multiple subplots
% h3 = figure; %create new figure
% s1 = subplot(1,2,1); %create and get handle to the subplot axes
% s2 = subplot(1,2,2);
% fig1 = get(ax1,'children'); %get handle to all the children in the figure
% fig2 = get(ax2,'children');
% copyobj(fig1,s1); %copy children to new parent axes i.e. the subplot axes
% copyobj(fig2,s2);

return

%%
n=1;
q=sqrt(2)/n;
radius=max(abs(X),abs(Y));
color={'k.','b.','g.','m.'};
K=[3,4,5,10];
y=linspace(0,q,1000);
for n=15:2:15
    figure(3+(n-15)/2)
    for i=1:2
        subplot(1,2,i)
        k=K(i);
        psi=((2*pi*n)^2-8*pi^2/q^2).*phi_gefaltet(X).*phi_gefaltet(Y)+...
            pi^2/q^2*double(abs(X)<q).*double(abs(Y)<q).*(((q/(2*pi)*sin(2*pi/q*abs(X))+q-abs(X)).*phi_gefaltet(Y)+...
            (q/(2*pi)*sin(2*pi/q*abs(Y))+q-abs(Y)).*phi_gefaltet(X)));
        psi_hat=fftshift(real(ifft2(ifftshift(psi))));
        plot(radius(:),(psi((aufwand-1)/2+1,(aufwand-1)/2+1)-psi(:))...
            /psi_hat((aufwand-1)/2+1,(aufwand-1)/2+1),'rd')
        hold on
        phi_k=@(x) double(abs(x)<q/2).*cos(pi*x/q).^k;
        phi_sampled_k=phi_k(x);
        phi_hat_k=fftshift(ifft(ifftshift(phi_sampled_k)));
        [tensorprod1,tensorprod2]=meshgrid(phi_hat_k.^2,phi_hat_k.^2);
        psi_k=4*pi^2*real(fftshift(fft2(ifftshift((n^2-VX.^2-VY.^2).*tensorprod1.*tensorprod2))));
        plot(radius(:),(psi_k((aufwand-1)/2+1,(aufwand-1)/2+1)-psi_k(:))...
            /(n^2*4*pi^2*tensorprod1((aufwand-1)/2+1,(aufwand-1)/2+1)...
            *tensorprod2((aufwand-1)/2+1,(aufwand-1)/2+1)),color{i})
        hold on
        plot(y(y>=0 & y<=q/2),5/4*pi^2/q^2*y(y>=0 & y<=q/2).^2/(4*pi^2*n^2*(q/2)^4),'r--');
        hold on
        plot(y(y>=q/2 & y<=q),(9/16*pi^2*y(y>=q/2 & y<=q))/q/(4*pi^2*n^2*(q/2)^4),'r--');
        hold on
        plot(y(y>=0 & y<=q),9/16*pi^2*y(y>=0 & y<=q).^2/q^2/(4*pi^2*n^2*(q/2)^4),'b--');
        xlabel('$\|x\|_{T^2}$','Interpreter','Latex')
        ylabel('$\frac{\psi(0)-\psi(x)}{\hat{\psi}(0)}$','Interpreter','Latex','Fontsize',12)
        xlim([0,q])
        ylim([-50,200])
        legend('$\frac{\psi(0)-\psi(x)}{\hat{\psi}(0)} (k=2)$',...
            ['$\frac{\psi_k(0)-\psi_k(x)}{\hat{\psi_k}(0)}$ for $k=$' num2str(k)],...
            'Lower bound for $k=2$','Interpreter','Latex','Fontsize',14,'Location','Southeast','color','none')
    end
    sgtitle(['Localising function around zero for ' '$d=2, q=\frac{\sqrt{2}}{15}\approx$ ' num2str(q,3) '$, n=$'...
        num2str(n)],'Interpreter','Latex')
end
    




