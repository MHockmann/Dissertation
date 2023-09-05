%% Plotting Bessel functions
x=linspace(0,10,10000);
y0=besselj(0,x);
y1=besselj(1,x);
yhalf=besselj(0.5,x);
yminushalf=besselj(-0.5,x);

plot(x,y0,'r--',x,y1,'g-.',x,yhalf,'k',x,yminushalf,'m:')

legend('$J_0(x)$', '$J_1(x)$', '$J_{1/2}(x)$', '$J_{-1/2}(x)$','Interpreter','Latex','FontSize',18)
axis([-0.01 10.001 -0.5 2])
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

%% Beurling
x=linspace(-4,4,2000);
trunc=100;
k=(1:trunc)';
Beurling=@(x) (sin(pi*x)/pi).^2 .*(sum(1./(ones(trunc,1)*x-k*ones(size(x))).^2,1)-...
    sum(1./(ones(trunc,1)*x+k*ones(size(x))).^2,1)+2./x+1./x.^2);
B=Beurling(x);
plot(x,B)
ax = gca;
axis([-4 4 -1.35 1.35])
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
legend('$B(x)$','Interpreter','Latex','FontSize',18,'Location','Northwest')

%% Selberg
a=-1; b=1;
s_lower=-1/2*(Beurling(a-x)+Beurling(x-b));
s_upper=1/2*(Beurling(x-a)+Beurling(b-x));
s_Diederich=s_lower-(sin(pi*x)/pi).^2./(1-x.^2);
plot(x,s_lower,'r--',x,s_Diederich,'m',x,s_upper,'g-.',x,abs(x)<1,'k:')
legend('$s_{-}$','$\tilde{s}_{-}$','$s_{+}$','Interpreter','Latex','FontSize',18,'Location','Northwest')

%% Fourier transform of lower Selberg
s_lower_hat=(abs(x)<1).*(sin(2*pi*abs(x))/pi+1-abs(x));
s_Diederichs_hat=(abs(x)<1).*(sin(2*pi*abs(x))/(2*pi)+1-abs(x));
plot(x,s_lower_hat,'r--')
legend('$\hat{s}_{-}$','Interpreter','Latex','FontSize',18,'Location','Northwest')
ax = gca;
axis([-4 4 -0.2 1.2])
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
%% FT of Diederichs
plot(x,s_Diederichs_hat,'m')
legend('$\hat{\tilde{s}}_{-}$','Interpreter','Latex','FontSize',18,'Location','Northwest')
ax = gca;
axis([-4 4 -0.2 1.2])
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
