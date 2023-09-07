%%% Calculate in uncorrelated stacksm
%% Anscombe transfor on the data
im=2*sqrt(double(tiffreadVolume("sequence-as-stack-MT0.N1.HD-2D-Exp.tif"))+3/8); 
%% Apodisation in order to see just one spike
[kx,ky]=meshgrid(1:64,1:64);
tmp=im(:,:,1);
tmp=tmp/sum(sum(tmp));
apo=zeros(size(tmp));
apo(sqrt((kx-31).^2+(ky-20).^2)<15)=exp(1/3)*exp(-1./(3*(1-...
    sqrt((kx(sqrt((kx-31).^2+(ky-20).^2)<15)-31).^2+...
    (ky(sqrt((kx-31).^2+(ky-20).^2)<15)-20).^2)/15)));
figure(1)
subplot(2,3,1)
imagesc(tmp.*apo)
colorbar
axis square
title('One spike in the first image')
%% FT
[kx,ky]=meshgrid(-32:31,-32:31);
FT=ifftshift(fft2(tmp))/numel(tmp);
OTF=FT;
save('parameterOTF.mat','OTF');
subplot(2,3,2)

imagesc(abs(FT))
axis square
colorbar
title('Fourier transform of this spike')
%% Plot radial values
r=sqrt(kx.^2+ky.^2);
subplot(2,3,3)
semilogy(r(:),abs(FT(:)),'r.')
title('Radial values')
%% Take quadratic mean over every radius
[r,I]=sort(r(:));
FT=FT(:);
FT=FT(I);
numbins=30;
t=linspace(0,45,numbins);
valuebin=zeros(numbins,1);
for j=1:numbins
    J=find(abs(r-t(j))<23/(numbins-1));
    valuebin(j)=sqrt(sum(abs(FT(J).^2))/numel(J));
end
%%
% 1. logarithm of Ansatz function in order to avoid focus on large absolute
% errors for small radius
g = @(A,x) log((abs(x)<A(2)).*(A(1).*exp(-A(3)*x.^2)) + A(4));
% ---Parameters---
A0 = [1e0,30,0.0005,0.5e-2];   % Inital (guess) parameters
% ---Fit---
% Define lower and upper bounds 
lb = [0,0,0,0];
ub = [realmax('double'),realmax('double'),1,5e0];
% Fit sample data
[A_opt,~,~,~,output] = lsqcurvefit(g,A0,t,log(valuebin'),lb,ub);
disp(output); % display summary of LSQ algorithm
subplot(2,3,4)
semilogy(t,valuebin,'b.')
hold on
semilogy(t,exp(g(A_opt,t)),'r--')
legend('Values of radially averaged data','Fitted curve')
title('Radial average and its gaussian fit') 
%%
% 1. logarithm of Ansatz function in order to avoid focus on large absolute
% errors for small radius
g_airy = @(A,x) log(A(1)*(abs(x/(2*A(3)))<1).*...
    (pi/4-asin(x/(2*A(3)))/2-x/(2*A(3))/2.*sqrt(1-(x/(2*A(3))).^2))+ A(4));
% ---Parameters---
A0_airy = [1e2,2,8,0.5*1e0];   % Inital (guess) parameters
InterpMethod='nearest'; % 'nearest','linear','spline','cubic'
% ---Fit---
% Define lower and upper bounds 
lb = [0,0,3,0];
ub = [realmax('double'),realmax('double'),14,1e1];
% Fit sample data
[A_airy,resnorm,res,flag,output] = lsqcurvefit(g_airy,A0_airy,t,log(valuebin'),lb,ub);
disp(output); % display summary of LSQ algorithm
subplot(2,3,5)
semilogy(t,valuebin,'b.')
hold on
semilogy(t,exp(g_airy(A_airy,t)),'r--')
legend('Values of radially averaged data','Fitted curve')
title('Radial average and its fit by theoretical OTF') 
%% Resulting PSF
[X,Y]=meshgrid(linspace(-0.5,0.5,65),linspace(-0.5,0.5,65));
% h=(besselj(1,4*pi*A_airy(3)*sqrt(X.^2+Y.^2))./sqrt(X.^2+Y.^2)).^2;
% h(33,33)=4*pi^2*A_airy(3)^2;
h=exp(-(X.^2+Y.^2)/A_opt(3)*pi^2);
subplot(2,3,6)
imagesc(h)
axis square
colorbar
title('PSF given by Gaussian fit')