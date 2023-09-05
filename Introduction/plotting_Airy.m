sigma=1;
[X,Y]=meshgrid(linspace(-10,10,100),linspace(-10,10,100));
r=sqrt(X.^2+Y.^2);
Airy=@(r) 1/sigma^2/pi*(besselj(1,r/sigma)./(r/sigma)).^2;
% subplot(1,3,1)
imagesc(mat2gray(Airy(r).^(1/3.5)))
axis square
axis off
% subplot(1,3,2)
% surf(Airy(r))
% axis square
% axis off
%% Diffraction - Houston
delta=1.03*pi;
y=Airy(r)+Airy(sqrt((X-delta).^2+Y.^2));
% subplot(1,2,1)
% imshow(mat2gray(y.^(1/2.5)))
% axis square
% axis off
% % subplot(1,3,2)
% % surf(y)
% % axis square
% % axis off
% subplot(1,2,2)
x=linspace(-10,10,100);
plot(x/pi,Airy(x),'k')
hold on
plot(x/pi,Airy(x-delta),'b.-')
hold on
plot([1.03 1.03],[0,0.08],'k--')
hold on
plot([0 0],[0,0.08],'k--')
hold on
% an = annotation('doublearrow',[0.52 0.595],[0.51 0.51]);
% c = an.Color;
% an.Color = 'red';
% a = annotation('textarrow',[0.5 0.55],[0.7 0.51],'String','Full width at half maximum (FWHM)');
set(gca,'ytick',[])
axis square
xticks([0,delta/pi])
xticklabels({'0','1.03 \cdot n^{-1}'})
%% Diffraction - Sparrow
delta=0.94*pi;
y=Airy(r)+Airy(sqrt((X-delta).^2+Y.^2));
% subplot(1,2,1)
% imshow(mat2gray(y.^(1/2.5)))
% axis square
% axis off
% subplot(1,2,2)
x=linspace(-10,10,100);
plot(x/pi,Airy(x),'k')
hold on
plot(x/pi,Airy(x-delta),'b.-')
hold on
plot(x/pi,Airy(x)+Airy(x-delta),'r--')
hold on
plot([0 0],[0,0.08],'k--')
hold on
plot([0.94 0.94],[0,0.08],'b--')
set(gca,'ytick',[])
axis square
xticks([0,delta/pi])
xticklabels({'0','0.94 \cdot n^{-1}'})
%% Diffraction Rayleigh
delta=1.22*pi;
y=Airy(r)+Airy(sqrt((X-delta).^2+Y.^2));
% subplot(1,2,1)
% imshow(mat2gray(y.^(1/2.5)))
% axis square
% axis off
% subplot(1,2,2)
x=linspace(-10,10,100);
plot(x/pi,Airy(x),'k')
hold on
plot(x/pi,Airy(x-delta),'b.-')
hold on
plot([0 0],[0,0.08],'k--')
hold on
plot([delta/pi delta/pi],[0,0.08],'b--')
set(gca,'ytick',[])
axis square
xticks([0,delta/pi])
xticklabels({'0','1.22 \cdot n^{-1}'})
%% Diffraction perfect
delta=0.8*pi;
y=Airy(x)+Airy(x-delta);
x=linspace(-10,10,100);
plot(x/pi,max(y)/max(Airy(x))*Airy(x-delta/2),'k--')
hold on
x=linspace(-10,10,100);
plot(x/pi,y,'r')
%hold on
%plot(x/pi,Airy(x-delta),'k--')
legend('One Airy pattern','Superposition of two patterns')
title('Separation below Sparrows limit')
%% Composition of PSFs
R=3;
t=15*(rand(2,R)-0.5);
alpha=0.5+0.5*rand(1,R);
sigma=0.00001;
tmp=zeros(size(X));
for j=1:R
    r=sqrt((X-t(1,j)).^2+(Y-t(2,j)).^2);
    tmp=tmp+alpha(j)*Airy(r);
end
subplot(1,3,3)
imshow(mat2gray(tmp.^(1/2.5))')
axis square
axis off
