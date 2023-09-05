nn=1025;

n=29;

xx=0:1/nn:1;
xx(end)=[];
xx=xx-1/2;

r0=1/3;

k=-(nn-1)/2:(nn-1)/2;
[kx,ky]=meshgrid(k);

muh=2*pi*r0*besselj(0,-2*pi*r0*sqrt(kx.^2+ky.^2));

Fnh=1-abs(k)/(n+1);
Fnh(Fnh<0)=0;
Fnh=Fnh' * Fnh;
ph=Fnh.*muh;
p=real(fftshift(ifft2(ifftshift(ph))));

Dnh=ones(1,size(muh,1));
whos
Dnh(abs(k)>n)=0;
Dnh=Dnh' * Dnh;
ph0=muh.*Dnh;
p0=real(fftshift(ifft2(ifftshift(ph0))));

%c=3/4;
%mu(nn*x0)=0;

figure(1)
imagesc(xx,xx,p0)
colormap(1-gray(256))
%shading interp
axis equal
axis tight
set(gca,'FontSize',10);
%export_fig(gcf,'ExampleMeasure2d1','-pdf','-transparent')

figure(2)
imagesc(xx,xx,p)
colormap(1-gray(256))
%shading interp
%caxis(c)
axis equal
axis tight
set(gca,'FontSize',10);
%export_fig(gcf,'ExampleMeasure2d2','-pdf','-transparent')