nn=2048;

n=19;

xx=0:1/nn:1;
xx(end)=[];

mu=zeros(nn,1);
mu(2*nn/8:5*nn/8-1)=8/9;

x0=1/8;
mu(nn*x0)=nn/3;

xx1=(6*nn/8:8*nn/8-1)/nn;
mu(6*nn/8:8*nn/8-1)=(1./sqrt(abs(xx1-7/8))-sqrt(8))*(sqrt(2)/3);
mu(7/8*nn)=mu(7/8*nn-1);

muh=fft(mu);
k=0:nn-1;
k=min(k',nn-k');
Fnh=1-k/(n+1);
Fnh(Fnh<0)=0;
ph=muh.*Fnh;
p=ifft(ph);

Dnh=ones(size(muh));
Dnh(k>n)=0;
ph0=muh.*Dnh;
p0=ifft(ph0);

c=1/3;
mu(nn*x0)=0;

figure(1)
plot(xx,mu,'k',[x0,x0],[0,c*n],'k',[x0-0.005,x0],[0.95*c*n,c*n],'k',[x0+0.005,x0],[0.95*c*n,c*n],'k',xx,p0,'k:','Linewidth',2);

axis([0,1,-2.1,15])
set(gca,'FontSize',10);
%export_fig(gcf,'ExampleMeasure1','-pdf','-transparent')

figure(2)
plot(xx,mu,'k',[x0,x0],[0,c*n],'k',[x0-0.005,x0],[0.95*c*n,c*n],'k',[x0+0.005,x0],[0.95*c*n,c*n],'k',xx,p,'b-.','Linewidth',2);

axis([0,1,-2.1,15])
set(gca,'FontSize',10);
%export_fig(gcf,'ExampleMeasure2','-pdf','-transparent')