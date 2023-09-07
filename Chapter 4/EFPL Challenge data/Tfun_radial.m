function z = Tfun_radial(x,~,mu,n)
%Tfun Computes the matrix vector product Tx=y
%   mu: single image with size LxK
%   x: vector from C^{N}\approx C^{\pi*n^2}
[L,K]=size(mu);
[vx,vy]=meshgrid(-n:n,-n:n);
X=reshape(x,2*n+1,2*n+1);
% if strcmp(tflag,'notrans') % Matrix T ist Hermitian
y=fft2(mu.*ifft2(X,L,K));
y=y(1:2*n+1,1:2*n+1);%.*(vx.^2+vy.^2<=n^2);
z=y.*(vx.^2+vy.^2<=n^2);
z=z(:);
end