function h = otf(kx,ky,k_0)
%OTF Computes the optical transfer function as in Ingerman
%   kx,ky: matrices of same size representing grid in Fourier space
%   k0: Cut-off frequency or support of the OTF
h=zeros(size(kx));
reg=find(sqrt(kx.^2+ky.^2)<k_0); %support of OTF
b=acos(sqrt(kx(reg).^2+ky(reg).^2)/k_0);
h(reg)=(2*b-sin(2*b))/pi;
end

