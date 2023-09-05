function f_hat = exp_sum(Z,kX,kY)
%Exp_sum Computes the Fourier transform of a Dirac ensemble at grid points
%   Z: Matrix containing locations of Dirac measures and their weights
%   kX,kY: matrices of the same size representing a grid in Fourier space
[r,~]=size(Z);
tmp=zeros([size(kX) r]);
for l=1:r
    tmp(:,:,l)=Z(l,3).*exp(-2*pi*1i*kX*Z(l,1)-2*pi*1i*kY*Z(l,2));
end
f_hat=sum(tmp,3);
end

