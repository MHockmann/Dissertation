%%
A=tiffreadVolume('sequence-as-stack-MT0.N1.HD-2D-Exp.tif');
%%
original_SR=imread('sample.png');
%% Load fitted PSF
dimens=size(A);
load('parameterOTF.mat','OTF'); %Estimate of OTF from estimating_PSF.mat 
considered_set_length=2500;
% Anscombe transform
a=double(A(:,:,1:considered_set_length));
a=a./sum(sum(a,1),2);

%% Deconvolve data heuristically
[x,y]=meshgrid((-32:31));
r=sqrt(x.^2+y.^2);
gh=exp(-3.9e-1/64*r.^2);

figure(1)

subplot(2,3,1)
imagesc(a(:,:,1))
axis square
colorbar

subplot(2,3,2)
imagesc(gh);
axis square
colorbar

subplot(2,3,3)
imagesc(abs(fftshift(ifft2(gh))));
axis square
colorbar

% deconvolve data, restrict to low-frequencies
n=20;% makes sense because sqrt(20^2+20^2)\approx 28 is the expected cutoff

muh=ifftshift(fft2(a(:,:,1)))./gh;

muh=muh(33-n:32+n,33-n:32+n);

subplot(2,3,4)
imagesc(real(ifft2(fftshift(muh))));
axis square
colorbar

subplot(2,3,5)
imagesc(real(muh))
axis square
colorbar
%% ParFor Loop
tic
r_max=8;
trunc=0.75;
%updateWaitbar = waitbarParfor(considered_set_length, "Calculation on individual images...");
[vx,vy]=meshgrid(-n/2:n/2,-n/2:n/2);
%[set_index_row,set_index_col]=find(vx.^2+vy.^2<=n^2/4);
%N=length(set_index);
N=numel(vx);
N_eff=numel(find(r<=n/2));

%%
nn=2048;
q=zeros(nn,nn,considered_set_length);
parfor j=1:considered_set_length
    muh=ifftshift(fft2(a(:,:,j)))./gh;
    muh=muh(33-n:32+n,33-n:32+n);
    mu=ifft2(fftshift(muh));

    [P,S,~]=svds(@(x,tflag) Tfun_radial(x,tflag,mu,n/2), [N N],r_max); % sparse SVD
    s=diag(S);
    R=find(s>trunc*s(1),1,'last');
    u1=zeros(2*n+1,2*n+1);
    ujx=cell(R,1);
    epsilon=0;
    for l=1:R
        u1=reshape(P(:,l),n+1,n+1);
        ujx=abs(ifft2(u1,nn,nn)).^2*nn^4;
        q(:,:,j)=q(:,:,j)+(1-epsilon/s(l))*ujx;
    end
    
    %updateWaitbar();
end

%%
epsilon=1;
q_eps=epsilon./(N_eff-q);
t_tries=toc;
%% Putting images together
figure(2)
rmap = [linspace(0,1,30) ones(1,70)];
gmap= [linspace(0,0.6^2,30) sqrt(linspace(0.6^2,1,70))];
bmap= [linspace(0,0.1250^2,30) linspace(sqrt(0.125),0,70).^2];
map=zeros(100,3);
map(:,1)=rmap;
map(:,2)=gmap;
map(:,3)=bmap;
imagesc(max(q_eps,[],3))
colormap(map)
axis square
colorbar
axis off
%% Single image
figure(3)
imagesc(a(:,:,1))
colormap(map)
axis square
colorbar
axis off
figure(4)
imagesc(q_eps(:,:,1))
colormap(map)
axis square
colorbar
axis off