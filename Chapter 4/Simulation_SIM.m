%% Setup
close all
noise=0.3; %max norm of noise
% OTF as in Ingerman
k0=10;
kx=-k0:k0; ky=-k0:k0;
[KX,KY]=meshgrid(kx,ky);
h_hat=otf(KX,KY,k0); 

% Exact Points in unit intervall away from the boundary for simplicity
r=10; %number of points
Z=0.1+0.8*rand(r,3); % first and second column for the points and third column ...
%for % weights
Z(:,3)=1;

% Easier Case
%Z=[[0.5 0.5] 1]; 

% Exponential sum
f_hat=exp_sum(Z,KX,KY);

%% Fourier transform of conventional image
D_hat=f_hat.*h_hat;

% estimate for conventional image
K=255; L=255; %image size

%extend D_hat and h_hat with zeros if necessary and invert FFT
h_hat_extended=h_hat;
if 2*k0+1<K
    D_hat=[zeros(2*k0+1,(K-1)/2-k0) D_hat zeros(2*k0+1,(K-1)/2-k0)];
    h_hat_extended=[zeros(2*k0+1,(K-1)/2-k0) h_hat_extended zeros(2*k0+1,(K-1)/2-k0)];
end
if 2*k0+1<L
    D_hat=[zeros((L-1)/2-k0,K); D_hat; zeros((L-1)/2-k0,K)];
    h_hat_extended=[zeros((L-1)/2-k0,K); h_hat_extended; zeros((L-1)/2-k0,K)];
end

D=ifft2(ifftshift(D_hat))+noise*randn(size(D_hat));
% Padding with zeros on the boundary!
padsize=0;
mask=ones(L-2*padsize,K-2*padsize);
mask=padarray(mask,padsize);
D=D.*mask;
%% For Plotting
x=linspace(0,1,K+1); y=linspace(0,1,L+1);
x=x(1:K);
y=y(1:L);
[X,Y]=meshgrid(x,y);
%% Simulating SIM images
D_SIM=cell(3,3); % 3 phases, 3 orientations
D_hat_SIM=cell(3,3);
phi=[0 2*pi/3 4*pi/3]; % phases
angle=[0 2*pi/6 2*2*pi/6];
p=k0*[cos(angle);
    sin(angle)]; % illumination directions
%p=floor(p);
c0=1; % modulation depth
for l=1:3
    for j=-1:1
        f_hat=exp_sum(Z,KX,KY)*exp(j*1i*phi(2)*0)+...
            c0/2*exp(j*1i*phi(2)*(1))*exp_sum(Z,KX-p(1,l),KY-p(2,l))+...
            c0/2*exp(j*1i*phi(2)*(-1))*exp_sum(Z,KX+p(1,l),KY+p(2,l));
        D_temp=f_hat.*h_hat;
        % extend with zeros if necessary 
        if 2*k0+1<K
            D_temp=[zeros(2*k0+1,(K-1)/2-k0) D_temp zeros(2*k0+1,(K-1)/2-k0)];
        end
        if 2*k0+1<L
            D_temp=[zeros((L-1)/2-k0,K); D_temp; zeros((L-1)/2-k0,K)];
        end
        % Invert FFT and add noise
        D_temp=ifftshift(D_temp);
        D_SIM{l,j+2}=ifft2(D_temp);
        D_SIM{l,j+2}=(D_SIM{l,j+2}+noise*max(max(D_SIM{l,j+2}))*rand(size(D))).*mask;
        D_hat_SIM{l,j+2}=fftshift(fft2(D_SIM{l,j+2}));
    end
end
%% Plotting single SIM image
figure(1)

subplot(2,2,1)
imagesc(x,y,real(D_SIM{1,1}))
title('angle=0, phase=0')
hold on
plot(Z(:,1),Z(:,2),'r*')
axis square
hold off
subplot(2,2,2)
imagesc(x,y,real(D_SIM{1,2}))
title('angle=0, phase=2\pi/3')
hold on
plot(Z(:,1),Z(:,2),'r*')
axis square
hold off
subplot(2,2,3)
imagesc(x,y,real(D_SIM{1,3}))
title('angle=0, phase=4\pi/3')
hold on
plot(Z(:,1),Z(:,2),'r*')
axis square
hold off
subplot(2,2,4)
imagesc(x,y,real(D_SIM{2,1}))
title('angle=\pi/3, phase=0')
hold on
plot(Z(:,1),Z(:,2),'r*')
axis square
hold off
sgtitle('Image using SIM-mode')
%% SIM Processing by Gustafsson
R_hat_translated=cell(3,3);
R_hat=cell(3,3);
h_hat_translated=cell(3,3);
product=cell(3,3); %F(exp(-2\pi\ii mp_lx) F^{-1} [h_hat*R_hat_translated*Apo_translated])

% Apodization
Kx=(-K+1)/2:(K-1)/2; Ky=(-L+1)/2:(L-1)/2;
[Kx,Ky]=meshgrid(Kx,Ky);
cutoff=1.9;
supp=find(sqrt(Kx.^2+Ky.^2)<cutoff*k0);
apo=zeros(L,K);
apo(supp)=exp(1./cutoff/k0)*exp(-1./(0.1*(cutoff*k0-sqrt(Kx(supp).^2+Ky(supp).^2))));
apo_translated=cell(3,3);

for m=-1:1
    for l=1:3 
        % equation 13 in Ingerman
        R_hat_translated{m+2,l}=1/3*(D_hat_SIM{l,1}*exp(-1i*m*phi(2)*(-1))...
            +D_hat_SIM{l,2}*exp(-1i*m*phi(2)*0)...
            +D_hat_SIM{l,3}*exp(-1i*m*phi(2)*1));
        
        % Translation in Fourierspace as a modulation in real space 
        R_hat{m+2,l}=fftshift(fft2(ifftshift((exp(-1i*2*pi/L*m*p(1,l)...
            *transpose((-L+1)/2:(L-1)/2))*...
            exp(-1i*2*pi/K*m*p(2,l)*((-K+1)/2:(K-1)/2)))...
            .*fftshift(ifft2(ifftshift(R_hat_translated{m+2,l})))))); 
        
        interpolant=scatteredInterpolant(Kx(:),Ky(:),...
            reshape(conj(h_hat_extended).*R_hat_translated{m+2,l},[],1),'natural','none');
        
        kxshifted=Kx+m*p(1,l); kyshifted=Ky+m*p(2,l);
        
        temp=interpolant(kxshifted,kyshifted);
        
        % Manche Bereiche sollten wir später vermeiden, dort lässt sich die
% Fourier-Transformierte nicht rekonstruieren
        temp(isnan(temp))=0; 
        
        R_hat{m+2,l}=temp;
        
        % Translation in Fourierspace as a modulation in real space
        h_hat_translated{m+2,l}=transpose(fftshift(fft2(ifftshift((exp(-1i*2*pi/L*m*p(1,l)...
            *transpose((-L+1)/2:(L-1)/2))*...
            exp(-1i*2*pi/K*m*p(2,l)*((-K+1)/2:(K-1)/2)))...
            .*fftshift(ifft2(ifftshift(h_hat_extended)))))));
        
        % Translation in Fourierspace as a modulation in real space
        apo_translated{m+2,l}=fftshift(fft2(ifftshift((exp(1i*2*pi/L*m*p(1,l)...
            *transpose((-L+1)/2:(L-1)/2))*...
            exp(1i*2*pi/K*m*p(2,l)*((-K+1)/2:(K-1)/2)))...
            .*fftshift(ifft2(ifftshift(apo))))));
        
        % Both translated together for Artefact-free reconstruction of the
        % nominator
        product{m+2,l}=fftshift(fft2(ifftshift((exp(-1i*2*pi/L*m*p(1,l)...
            *transpose((-L+1)/2:(L-1)/2))*...
            exp(-1i*2*pi/K*m*p(2,l)*((-K+1)/2:(K-1)/2)))...
            .*fftshift(ifft2(ifftshift(h_hat_extended.*R_hat_translated{m+2,l}.*apo_translated{m+2,l}))))));
    end
end
% Wiener filter
alpha=0.001;
b=[c0/2 1 c0/2];
temp1=zeros(size(D));
temp2=zeros(size(D));
temp_modified=zeros(size(D));
for m=-1:1
    for l=1:3
            temp1=temp1+conj(b(m+2))*R_hat{m+2,l};
            temp_modified=temp_modified+b(m+2)*product{m+2,l};
            temp2=temp2+abs(b(m+2)*h_hat_translated{m+2,l}).^2;
    end
end
f_hat_Wiener=temp1./(temp2+alpha);
f_hat_Wiener_modified=temp_modified./(temp2+alpha);

% Apodization bei greedy Implementierung
f_hat_Wiener=f_hat_Wiener.*apo;

% Inverse transform
f_Wiener=ifft2(ifftshift(f_hat_Wiener));
f_Wiener_modified=ifft2(ifftshift(f_hat_Wiener_modified));
%% Plotting result of SIM processing
figure(2)
subplot(1,2,1)
imagesc(x,y,real(f_Wiener))
hold on
plot(Z(:,1),Z(:,2),'r*')
axis square
hold off
title('Image with Wiener filter')
colorbar
subplot(1,2,2)
imagesc(real(f_hat_Wiener))
axis square
colorbar
hold off
title('f_hat_Wiener')
%% Testing
figure(3)
subplot(2,2,1)
test1=R_hat{3,2}./(h_hat_translated{3,2}+alpha);
imagesc(abs(test1));
title('temp1')
axis square
subplot(2,2,2)
test2=R_hat{3,2};
imagesc(real(test2));
title('R_hat_translated')
axis square
subplot(2,2,3)
imagesc(real(h_hat_translated{3,3}))
title('translated h')
axis square
subplot(2,2,4)
ta=f_hat_Wiener./(exp_sum(Z,Kx,Ky)+alpha);
imagesc(real(ta),[0 0.999])
colorbar
title('effective OTF')
axis square
%% Plotting the problem
figure(4)
subplot(2,2,1)
imagesc(110:145,110:145,real(h_hat_extended(110:145,110:145)))
title('$\hat{h}(k)$','Interpreter','Latex')
subplot(2,2,2)
imagesc(110:145,110:145,real(h_hat_translated{1,2}(110:145,110:145)))
title('$\hat{h}(k+mp_l)$','Interpreter','Latex')
subplot(2,2,3)
imagesc(real(R_hat_translated{3,2}))
title('$\hat{R}_{m,p_l}(k-mp_l)$','Interpreter','Latex')
subplot(2,2,4)
imagesc(real(R_hat{3,2}))
title('$\hat{R}_{m,p_l}(k)$','Interpreter','Latex')
sgtitle('Translation in Fourier space') 

%% Using a Pencil without SIM and Pencil-SIM in all three directions
%[Z_sim,~]=prony_SIM_2D(R_hat_translated,h_hat,b,k0,r);
%% Using Pencil-SIM in one direction
[Z_sim_new,Z_pencil]=pencil_SIM_2D(R_hat_translated,h_hat,b,k0,r);
%% Plotting
figure(5)
imagesc(x,y,real(D))
hold on
plot(Z(:,1),Z(:,2),'rs')
axis square
hold on
plot(Z_sim_new(:,1),Z_sim_new(:,2),'k+')
hold on
plot(Z_pencil(:,1),Z_pencil(:,2),'md')
legend('exact position',...
    'SIM in one direction','Matrix-Pencil without SIM','Location','Best')
title('2D-SIM','Interpreter','Latex')
%% Error control
%error_Pencil=find_distance(transpose(Z(:,1:2)),transpose(Z_pencil));
%error_sim=find_distance(transpose(Z(:,1:2)),transpose(Z_sim));
%error_sim_new=find_distance(transpose(Z(:,1:2)),transpose(Z_sim_new));

%% New SIM approach
% Separating components
f_ml=cell(3,3);
for l=1:3
    for m=-1:1
        f_ml{l,m+2}=1/3*(D_SIM{l,1}*exp(-1i*m*phi(2)*(-1))...
            +D_SIM{l,2}*exp(-1i*m*phi(2)*0)...
            +D_SIM{l,3}*exp(-1i*m*phi(2)*1));
    end
end

% Summing the components appropriately
f_total=zeros(size(D_SIM{1,1}));
h_total=f_total;
for l=1:3
    for m=-1:1
        f_total=f_total+exp(-2*pi*1i*m*(p(1,l)*X+p(2,l)*Y)).*f_ml{l,m+2};
        h_total=h_total+b(m+2)*transpose(h_hat_translated{m+2,l});
    end
end

% Recombination in Fourier domain
nominator_alternative=zeros(size(f_total));
denominator_alternative=zeros(size(f_total));
for l=1:3
    for m=-1:1
        f_alternative=ifftshift(fft2(fftshift(exp(-2*pi*1i*m*(p(1,l)*X+p(2,l)*Y)).*f_ml{l,m+2})));
        nominator_alternative=nominator_alternative+b(m+2)*h_hat_translated{m+2,l}.*f_alternative;
        denominator_alternative=denominator_alternative+abs(b(m+2)*h_hat_translated{m+2,l}).^2;
    end
    
end

f_hat_Wiener_alternative=nominator_alternative./(denominator_alternative+alpha).*apo;
f_Wiener_alternative=ifftshift(ifft2(fftshift(f_hat_Wiener_alternative)));

%%
figure(6)
imagesc(x,y,real(D_SIM{1,1}))
hold on
plot(Z(:,1),Z(:,2),'r*')
axis square
hold off
title('Noisy SIM Data','Interpreter','Latex','Fontsize',16)
colorbar
%%
imagesc(x,y,real(f_Wiener))
hold on
plot(Z(:,1),Z(:,2),'r*')
axis square
hold off
title('Result with Gustafsson method','Interpreter','Latex','Fontsize',16)
colorbar
%%
imagesc(x,y,real(ifft2(ifftshift(fftshift(fft2(f_total)).*h_total./(h_total.^2+alpha).*apo))))
hold on
plot(Z(:,1),Z(:,2),'r*')
colorbar
axis square
hold off
title('Algorithm 2', 'Interpreter','Latex','Fontsize',16)
%%
imagesc(x,y,real(f_Wiener_alternative))
hold on
plot(Z(:,1),Z(:,2),'r*')
colorbar
axis square
hold off
title('Algorithm 3','Interpreter','Latex','Fontsize',16)

