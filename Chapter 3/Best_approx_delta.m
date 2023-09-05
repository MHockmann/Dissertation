N=1001;
x=linspace(-0.5,0.5,N);
n=10;
D_n=sin((2*n+1)*pi*x)./sin(pi*x);
D_n(floor(N/2))=2*n+1;
F_n=1/(n+1)*sin(pi*(n+1)*x).^2./sin(pi*x).^2;
F_n(floor(N/2)+1)=n+1;
m=floor(n/2);
J_n=3/((m+1)*(2*(m+1)^2+1))*sin(pi*(m+1)*x).^4./sin(pi*x).^4;
J_n(floor(N/2)+1)=3*(m+1)^3/(2*(m+1)^2+1);
%% Best approximation of \delta_0
k=transpose(0:(2*n+1));
B1_eval=0.5-k/(2*n+2);
B1_eval(1)=0;
pj_tilde_shifted=fft(B1_eval)/(2*n+2);
pj_tilde=[fliplr(pj_tilde_shifted(n+3:end));pj_tilde_shifted(1:n+1)];
pj=2*pi*1i*(-n:n)'.*pj_tilde; % Neglecting the coefficient p_0
p_tilde=2*real(ifft(pj_tilde(n+1:end),N)*N); % Interpolation polynomial
p_opt=1+2*fftshift(real(ifft(pj(n+1:end),N)))*N; % Best approximation to \delta_0

figure(1)

subplot(1,2,1)
X=linspace(0,1,100);
plot(x+0.5,p_tilde,'m--','LineWidth',2)
hold on
B1=0.5-X;
B1(1)=0;
B1(end)=0;
plot(X,B1,'k.','MarkerSize',4)
hold on
plot(k(2:end)/(2*n+2),B1_eval(2:end),'g*','MarkerSize',10)
title('Interpolation of Bernoulli spline $\mathcal{B}_1$','Interpreter','Latex')
legend(['Interpolation $\tilde{p}$ for $n=$' num2str(n)],'$\mathcal{B}_1$',['$2n+1=$' num2str(2*n+1) ' interpolation points'],'Interpreter','Latex')

subplot(1,2,2)
plot(x,D_n,'k:',x,F_n,'b-.',x,J_n,'r--',x,p_opt,'c-')
legend('$D_{10}$','$F_{10}$','$J_{10}$','$p^*$','Interpreter','Latex')
title('Different approximations to $\delta_0$','Interpreter','Latex')

% subplot(1,3,3)
% nn=linspace(-n,n,100)';
% coeffs_closed_form=pi*nn.*cot(nn*pi/(2*n+2))/(2*n+2);
% coeffs_closed_form(nn==0)=1;
% plot(-n:n,real(pj),'rx',nn,coeffs_closed_form,'b--')

%% Find behaviour of \|C_n\|_1
N_max=2000;
nn=1:N_max;
N_eval=2501;
C_n=zeros(N_max,N_eval);
coeffs=zeros(N_max,N_max);
for j=1:N_max
    coeffs(j,1:j)=pi*nn(1:j).*cot(nn(1:j)*pi/(2*j+2))/(2*j+2);
    C_n(j,:)=1+2*fftshift(real(ifft([0 coeffs(1:j)],N_eval)))*N_eval;
end
L1_C_n=sum(abs(C_n),2)/N_eval;
figure(2)
subplot(1,2,1)
plot(nn,L1_C_n./log(nn')) %-> looks like logarithmic growth

%% Find constant for convergence of Jackson kernel
x=linspace(-0.5,0.5,N_eval);
m=1:N_max;
J_int=zeros(N_max,1);
for l=1:N_max
    J_int(l)=2*integral(@(x) x.*Jackson(l,x),0,0.5);
end
subplot(1,2,2)
plot(2*m-2,J_int.*(2*m-1)')

