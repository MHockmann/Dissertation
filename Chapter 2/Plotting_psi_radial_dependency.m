%% Finding Bessel zeros
j=[pi,0,pi];
index=1;
start=max(index+sign(index)*1.8*abs(index)^(1/3),1);
j(2) = fzero(@(z) besselj(1, z), start);
j(3) = fzero(@(z) besselj(1.5, z), start);

%% Evaluating psi
N=100;
v=linspace(0,12,N);
varphi_hat_samples=zeros(3,length(v));
y=linspace(0,1.5,N);
psi=zeros(3,N,2);
tau=[0 0.1];

for d=1:3
    varphi_hat_samples(d,:)=varphi_hat(v,d,j);
    for m=1:2
        for l=1:N
            if y(l)==0
                psi(d,l,m)=2*pi^(d/2)/gamma(d/2)*...
                    integral(@(omega) psi_hat(omega,v,tau(m),varphi_hat_samples,d).*omega.^(d-1),0,10);
            else
                psi(d,l,m)=2*pi*y(l)^(-d/2+1)*...
                    integral(@(omega) psi_hat(omega,v,tau(m),varphi_hat_samples,d)...
                    .*omega.^(d/2).*besselj(d/2-1,2*pi*omega*y(l)),0,10);
            end
        end
    end
end

%% Plotting psi
set(groot,'defaultAxesTickLabelInterpreter','latex'); 

plot(y,psi(1,:,1),'r--',y,psi(2,:,1),'m--',y,psi(3,:,1),'b--')

hold on 

plot(y,psi(1,:,2),'r:',y,psi(2,:,2),'m:',y,psi(3,:,2),'b:')

legend('$\psi_0$ for $d=1$','$\psi_0$ for $d=2$','$\psi_0$ for $d=3$',...
    '$\psi_\tau$ for $d=1$','$\psi_\tau$ for $d=2$','$\psi_\tau$ for $d=3$',...
    'Interpreter','Latex','FontSize',16)

ax = gca;
axis([0 1.5 0 70])
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

xticks([0 j(1)/pi j(2)/pi j(3)/pi])
xticklabels({'0','$\frac{j_{1/2,1}}{\pi}=1$','$\frac{j_{1,1}}{\pi}\approx 1.22$',...
    '$\frac{j_{3/2,1}}{\pi}\approx 1.43$'})

xlabel('$\|x\|_2$','Interpreter','Latex')

%% Functions
function y=varphi(r,d,j)
    y=(r<=j(d)/(2*pi)).*...
        (1-(j(d)./(2*pi*r)).^(d/2-1).*besselj(d/2-1,2*pi*r)/besselj(d/2-1,j(d)));
    y(r==0)=1-(j(d)/2).^(d/2-1)/gamma(d/2)/besselj(d/2-1,j(d));
end

function y=varphi_hat(v,d,j)
    n=length(v);
    y=zeros(1,n);
    for l=1:n
        if v(l)==0
            y(l)=2*pi^(d/2)/gamma(d/2)*integral(@(x) varphi(x,d,j).*x.^(d-1),0,2);
        else
            y(l)=2*pi*v(l)^(-d/2+1)*...
                integral(@(x) varphi(x,d,j).*x.^(d/2).*besselj(d/2-1,2*pi*x*v(l)),0,2);
        end
    end
end

function y=psi_hat(omega,v,tau,varphi_hat_samples,d)
    y=interp1(v,varphi_hat_samples(d,:),sqrt(1+tau)*omega);
    y=4*pi^2*(1+tau)*(1-omega.^2).*y.^2;
end
