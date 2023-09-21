close all
%% Displaying Christoffel function in univariate example
x=[0.16;0.33;0.75];
c=[1/6;1/3;1/2];

N_max=6;
for n=4:N_max
 k=0:n;
 A=exp(2*pi*1i*x*k);
 T=A'*diag(c)*A;
 Trho=0.05*randn(size(T))*max(abs(T(:)));
 %T=T+Trho; %Comment out in the exact case
 [U,S,V]=svd(T);
 U0=U(:,4:end);
 U1=U(:,1:3);
 V0=V(:,4:end);
 V1=V(:,1:3);
 s=diag(S);
 
 xx=(0:0.001:1-0.001)';

 AA=exp(2*pi*1i*xx*k);
 
 % Approximation by convolution
 p=diag(AA*T*AA')/ (n+1);
 
 % Christoffel
 epsilon=1/(n+1);
 a=@(xx) exp(2*pi*1i*reshape(xx,[length(xx),1])*k);
 q=@(xx,epsi) 1./ ((abs(a(xx)*U1).^2)*(1./(s(1:3))) + 1/epsi*sum(abs(a(xx)*U0).^2,2));
 q_eval=q(xx,epsilon);
 
 N=n+1;
 
 C_fix=(N*integral(@(x) q(x,epsilon),0,1,'ArrayValued',true)-epsilon).^(-1);
 
 % Signal polynomial
 p1=1/N*sum(abs(AA*U1).^2,2);
 
 % Upper bound on Christoffel function 
 q_approx=epsilon./(N*(1-p1)+epsilon/c(3)*N/max(svd(A)).^2);
 
 % Lower bound on Christoffel function
 q_lower=epsilon./(N*(1-p1)+epsilon/c(1)*N/min(svd(A)).^2);
 
 figure(1)
 plot(xx,real(p)/N,'b-.')
 hold on
 plot(xx,real(q_eval),'g:','LineWidth',4)
 hold on
 plot(xx,real(q_approx),'r-', xx,real(q_lower), 'm--')
 hold on 
 pl=stem(x,c,'kd');
 pl.MarkerSize = 12;
 axis([0,1,0,7/12])
 yticks([epsilon/N c'])
 yticklabels({'$\varepsilon/N$','$\alpha_1 $','$\alpha_2 $','$\alpha_3 $'})
 ax=gca;
 ax.TickLabelInterpreter ='Latex';
 hold off
 legend('$N^{-1} \cdot F_{n}*\mu$','$q_{\varepsilon,n}$',...
     'Upper bound from Lemma 3.3.5','Lower bound from Lemma 3.3.5','$\mu$','Interpreter','Latex',...
     'Location','Northwest','Fontsize',12)
end

%% Analyse influence of noise and choice of epsilon
N_max=24;
eps_number=100;
eps=1e-12+linspace(0,6,eps_number);
T_tries=100;
T_values=linspace(0,6,T_tries);
h=ones(eps_number,T_tries);
test=zeros(eps_number,T_tries);
for n=N_max:N_max
    N=n+1;
    for l=1:T_tries
        k=0:n;
        A=exp(2*pi*1i*x*k);
        T=A'*diag(c)*A;
        rho=randn(n+1,1)+1i*randn(n+1,1);
        Trho=toeplitz(rho); %randn(size(T))+1i*randn(size(T));
        %Trho=Trho-diag(Trho); % zero noise on absolute mass
        %Trho=0.5*(Trho'+Trho);
        Trho=T_values(l)*Trho/norm(Trho);
        T=T+Trho; %Comment out in the exact case
        [U,S,V]=svd(T);
        for i=1:eps_number 
            s=diag(S);          
            r=find(s>eps(i),1,'last');
            if r>0
               if r<N 
                   U0=U(:,r+1:end);
                   V0=V(:,r+1:end);
               elseif r==N
                   U0=zeros(N,1);
                   V0=zeros(N,1);
               end
               U1=U(:,1:r);
               V1=V(:,1:r);
               xx=(0:0.001:1-0.001)';
               AA=exp(2*pi*1i*xx*k);
               a=@(xx) exp(2*pi*1i*reshape(xx,[length(xx),1])*k);
               q=@(xx,epsi) epsi./ ((abs(a(xx)*U1).^2)*(epsi./(s(1:r))) ...
                + sum(abs(a(xx)*U0).^2,2));
               q_eval=q(xx,eps(i));

               q_eval=q_eval-min(q_eval);

               C_total=sum(q_eval);
               C_int=sum(q_eval(1:80))+sum(q_eval(450:580))+sum(q_eval(900:1000));
               h(i,l)=C_int/C_total;
               test(i,l)=sum(abs(diff(q_eval)));
            end
            if l==25 && i==25
                figure(3)
                plot(xx,q_eval,'g:','LineWidth',4)
                legend('$\tilde q_{\varepsilon,n}(x) -\min_y \tilde q_{\varepsilon,n}(y)$'...
                    ,'Interpreter','Latex','Location','Northwest','Fontsize',14)
            end
        end     
    end     
end

%% Find optimal epsilon

[h_values,eps_op]=min(h,[],1);
eps_op=eps(eps_op);
%% Plotting the dependency on noise, eps and n
figure(2)
[EPS,TRho]=meshgrid(eps,T_values);
%stem3(log(EPS)/log(10),Trho,h')
s=surf(EPS',TRho',h,'FaceColor','Interp');
map=colormap('bone');
colormap(flip(map));
colorbar
hold on
stem3(eps_op,T_values,max(max(h))*ones(size(eps_op)),'rx','MarkerSize',18)
hold on
stem3(eps,max(4.14-eps,0),max(max(h))*ones(size(eps)),'m.','MarkerSize',18)
xlabel('Regularisation parameter $\varepsilon$','Interpreter','Latex','Fontsize',14)
xlim([0,6])
ylabel('Noise level $\|T_{\varrho,n}\|_2$','Interpreter','Latex','Fontsize',14)
axis square
s.EdgeColor = 'none';
%% Analysing convergence rate in n
N_max=150;
N_min=5;
h1=zeros(N_max-N_min+1,1);
for n=N_min:N_max
    k=0:n;
    A=exp(2*pi*1i*x*k);
    T=A'*diag(c)*A;
    rho=randn(n+1,1)+1i*randn(n+1,1);
    Trho=toeplitz(rho);
    Trho=1.5*Trho/norm(Trho); % Norm of noise always 1.5
    T=T+Trho; %Comment out in the exact case
    [U,S,V]=svd(T);    
    s=diag(S); 
    epsilon=1.5;
    r=find(s>epsilon,1,'last');
    if r>0
       if r<N 
           U0=U(:,r+1:end);
           V0=V(:,r+1:end);
       elseif r==N
           U0=zeros(N,1);
           V0=zeros(N,1);
       end
       U1=U(:,1:r);
       V1=V(:,1:r);
       xx=(0:0.001:1-0.001)';
       AA=exp(2*pi*1i*xx*k);
       a=@(xx) exp(2*pi*1i*reshape(xx,[length(xx),1])*k);
       q=@(xx,epsi) epsi./ ((abs(a(xx)*U1).^2)*(epsi./(s(1:r))) ...
        + sum(abs(a(xx)*U0).^2,2));
       q_eval=q(xx,epsilon);

       q_eval=q_eval-min(q_eval);

       C_total=sum(q_eval);
       C_int=sum(q_eval(1:80))+sum(q_eval(450:580))+sum(q_eval(900:1000));
       h1(n-N_min+1)=C_int/C_total;
    end  
end

figure(4)
loglog(N_min:N_max,h1,'r--',N_min:N_max,(N_min:N_max).^(-3/2),'b:')
legend('Fraction of mass outside support','Rate of order $n^{-3/2}$',...
    'Interpreter','Latex','Fontsize',14)
xlabel('Bandlimit $n$','Interpreter','Latex')