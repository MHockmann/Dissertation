%% Plotting Wasserstein rates
% Load results obtained with ExampleProxy.m (stepsize=2)
WF=readmatrix('results/W1_prox0-F_d-15_n1-160.csv'); % Results for Fejer
WK2=readmatrix('results/W1_prox0-K2_d-15_n1-160.csv'); % Jackson kernel 
Wp1=readmatrix('results/W1_prox1--_d-15_n1-160.csv'); % signal polynomial
n_max=160;
n_min=1;
n=2:2:n_max;
% Plotting
loglog(n,WF,'b:','LineWidth',4)
hold on
loglog(n,WK2,'r--',n,Wp1,'m-')
hold on
loglog(n,2/pi^2*(log(n+1)+3)./n,'k',n,0.75./n,'k-.')
hLg=legend('$W_1(F_n*\mu,\mu)$','$W_1(J_n*\mu,\mu)$','$W_1(p_{1,n},\tilde\mu)$',...
    '$\frac{2}{\pi^2} \cdot \frac{\log(n+1)+3}{n}$','$\frac34 \cdot n^{-1}$','Interpreter','Latex',...
    'Fontsize',16);
pos=hLg.Position;        % retrieve existing position
pos(4)=1.25*pos(4);       % increase width value 50% in position 4-vector
hLg.Position=pos;        % set new position
xlim([1.5,165])
xlabel('$n$','Interpreter','Latex',...
    'Fontsize',16)
ylabel('$W_1$','Interpreter','Latex',...
    'Fontsize',16)