
%%



file='/path/to/SSA/calculations.jld2';

%Variances for N=25:25:100
variancesS_25=h5read(file,'/variancesS_25');
variancesS_50=h5read(file,'/variancesS_50');
variancesS_75=h5read(file,'/variancesS_75');
variancesS_100=h5read(file,'/variancesS_100');
variancesSL_25=h5read(file,'/variancesSL_25');
variancesSL_50=h5read(file,'/variancesSL_50');
variancesSL_75=h5read(file,'/variancesSL_75');
variancesSL_100=h5read(file,'/variancesSL_100');
alpha_vec=h5read(file,'/a_vec');
nS=h5read(file,'/nS');
nL=h5read(file,'/nL');

%% 

N_vec=[25,50,75,100];%carrying capacities to be considered

stoichiometricMatrix=[1,-1,1,0;0,1,-1,-2]; %CRN Stoichiometric matrix
Omega=1;%system volume

%memory allocation for variances of S, and S+L
varS_LNA=nan(numel(N_vec),numel(alpha_vec));
varSL_LNA=nan(numel(N_vec),numel(alpha_vec));

tic
for r1=1:numel(N_vec)
    N=N_vec(r1);%carrying capacity
    phiS=nS/N; phiP=nL/N;%homeostatic concentrations 
    %stochastic propensity functions
    propensityFunctions= @(x) [x(1)*(N-x(1))/N,x(1),x(2)*(N-x(1))/N,x(2)*(x(2)-1)/N];
    %deterministic propensity functions
    propensityFunctionsDet=@(x) [x(1)*(1-x(1)),x(1),x(2)*(1-x(1)),x(2)^2];
    disp(r1)
    for r2=1:numel(alpha_vec)
        
        disp(r2)
        %switching rate setting
        k2=alpha_vec(r2);
        %calculates the other rates so in steady state nS and nL are the
        %average number of stem cells
        k3=phiS/phiP*(k2/(1-phiS)-1);
        k4=phiS*(1-phiS)/2/(phiP^2);
        
        %building the RxNet object
        rateVector=[1,k2,k3,k4];
        rxn= RxnNet(stoichiometricMatrix,rateVector,propensityFunctions);
        
        %% LNA
        fluctuationsMatrix=rxn.LNAnumerical([phiS,phiP],N,propensityFunctionsDet);
        varS_LNA(r1,r2)=fluctuationsMatrix(1,1);
        varSL_LNA(r1,r2)=varS_LNA(r1,r2)+fluctuationsMatrix(2,2)+2*fluctuationsMatrix(1,2);

    end
 end
 toc

%This can be commented after one run
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex'); 


close all
figure, 
errorbar(alpha_vec,mean(variancesS_100),std(variancesS_100),'.k','MarkerSize',20),hold on,
plot(alpha_vec,varS_LNA(4,:),'k','LineWidth',2)
errorbar(alpha_vec,mean(variancesS_25),std(variancesS_25),'MarkerSize',20)
plot(alpha_vec,varS_LNA(1,:),'LineWidth',2)
errorbar(alpha_vec,mean(variancesS_50),std(variancesS_50),'.','MarkerSize',20), 
plot(alpha_vec,varS_LNA(2,:),'LineWidth',2)
errorbar(alpha_vec,mean(variancesS_75),std(variancesS_75),'.','MarkerSize',20),
plot(alpha_vec,varS_LNA(3,:),'LineWidth',2)

xlabel('licensing/proliferation rate, $\alpha$','FontSize',24,'Interpreter','latex');
ylabel('$\sigma^2_S$','FontSize',24,'Interpreter','latex');
ax=gca; ax.FontSize=22; ax.LineWidth=1.5;
title('$n_S^*=18 \quad n_L^*=26$','FontSize',22,'Interpreter','latex');
legend('$SSA$','$LNA$','FontSize',20,'Interpreter','Latex')
legend boxoff
xlim([1,10])


figure,
errorbar(alpha_vec,mean(variancesSL_100),std(variancesSL_100),'.k','MarkerSize',20),hold on,
plot(alpha_vec,varSL_LNA(4,:),'k','LineWidth',2)
errorbar(alpha_vec,mean(variancesSL_25),std(variancesSL_25),'MarkerSize',20)
plot(alpha_vec,varSL_LNA(1,:),'LineWidth',2)
errorbar(alpha_vec,mean(variancesSL_50),std(variancesSL_50),'.','MarkerSize',20), 
plot(alpha_vec,varSL_LNA(2,:),'LineWidth',2)
errorbar(alpha_vec,mean(variancesSL_75),std(variancesSL_75),'.','MarkerSize',20),
plot(alpha_vec,varSL_LNA(3,:),'LineWidth',2)
xlabel('licensing/prokliferation rate, $\alpha$','FontSize',24,'Interpreter','latex');
ylabel('$\sigma^2_{SL}$','FontSize',24,'Interpreter','latex');
ax=gca; ax.FontSize=22; ax.LineWidth=1.5;
title('$n_S^*=18 \quad n_L^*=26$','FontSize',22,'Interpreter','latex');
legend('$SSA$','$LNA$','FontSize',20,'Interpreter','Latex')
legend boxoff
xlim([1,10])

