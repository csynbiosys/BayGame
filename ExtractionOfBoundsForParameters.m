% Definition of bounds for the transformed parameter 
%% 
clc, clear all, close all;

%% % Mean estimates for the parameters (same order of magnitude of estimates in Lugagne )
k_Lm0 = 3e-2;
k_Lm = 10;
k_Tm0 = 1e-1;
k_Tm = 1;
k_Lp = 1;
k_Tp = 1;

g_Lm = 1.386e-1;
g_Tm = g_Lm;

g_Lp = 1.65e-2;
g_Tp = g_Lp;

Theta_L = 30;
Theta_T = 30;
Theta_IPTG = 0.5;
Theta_aTc = 50;
n_L = 2;
n_T = 2;
n_aTc = 2;
n_IPTG = 2;

k_IPTG = 4e-2;
k_in_IPTG = 4e-2;
k_out_IPTG = 4e-2;

k_aTc = 1e-1;
k_in_aTc = 1e-1;
k_out_aTc = 1e-1;

%% For the unknown parameters, build normal prior centered on the mean and covering 0.5-1.5 times their value
Params = [k_Lm0*k_Lp k_Lm*k_Lp Theta_T Theta_aTc n_aTc n_T k_Tm0*k_Tp k_Tm*k_Tp Theta_L Theta_IPTG n_IPTG n_L k_IPTG k_in_IPTG k_out_IPTG k_aTc k_in_aTc k_out_aTc];
Start = zeros(1,length(Params)+4);
Stop = Start;
Type = {};
Parameter = {'k_L_pm0','k_L_pm','Theta_T','Theta_aTc','n_aTc','n_T','k_T_pm0','k_T_pm','Theta_L','Theta_IPTG','n_IPTG','n_L','k_IPTG','k_in_IPTG','k_out_IPTG','k_aTc','k_in_aTc','k_out_aTc','g_Lm','g_Tm','g_Lp','g_Tp'};
hill = [5,6,11,12];
for i=1:length(Params)
    if ismember(i,hill)
        Parameter{i}
        Start(i) = 0;
        Stop(i) = 5;        
    else
        Start(i) = 0.1*Params(i);
        Stop(i) = 10*Params(i);
    end
    
end
%% Theta_aTc
Start(4) = 0;
Stop(4) = 100;

% Theta_IPTG
Start(10) = 0; 
Stop(10) = 1;
%%
Start(i+1:end) = [g_Lm, g_Tm,g_Lp,g_Tp];
Stop(i+1:end) = [g_Lm, g_Tm,g_Lp,g_Tp];

%% Creating the table
index = linspace(1,length(Parameter),length(Parameter));
rowsi = strread(num2str(index),'%s');
varNames = {'Parameter_Name','StartBound','StopBound'};

Data= table(Parameter',Start',Stop','RowNames',rowsi,'VariableNames',varNames);
writetable(Data,'BoundsTransformedParameters.csv');    



