% Definition of priors based on parameter estimates in Lugagne et al.
% Parameters obtained upon fitting on the calibration data in Nat.Com 2017
k_Lm0 = 3.045e-1;
k_Lm = 13.01;
k_Tm0 = 3.313e-1;
k_Tm = 5.055;
k_Lp = 6.606e-1;
k_Tp = 5.098E-1;

g_Lm = 1.386e-1;
g_Tm = g_Lm;

g_Lp = 1.65e-2;
g_Tp = g_Lp;

Theta_L = 124.9;
Theta_T = 76.40;
Theta_IPTG = 2.926e-1;
Theta_aTc = 35.98;
n_L = 2;
n_T = 2.152;
n_aTc = 2;
n_IPTG = 2;

k_IPTG = 4e-2;

%% For the unknown parameters, build normal prior centered on the mean and covering 0.5-1.5 times their value
Params = [k_Lm0*k_Lp k_Lm*k_Lp Theta_T Theta_aTc n_aTc n_T k_Tm0*k_Tp k_Tm*k_Tp Theta_L Theta_IPTG n_IPTG n_L k_IPTG];
Mean = zeros(1,length(Params)+4);
Std = Mean;
Type = {};
for i=1:length(Params)
    start = 0.5*Params(i);
    stop = 1.5*Params(i);
    Mean(i) = Params(i);
    Std(i) = (stop-start)/6;
    Type = [Type; 'unknown'];
    
end
%% Adding info on the fixed parameters
Mean(i+1:end) = [g_Lm, g_Tm,g_Lp,g_Tp];
Type{14,1} = 'fixed';
Type{15,1} = 'fixed';
Type{16,1} = 'fixed';
Type{17,1} = 'fixed';
Parameter = {'k_L_pm0','k_L_pm','Theta_T','Theta_aTc','n_aTc','n_T','k_T_pm0','k_T_pm','Theta_L','Theta_IPTG','n_IPTG','n_L','k_IPTG','g_Lm','g_Tm','g_Lp','g_Tp'};

%% Creating the table
index = linspace(1,length(Parameter),length(Parameter));
rowsi = strread(num2str(index),'%s');
varNames = {'Parameter_Name','Type','Mean','Std'};

Data= table(Parameter',Type,Mean',Std','RowNames',rowsi,'VariableNames',varNames);
writetable(Data,'InitialPriors_Lugagne_InducerExchange.csv');    



