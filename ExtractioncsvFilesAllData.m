% Extract data from Calibration experiments and create csv file to be used

clear all, clc, close all;
pathToData = '/Users/lucia/Documents/Projects/H2020/H2020_code/AMIGO2R2016d/Examples/Toggle_switch/Lugagne_model/Examining_Lugagne_experimental_data/Data_from_Lugagne_paper/CyberSwitch-master/Data/Experimental/';
data_to_compare = {'Calibration_1','Calibration_2','Calibration_3','Calibration_4','Calibration_5','Calibration_6',...
                    'PI_1','PI_2','PI_3',...
                    'BangBang_1','BangBang_2',...
                    'DynStim_1','DynStim_2','DynStim_3','DynStim_4','DynStim_5','DynStim_6','DynStim_7','DynStim_8','DynStim_9','DynStim_10','DynStim_11','DynStim_12','DynStim_13','DynStim_14',...
                    'DynStimTooFast'};
%% CALIBRATION %% Load up the datasets from Lugagne et al.
% FittingData are loaded and corrected for the lamp factor
% Only the Inputsbefore structure from FittingInputs is used. tcon, as well
% as the values of the inducers, are extracted from Media.csv for each
% experiment

for data_set=1:6
    data_comp= data_to_compare{data_set};
    disp(['Loading up data from: ' data_comp]);
    data(data_set).name = data_comp;
    cd(strcat(pathToData,data_comp));
    data(data_set).val = load('FittingData.mat');
    data(data_set).ind = load('FittingInputs.mat');
    
    
    % Adapt the levels with the lamp_factor
    data(data_set).val.rfpMothers=data(data_set).val.rfpMothers*data(data_set).val.lamp_mult;
    data(data_set).val.gfpMothers=data(data_set).val.gfpMothers*data(data_set).val.lamp_mult;
    % If time data starts negative, just shift it.
    data(data_set).val.timerfp = data(data_set).val.timerfp - data(data_set).val.timerfp(1);
    data(data_set).val.timegfp = data(data_set).val.timegfp - data(data_set).val.timegfp(1);
    cd('/Users/lucia/Documents/Projects/BayGame')
end
%% REMAINING DATA
for data_set=7:length(data_to_compare)
    data_comp= data_to_compare{data_set};
    disp(['Loading up data from: ' data_comp]);
    data(data_set).name = data_comp;
    cd(strcat(pathToData,data_comp));
    data(data_set).val = load('Data.mat');
    
    data(data_set).val.rfpMothers=data(data_set).val.XP.rfp;
    data(data_set).val.gfpMothers=data(data_set).val.XP.gfp;
    % If time data starts negative, just shift it.
    data(data_set).val.timerfp = data(data_set).val.XP.timepoints - data(data_set).val.XP.timepoints(1);
    data(data_set).val.timegfp = data(data_set).val.timerfp;
    cd('/Users/lucia/Documents/Projects/BayGame')
end
%% Extract switching times, IPTG and aTc values for each experiment

for data_set=1:length(data_to_compare)
    data_comp= data_to_compare{data_set};
    disp(['Loading up data from: ' data_comp]);
    cd(strcat(pathToData,data_comp));
    DataInd = readtable('Media.csv');
    IPTGLev = DataInd.Var1;
    aTcLev = DataInd.Var2;
    SwitchT = DataInd.Var3/60;
    if SwitchT(1)==0
        EXP_data{1,data_set}.IPTGext = IPTGLev(2:end);
        EXP_data{1,data_set}.aTcext = aTcLev(2:end);
        EXP_data{1,data_set}.switchingtimes = SwitchT;
    else 
        EXP_data{1,data_set}.IPTGext = IPTGLev;
        EXP_data{1,data_set}.aTcext = aTcLev;
        EXP_data{1,data_set}.switchingtimes = [0;SwitchT];
    end
end
cd('/Users/lucia/Documents/Projects/BayGame')

%% Note; experiments 18 and 26 had negative values, probably due to a loss of focus. Decided to exclude those points in reconstruction. 
Data.expName = data_to_compare;
i = 0;
for data_set=1:length(Data.expName)
    i = i+1;
    
    positiveSearch = nanmean(data(data_set).val.rfpMothers,2)>0;
    Data.t_samples{1,i} = [data(data_set).val.timerfp(1,positiveSearch)/60; data(data_set).val.timegfp(1,positiveSearch)/60]; % sampling time in minutes
    Data.n_samples{1,i} = [length(Data.t_samples{1,i}) length(Data.t_samples{1,i})]';
    Data.exp_data{1,i} = [nanmean(data(data_set).val.rfpMothers(positiveSearch,:),2)'; nanmean(data(data_set).val.gfpMothers(positiveSearch,:),2)'];
    Data.standard_dev{1,i} = [nanstd(data(data_set).val.rfpMothers(positiveSearch,:),[],2)';nanstd(data(data_set).val.gfpMothers(positiveSearch,:),[],2)'];
   
    if length(EXP_data{1,data_set}.switchingtimes(find(EXP_data{1,data_set}.switchingtimes<Data.t_samples{1,i}(end))))== length(EXP_data{1,data_set}.switchingtimes)
        Data.t_con{1,i} =[[EXP_data{1,data_set}.switchingtimes(1:end-1,1); max(Data.t_samples{1,data_set}(:,end))] [EXP_data{1,data_set}.switchingtimes(1:end-1,1); max(Data.t_samples{1,data_set}(:,end))]]';
    else
        Data.t_con{1,i} = [[EXP_data{1,data_set}.switchingtimes(find(EXP_data{1,data_set}.switchingtimes<Data.t_samples{1,i}(end)));Data.t_samples{1,i}(end)] [EXP_data{1,data_set}.switchingtimes(find(EXP_data{1,data_set}.switchingtimes<Data.t_samples{1,i}(end)));Data.t_samples{1,i}(end)]]'; 
    end
        Data.input{1,i} = [EXP_data{1,data_set}.IPTGext(1:length(Data.t_con{1,i})-1,1) EXP_data{1,data_set}.aTcext(1:length(Data.t_con{1,i})-1,1)]';
    if data_set<7
        Data.Initial_IPTG{1,i} = data(data_set).ind.inputsbefore(1,2);
        Data.Initial_aTc{1,i} = data(data_set).ind.inputsbefore(1,1);
    elseif strcmp(data_to_compare{data_set},'DynStimTooFast')
        Data.Initial_IPTG{1,i} = 0.25;
        Data.Initial_aTc{1,i} = 20;
    else
        Data.Initial_IPTG{1,i} = 1;
        Data.Initial_aTc{1,i} = 0;
    end
        
end

%% Check if the correction is fine
for data_set=1:length(Data.expName)
    Data.t_con{1,data_set}(1,end) == Data.t_samples{1,data_set}(1,end)
end
    
%% Extract a csv file with time, inputs_before and inputs during the experiment for each dataset

for data_set=1:length(Data.expName)
    data_comp= Data.expName{1,data_set};
    disp(['Considering data from: ' data_comp]);
    IPTGLev = Data.input{1,data_set}(1,:);
    aTcLev = Data.input{1,data_set}(2,:);
    SwitchT = Data.t_con{1,data_set}(1,:);
    time_input = [0:1:round(SwitchT(end))]';
    IPTG = zeros(1,length(time_input));
    aTc = zeros(1,length(time_input));

    for i=1:length(IPTGLev)
        if i== length(IPTGLev)
            IPTG(round(SwitchT(i)):1:end) = IPTGLev(i);
            aTc(round(SwitchT(i)):1:end) = aTcLev(i);

        elseif i== 1
            IPTG(1:1:round(SwitchT(i+1))) = IPTGLev(i);
            aTc(1:1:round(SwitchT(i+1))) = aTcLev(i);
        else
            IPTG(round(SwitchT(i)):1:round(SwitchT(i+1))) = IPTGLev(i);
            aTc(round(SwitchT(i)):1:round(SwitchT(i+1))) = aTcLev(i);
            
        end
    end
    
    aTc_pre = Data.Initial_aTc{1,data_set}.*ones(length(aTc),1);
    IPTG_pre = Data.Initial_IPTG{1,data_set}.*ones(length(IPTG),1);
    index = linspace(1,length(time_input),length(time_input));
    rowsi = strread(num2str(index),'%s');
    varNames = {'time','IPTGpre','aTcpre','IPTG','aTc'};
    Data_table= table(time_input,IPTG_pre,aTc_pre,IPTG',aTc','RowNames',rowsi,'VariableNames',varNames);
    cd('/Users/lucia/Documents/Projects/BayGame/ExperimentalDatacsv')

    writetable(Data_table,strcat(data_to_compare{data_set},'_Inputs.csv'));   
end

%% Extract a csv file with time, inputs_before and inputs during the experiment for each dataset using
%% the Events-based representation

for data_set=1:length(Data.expName)
    data_comp= Data.expName{1,data_set};
    disp(['Considering data from: ' data_comp]);
    IPTGLev = Data.input{1,data_set}(1,:);
    aTcLev = Data.input{1,data_set}(2,:);
    SwitchT = Data.t_con{1,data_set}(1,:);
    IPTG = IPTGLev';
    aTc = aTcLev';
    SwitchingTimes = SwitchT(1:end-1)';
    
    IPTG_pre = Data.Initial_IPTG{1,data_set}.*ones(length(IPTG),1);
    aTc_pre = Data.Initial_aTc{1,data_set}.*ones(length(aTc),1);
    FinalTime = SwitchT(end).*ones(length(SwitchingTimes),1);
    
    index = linspace(1,length(IPTG),length(IPTG));
    rowsi = strread(num2str(index),'%s');
    varNames = {'Switchingtimes','FinalTime','IPTGpre','aTcpre','IPTG','aTc'};
    Data_table= table(SwitchingTimes,FinalTime,IPTG_pre,aTc_pre,IPTG,aTc,'RowNames',rowsi,'VariableNames',varNames);
    cd('/Users/lucia/Documents/Projects/BayGame/ExperimentalDatacsv')

    writetable(Data_table,strcat(data_to_compare{data_set},'_Events_Inputs.csv'));   
end
%% Extract a csv file with time (min), mean and std for each fluorescent reporter, for each experiment

for data_set=1:length(Data.expName)
    time_rfp_min = Data.t_samples{1,data_set}(1,:);
    time_gfp_min = Data.t_samples{1,data_set}(2,:);
    gfp_mean = Data.exp_data{1,data_set}(1,:);
    gfp_std = Data.standard_dev{1,data_set}(1,:);
    rfp_mean = Data.exp_data{1,data_set}(2,:);
    rfp_std = Data.standard_dev{1,data_set}(2,:);
    index = linspace(1,length(time_rfp_min),length(time_rfp_min));
    rowsi = strread(num2str(index),'%s');
    varNames = {'timeGFP','GFPmean','GFPstd','timeRFP','RFPmean','RFPstd'};
    Data_rep= table(time_gfp_min',gfp_mean',gfp_std',time_rfp_min',rfp_mean',rfp_std','RowNames',rowsi,'VariableNames',varNames);
    writetable(Data_rep,strcat(data_to_compare{data_set},'_Observables.csv'));    

end
%%

% %% Extract a csv file with time aTc and IPTG for each experiment
% for data_set=1:length(data_to_compare)
%     data_comp= data_to_compare{data_set};
%     disp(['Loading up data from: ' data_comp]);
%     cd(strcat(pathToData,data_comp));
%     DataInd = readtable('Media.csv');
%     IPTGLev = DataInd.Var1;
%     aTcLev = DataInd.Var2;
%     SwitchT = DataInd.Var3;
%     IPTG = [];
%     aTc = [];
%     for i=1:length(IPTGLev)
%         if i== length(IPTGLev)
%             t = max([data(data_set).val.timerfp(end)/60,data(data_set).val.timegfp(end)/60]); % last acquisition time point 
%             IPTG = [IPTG repmat(IPTGLev(i),[1,length(length(IPTG):1:round(t))])];
%             aTc = [aTc repmat(aTcLev(i),[1,length(length(aTc):1:round(t))])];
% 
%         else  
%             IPTG = [IPTG repmat(IPTGLev(i),[1,length(length(IPTG):1:round(SwitchT(i)/60))])];
%             aTc = [aTc repmat(aTcLev(i),[1,length(length(aTc):1:round(SwitchT(i)/60))])];
%         end
%     end
%     time_input = [0:1:round(t)]';
%     aTc_pre = data(data_set).ind.inputsbefore(1,1).*ones(length(aTc),1);
%     IPTG_pre = data(data_set).ind.inputsbefore(1,2).*ones(length(IPTG),1);
%     index = linspace(1,length(time_input),length(time_input));
%     rowsi = strread(num2str(index),'%s');
%     Data= table(time_input,IPTG_pre,aTc_pre,IPTG',aTc','RowNames',rowsi);
%     cd('/Users/lucia/Documents/Projects/H2020/H2020_code/AMIGO2R2016d/Examples/Toggle_switch/ParameterEstimation/ExtractionAllExperimentalData')
% 
%     writetable(Data,strcat(data_to_compare{data_set},'_Inputs.csv'));   
%   
%     
% end
% 
% 
% 
% %% Extract a csv file with time (min), mean and std for each fluorescent reporter, for each experiment
% GFP_start_mean = [];
% RFP_start_mean = [];
% GFP_start_std = [];
% RFP_start_std = [];
% for data_set=1:length(data)
%     time_rfp_min = data(data_set).val.timerfp/60;
%     time_gfp_min = data(data_set).val.timegfp/60;
%     gfp_mean = nanmean(data(data_set).val.gfpMothers,2);
%     gfp_std = nanstd(data(data_set).val.gfpMothers,[],2);
%     rfp_mean = nanmean(data(data_set).val.rfpMothers,2);
%     rfp_std = nanstd(data(data_set).val.rfpMothers,[],2);
%     index = linspace(1,length(time_rfp_min),length(time_rfp_min));
%     rowsi = strread(num2str(index),'%s');
%     Data= table(time_gfp_min',gfp_mean,gfp_std,time_rfp_min',rfp_mean,rfp_std,'RowNames',rowsi);
%     writetable(Data,strcat(data_to_compare{data_set},'_Expression.csv'));    
%     GFP_start_mean = [GFP_start_mean,gfp_mean(1)];
%     GFP_start_std = [GFP_start_std,gfp_std(1)];
%     RFP_start_mean = [RFP_start_mean,rfp_mean(1)];
%     RFP_start_std = [RFP_start_std,rfp_std(1)];
% end
% %% Plot to verify which are the initial conditions
% figure(1);
% errorbar(1:1:length(GFP_start_mean),GFP_start_mean,GFP_start_std,'o'); hold on;
% figure(2);
% errorbar(1:1:length(RFP_start_mean),RFP_start_mean,RFP_start_std,'o'); hold on;
% % Apparently, the info.txt in the calibration experiment is wrong as to
% % what concerns the ON incubation. You should consider
% % data(i).ind.inputsbefore, where 1: aTC and 2: IPTG
% 


