%% This script contains code to generate temporal properties error bar plots


clear
clear test
% load('C:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis\OFF\Whole_brain_stn_lfp_medication_OFF_06_Jan_2020_18_39_55_HMM_model_pca_NO_MAR_Motor_cortex_LFP_all_embed_lags\LifeTimes.mat')
% lifetimesoff = LifeTimes;
% 
% clearvars -except lifetimesoff
% load('C:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis\ON\Whole_brain_stn_lfp_medication_ON_06_Jan_2020_18_45_18_HMM_model_pca_NO_MAR_Motor_cortex_LFP_all_embed_lags\LifeTimes.mat')
% lifetimeson = LifeTimes;

%% Path to temporal properties data
load('C:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis_2\OFF\LifeTimes_with_S032')
 
load('C:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis_2\ON\LifeTimes_with_S032')

%% Error bar plot
offstate = [1,2,3];
onstate = [2,4,1];
offinx = [1,3,5];
oninx = [2,4,6];

for j = 1:1:3
    
    test{1,1} = cell2mat(lifetimesoff(:,offstate(j))');
    test{1,2} = cell2mat(lifetimeson(:,onstate(j))');
    test{1,1}(find(test{1,1} < 100)) = [];
    test{1,2}(find(test{1,2} < 100)) = [];
%     test{1,1}(find(test{1,1} > 1000)) = [];
%     test{1,2}(find(test{1,2} > 1000)) = [];
    mean_off = mean(test{1,1});
    mean_on = mean(test{1,2});
    stderror_off = std( test{1,1} ) / sqrt( length( test{1,1} ));
    stderror_on = std( test{1,2} ) / sqrt( length( test{1,2} ));
    M_off(j) =  mean_off;
    M_on(j) = mean_on;
    err_off(j) = stderror_off;
    err_on(j) = stderror_on;
    
end

number = 1;
fig_handle = figure(1);
hold on
pl = subplot(1,3,number);
pl.LineWidth = 1.5;
pl.XLim = [0.5 2.5];
max_val = max(M_off(number)+ err_off(number), M_on(number)+ err_on(number));
max_val = ceil(max_val) + 1;
min_val = min(M_off(number)- err_off(number), M_on(number)- err_on(number));
min_val = floor(min_val) - 2;
numticks = 4;
ytickval = linspace(min_val,max_val,numticks);
pl.YTick = ceil(ytickval);

% pl.YLim = [min_val max_val];
pl.YLim = [155 175];

pl.XTick = [1, 2];
pl.XTickLabel = {'OFF', 'ON'};
% pl.FontWeight = 'bold';
pl.FontSize = 10;
pl.TickLength = [0.03,0.025];
hold on

e1 = errorbar(1,M_off(number),err_off(number),'-s','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',1.5,'CapSize',8);

e2 = errorbar(2,M_on(number),err_on(number),'-s','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',1.5,'CapSize',8);

if number == 3
a2 = axes();
a2.Position = [0.782499 0.2647619 0.083571428571429 0.223398864132803];
a2.XTick = [];
errorbar(1,M_off(number),err_off(number),'-s','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',1.5,'CapSize',8);
end


%% For testing fractional occupancy difference between state across subjects

% load('C:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis\OFF\Whole_brain_stn_lfp_medication_OFF_06_Jan_2020_18_39_55_HMM_model_pca_NO_MAR_Motor_cortex_LFP_all_embed_lags\FO.mat')
% offFO = FO;
% clear FO
% load('C:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis\ON\Whole_brain_stn_lfp_medication_ON_06_Jan_2020_18_45_18_HMM_model_pca_NO_MAR_Motor_cortex_LFP_all_embed_lags\FO.mat')
% onFO = FO;
% clear FO

load('C:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis_2\OFF\FO_with_S032')

load('C:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis_2\ON\FO_with_S032')

%%

M_off = mean(offFO);
err_off = std(offFO)/sqrt(length(offFO));
% errorbar(M_off,err_off)
M_off = M_off(:,[1,2,3]);
err_off = err_off(:,[1,2,3]);

M_on = mean(onFO);
M_on = M_on(:,[2,4,1]);
err_on = std(onFO)/sqrt(length(onFO));
err_on = err_on(:,[2,4,1]);
% errorbar(M_on,err_on)


% fig_handle = figure(3);
% hold on
% e1 = errorbar(m_off,err_off,'-s','MarkerSize',10,...
%     'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',4);
% 
% e2 = errorbar(m_on,err_on,'-s','MarkerSize',10,...
%     'MarkerEdgeColor','green','MarkerFaceColor','green','LineWidth',4);

number = 2;
fig_handle = figure(2);
hold on
pl = subplot(1,3,number);
pl.LineWidth = 1.5;
pl.XLim = [0.5 2.5];
max_val = max(M_off(number)+ err_off(number), M_on(number)+ err_on(number));
max_val = (max_val) + 0.01;
min_val = min(M_off(number)- err_off(number), M_on(number)- err_on(number));
min_val = (min_val) - 0.01;
numticks = 4;
ytickval = linspace(min_val,max_val,numticks);
pl.YTick = (ytickval);
pl.YLim = [min_val max_val];
pl.XTick = [1, 2];
pl.XTickLabel = {'OFF', 'ON'};
% pl.FontWeight = 'bold';
pl.FontSize = 10;
pl.TickLength = [0.03,0.025];
hold on

e1 = errorbar(1,M_off(number),err_off(number),'-s','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',1.5,'CapSize',8);

e2 = errorbar(2,M_on(number),err_on(number),'-s','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',1.5,'CapSize',8);

% saveas(fig_handle,['error_plot_FO']);
% saveas(fig_handle,['error_plot_FO'],'png');
% close(fig_handle)

%% For testing intervals difference between state across subjects

% load('C:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis\OFF\Whole_brain_stn_lfp_medication_OFF_06_Jan_2020_18_39_55_HMM_model_pca_NO_MAR_Motor_cortex_LFP_all_embed_lags\Intervals.mat')
% intervalsoff = Intervals;
% 
% load('C:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis\ON\Whole_brain_stn_lfp_medication_ON_06_Jan_2020_18_45_18_HMM_model_pca_NO_MAR_Motor_cortex_LFP_all_embed_lags\Intervals.mat')
% intervalson = Intervals;

load('C:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis_2\OFF\Intervals_with_S032')

load('C:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis_2\ON\Intervals_with_S032')

%% 
% Error bar plot
offstate = [1,2,3];
onstate = [2,4,1];
offinx = [1,3,5];
oninx = [2,4,6];

for j = 1:1:3
    
    test{1,1} = cell2mat(intervalsoff(:,offstate(j))');
    test{1,2} = cell2mat(intervalson(:,onstate(j))');
%     test{1,1}(find(test{1,1} < 100)) = [];
%     test{1,2}(find(test{1,2} < 100)) = [];
%     test{1,1}(find(test{1,1} > 1000)) = [];
%     test{1,2}(find(test{1,2} > 1000)) = [];
    mean_off = mean(test{1,1});
    mean_on = mean(test{1,2});
    stderror_off = std( test{1,1} ) / sqrt( length( test{1,1} ));
    stderror_on = std( test{1,2} ) / sqrt( length( test{1,2} ));
    M_off(j) =  mean_off;
    M_on(j) = mean_on;
    err_off(j) = stderror_off;
    err_on(j) = stderror_on;
    
end

number = 1;
fig_handle = figure(3);
hold on
pl = subplot(1,3,number);
pl.LineWidth = 1.5;
pl.XLim = [0.5 2.5];
max_val = max(M_off(number)+ err_off(number), M_on(number)+ err_on(number));
max_val = ceil(max_val) + 1;
min_val = min(M_off(number)- err_off(number), M_on(number)- err_on(number));
min_val = floor(min_val) - 8;
numticks = 4;
ytickval = linspace(min_val,max_val,numticks);
pl.YTick = ceil(ytickval);
pl.YLim = [min_val max_val];
pl.XTick = [1, 2];
pl.XTickLabel = {'OFF', 'ON'};
% pl.FontWeight = 'bold';
pl.FontSize = 10;
pl.TickLength = [0.03,0.025];
hold on

e1 = errorbar(1,M_off(number),err_off(number),'-s','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',1.5,'CapSize',8);

e2 = errorbar(2,M_on(number),err_on(number),'-s','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',1.5,'CapSize',8);

if number == 1
a2 = axes();
a2.Position = [0.21642857142857 0.212380952380959 0.083571428571429 0.223398864132803];
a2.XTick = [];
errorbar(1,M_off(number),err_off(number),'-s','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',1.5,'CapSize',8);
end
% gca
% ans.XLim = [0.5 6.5];
% ans.YLim = [20 35];
% saveas(fig_handle,['error_plot_LifeTimes']);
% saveas(fig_handle,['error_plot_LifeTimes'],'png');
% close(fig_handle)


%% Repeated measures ANOVA with multiple comparison

%% Fractional occupancy test
FO_comb(:,1) = [offFO(:,1);onFO(:,2)];
FO_comb(:,2) = [offFO(:,2);onFO(:,4)];
FO_comb(:,3) = [offFO(:,3);onFO(:,1)];

Treatment_off = cell(17,1);
Treatment_off(:) ={'off'};
Treatment_on = cell(17,1);
Treatment_on(:) ={'on'};
Treatment = [Treatment_off;Treatment_on];


t = table(Treatment,FO_comb(:,1),FO_comb(:,2),FO_comb(:,3),...
'VariableNames',{'Treatment','hyperda','comms','local'});

Meas = table([1 2 3]','VariableNames',{'States'});

rm = fitrm(t,'hyperda-local~Treatment','WithinDesign',Meas);
ranovatbl = ranova(rm);
States_OFFvsON = multcompare(rm, 'Treatment', 'By', 'States');
Treatment_StatevsState = multcompare(rm, 'States', 'By', 'Treatment');

%% FO with anovan

% hyper-da
hyperda_off = (offFO(1:17,1));
treatment_hyperda_off = cell((length(hyperda_off)),1);
treatment_hyperda_off(:) ={'off'};

hyperda_on = (onFO(1:17,2));
treatment_hyperda_on = cell((length(hyperda_on)),1);
treatment_hyperda_on(:) ={'on'};

hyperda = [hyperda_off;hyperda_on];
states_hyperda = cell((length(hyperda_off)+length(hyperda_on)),1);
states_hyperda(:) = {'hyperda'};
treatment_hyperda = [treatment_hyperda_off;treatment_hyperda_on];

% comms
comms_off = (offFO(1:17,2));
treatment_comms_off = cell((length(comms_off)),1);
treatment_comms_off(:) ={'off'};

comms_on = (onFO(1:17,4));
treatment_comms_on = cell((length(comms_on)),1);
treatment_comms_on(:) ={'on'};

comms = [comms_off;comms_on];
states_comms = cell((length(comms_off)+length(comms_on)),1);
states_comms(:) = {'comms'};
treatment_comms = [treatment_comms_off;treatment_comms_on];

% local
local_off = (offFO(1:17,3));
treatment_local_off = cell((length(local_off)),1);
treatment_local_off(:) ={'off'};

local_on = (onFO(1:17,1));
treatment_local_on = cell((length(local_on)),1);
treatment_local_on(:) ={'on'};

local = [local_off;local_on];
states_local = cell((length(local_off)+length(local_on)),1);
states_local(:) = {'local'};
treatment_local = [treatment_local_off;treatment_local_on];

y = [hyperda;comms;local];
y_states = [states_hyperda;states_comms;states_local];
y_treatment = [treatment_hyperda;treatment_comms;treatment_local];


% Testing
clearvars -except y y_states y_treatment

[p,tbl,stats,terms] = anovan(y,{y_states y_treatment},'model','full','varnames',{'States','Treatment'});

figure(2)
results_all = multcompare(stats,'Dimension',[1 2]);

figure(3)
results_states = multcompare(stats,'Dimension',[1]);

figure(4)
results_medication = multcompare(stats,'Dimension',[2]);

%% Intervals
intervalsoff = cellfun(@transpose,intervalsoff,'UniformOutput',false);
intervalson = cellfun(@transpose,intervalson,'UniformOutput',false);

% hyper-da
hyperda_off = cell2mat(intervalsoff(1:17,1));
treatment_hyperda_off = cell((length(hyperda_off)),1);
treatment_hyperda_off(:) ={'off'};

hyperda_on = cell2mat(intervalson(1:17,2));
treatment_hyperda_on = cell((length(hyperda_on)),1);
treatment_hyperda_on(:) ={'on'};

hyperda = [hyperda_off;hyperda_on];
states_hyperda = cell((length(hyperda_off)+length(hyperda_on)),1);
states_hyperda(:) = {'hyperda'};
treatment_hyperda = [treatment_hyperda_off;treatment_hyperda_on];

% comms
comms_off = cell2mat(intervalsoff(1:17,2));
treatment_comms_off = cell((length(comms_off)),1);
treatment_comms_off(:) ={'off'};

comms_on = cell2mat(intervalson(1:17,4));
treatment_comms_on = cell((length(comms_on)),1);
treatment_comms_on(:) ={'on'};

comms = [comms_off;comms_on];
states_comms = cell((length(comms_off)+length(comms_on)),1);
states_comms(:) = {'comms'};
treatment_comms = [treatment_comms_off;treatment_comms_on];

% local
local_off = cell2mat(intervalsoff(1:17,3));
treatment_local_off = cell((length(local_off)),1);
treatment_local_off(:) ={'off'};

local_on = cell2mat(intervalson(1:17,1));
treatment_local_on = cell((length(local_on)),1);
treatment_local_on(:) ={'on'};

local = [local_off;local_on];
states_local = cell((length(local_off)+length(local_on)),1);
states_local(:) = {'local'};
treatment_local = [treatment_local_off;treatment_local_on];

y = [hyperda;comms;local];
y_states = [states_hyperda;states_comms;states_local];
y_treatment = [treatment_hyperda;treatment_comms;treatment_local];

% Testing
clearvars -except y y_states y_treatment

[p,tbl,stats,terms] = anovan(y,{y_states y_treatment},'model','full','varnames',{'States','Treatment'});

figure(2)
results_all = multcompare(stats,'Dimension',[1 2]);

figure(3)
results_states = multcompare(stats,'Dimension',[1]);

figure(4)
results_medication = multcompare(stats,'Dimension',[2]);

%% Lifetimes
lifetimesoff = cellfun(@transpose,lifetimesoff,'UniformOutput',false);
lifetimeson = cellfun(@transpose,lifetimeson,'UniformOutput',false);

% hyper-da
hyperda_off = cell2mat(lifetimesoff(1:17,1));
treatment_hyperda_off = cell((length(hyperda_off)),1);
treatment_hyperda_off(:) ={'off'};

hyperda_on = cell2mat(lifetimeson(1:17,2));
treatment_hyperda_on = cell((length(hyperda_on)),1);
treatment_hyperda_on(:) ={'on'};

hyperda = [hyperda_off;hyperda_on];
states_hyperda = cell((length(hyperda_off)+length(hyperda_on)),1);
states_hyperda(:) = {'hyperda'};
treatment_hyperda = [treatment_hyperda_off;treatment_hyperda_on];

% comms
comms_off = cell2mat(lifetimesoff(1:17,2));
treatment_comms_off = cell((length(comms_off)),1);
treatment_comms_off(:) ={'off'};

comms_on = cell2mat(lifetimeson(1:17,4));
treatment_comms_on = cell((length(comms_on)),1);
treatment_comms_on(:) ={'on'};

comms = [comms_off;comms_on];
states_comms = cell((length(comms_off)+length(comms_on)),1);
states_comms(:) = {'comms'};
treatment_comms = [treatment_comms_off;treatment_comms_on];

% local
local_off = cell2mat(lifetimesoff(1:17,3));
treatment_local_off = cell((length(local_off)),1);
treatment_local_off(:) ={'off'};

local_on = cell2mat(lifetimeson(1:17,1));
treatment_local_on = cell((length(local_on)),1);
treatment_local_on(:) ={'on'};

local = [local_off;local_on];
states_local = cell((length(local_off)+length(local_on)),1);
states_local(:) = {'local'};
treatment_local = [treatment_local_off;treatment_local_on];

y = [hyperda;comms;local];
y_states = [states_hyperda;states_comms;states_local];
y_treatment = [treatment_hyperda;treatment_comms;treatment_local];

% Testing
clearvars -except y y_states y_treatment

[p,tbl,stats,terms] = anovan(y,{y_states y_treatment},'model','full','varnames',{'States','Treatment'});

figure(2)
results_all = multcompare(stats,'Dimension',[1 2]);

figure(3)
results_states = multcompare(stats,'Dimension',[1]);

figure(4)
results_medication = multcompare(stats,'Dimension',[2]);
