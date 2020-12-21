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

number = 3;
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

