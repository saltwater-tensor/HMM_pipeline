
%% TEMPORAL PERMUTATION TESTING
clear
clear test
load('C:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis\OFF\Whole_brain_stn_lfp_medication_OFF_06_Jan_2020_18_39_55_HMM_model_pca_NO_MAR_Motor_cortex_LFP_all_embed_lags\Intervals.mat')
% lifetimesoff = LifeTimes;
lifetimesoff = Intervals;
states_to_test_off = [1,2,3];
% within_state(lifetimesoff,states_to_test_off,'OFF')
clearvars -except lifetimesoff
load('C:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis\ON\Whole_brain_stn_lfp_medication_ON_06_Jan_2020_18_45_18_HMM_model_pca_NO_MAR_Motor_cortex_LFP_all_embed_lags\Intervals.mat')
% lifetimeson = LifeTimes;
lifetimeson = Intervals;
states_to_test_on = [2,4,1];
% within_state(lifetimeson,states_to_test_on,'ON')

clearvars -except lifetimeson lifetimesoff
offstate = input('Enter off state');
onstate = input('Enter on state');
test{1,1} = cell2mat(lifetimesoff(:,offstate)');
test{1,2} = cell2mat(lifetimeson(:,onstate)');
state_int = test;

% test{1,1}(find(test{1,1} < 100)) = [];
% test{1,2}(find(test{1,2} < 100)) = [];
% test{1,1}(find(test{1,1} > 1000)) = [];
% test{1,2}(find(test{1,2} > 1000)) = [];

figure(1)
hold on
histogram(test{1,1})
histogram(test{1,2})

[h_off_lesser_than_on,p_off_lesser_than_on] = ttest2(test{1,1},test{1,2},'VarType','unequal','Tail','left') %off < on
[h_off_greater_than_on,p_off_greater_than_on] = ttest2(test{1,1},test{1,2},'VarType','unequal','Tail','right') %off > on

P = NaN(2,2);
P(2,1) = p_off_greater_than_on; %on less than off
P(1,2) = p_off_lesser_than_on;  %on greater than off

mean_off_coh = mean(test{1,1});
mean_on_coh = mean(test{1,2});
stderror_off = std( test{1,1} ) / sqrt( length( test{1,1} ));
stderror_on = std( test{1,2} ) / sqrt( length( test{1,2} ));

Tt = {test{1,1}',test{1,2}'};
nperms = 500;

for st1 = 1:1:2
    sel_1 = Tt{1,st1};
    for st2 = 1:1:2
        sel_2 = Tt{1,st2};
        
%         sel_1 = T_1;
%         sel_2 = T_4;
        
        combined = [sel_1;sel_2];
        r = randperm(numel(combined));
        combined = combined(r);
        T_obs(st1,st2) = mean(sel_1) - mean(sel_2);
        diff = [];
        D = [];
        subset1 = [];
        subset2 = [];
        
        for n = 1:nperms
        
        % Without replacement  two unique susbets
        %     idx = randperm(numel(combined));
        %     subset1 = combined(1:length(sel_1));
        %     subset2 = combined(length(sel_1)+1:end);
        
        % With replacement two subsets
        subset1 = datasample(combined,length(sel_1),'Replace',false);
        subset2 = datasample(combined,length(sel_2),'Replace',false);
        diff(n) = mean(subset1) - mean(subset2);
        subset1 = [];
        subset2 = [];
        
        end
        h = histogram(diff);
        h.Normalization = 'probability';
        
        D_tail_neg = (diff) > (T_obs(st1,st2));
        D_tail_pos = diff < T_obs(st1,st2);
        prob_neg = sum(D_tail_neg)./length(diff);
        prob_pos = sum(D_tail_pos)./length(diff);
        
        D_neg(st1,st2) = prob_neg;
        D_pos(st1,st2) = prob_pos;
        
        % T_obs > all_possible_shuffle_differences
        % Then that means that our original observed difference is
        % meaningful 
        two_tailed = abs(T_obs(st1,st2)) >= abs(diff);
        two_tail_pval(st1,st2) = 1 - (sum(two_tailed)./length(two_tailed));
    end
end
test_str_1 = 'h_off_greater_than_on is for OFF > ON';
test_str_2 = 'D_neg : observed difference between ith row and jth column is lesser than the permutation sampled difference';
save(['LifeTimes_states_offvson_' num2str(offstate) 'vs' num2str(onstate)],'h_off_greater_than_on','h_off_lesser_than_on','p_off_lesser_than_on','p_off_greater_than_on'...
    ,'D_neg','D_pos','two_tail_pval','test_str_1','test_str_2','mean_off','mean_on','stderror_off','stderror_on')

close all
fig_handle = figure(2);
h = heatmap(P);
h.XLabel = 'HMM condition'; %columns
h.YLabel = 'HMM condition'; %rows
h.Title = ['ROW < COLUMN' num2str(offstate) 'vs' num2str(onstate)];
h.XDisplayLabels = {'OFF','ON'};
h.YDisplayLabels = {'OFF','ON'};
saveas(fig_handle,['map_states_offvson_' num2str(offstate) 'vs' num2str(onstate)]);
saveas(fig_handle,['map_states_offvson_' num2str(offstate) 'vs' num2str(onstate)],'png');
close(fig_handle)


% fig_handle = figure(2);
% h = heatmap(D_pos);
% h.XLabel = 'HMM condition'; %columns
% h.YLabel = 'HMM condition'; %rows
% h.Title = ['D positive ROW > COLUMN' num2str(offstate) 'vs' num2str(onstate)];
% h.XDisplayLabels = {'OFF','ON'};
% h.YDisplayLabels = {'OFF','ON'};
% saveas(fig_handle,['D_pos_map_states_offvson_' num2str(offstate) 'vs' num2str(onstate)]);
% saveas(fig_handle,['D_pos_map_states_offvson_' num2str(offstate) 'vs' num2str(onstate)],'png');
% close(fig_handle)


%%
% OFF Dataset
load('C:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis\OFF\Whole_brain_stn_lfp_medication_OFF_06_Jan_2020_18_39_55_HMM_model_pca_NO_MAR_Motor_cortex_LFP_all_embed_lags\fitmt_subj_fact_4b_filtered1to45.mat')
OFF_subject_specific_spectral = fitmt_subj_fact_4b;
clear fitmt_subj_fact_4b
% ON Dataset
load('C:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis\ON\Whole_brain_stn_lfp_medication_ON_06_Jan_2020_18_45_18_HMM_model_pca_NO_MAR_Motor_cortex_LFP_all_embed_lags\fitmt_subj_fact_4b_filtered1to45_type1.mat')
ON_subject_specific_spectral = fitmt_subj_fact_4b;
clear fitmt_subj_fact_4b
factor_number = 1;

mask = ones(48,48);

contacts = [1,2,3,4,5,6];
frontal = [7,8,11,12,21,22,23,24,25,26,33,34,35,36];
medial_PFC = [1,2,15,16];
temporal = [17,18,39,40];
sensory_motor = [27,28,29,30];
parietal = [5,6,19,20,31,32,37,38,41,42];
visual = [3,4,9,10,13,14];


%---------Select specific connections you want to run tests on
rows = contacts;
cols = sensory_motor; %1:6;
mask(rows,cols) = 0;

%Just for safety
mask(cols,rows) = 0;
mask = ~mask;
%------------------

mask = triu(mask,1);
mask = double(mask);
mask(find(mask == 0)) = NaN;

%Now mask 

OFF_coh =[];
OFF_psd = [];
ON_coh = [];
ON_psd = [];

%Collect across all
for offstate = 1:1:6
    for subject_number = 1:1:17
    
        data_off = mask.* squeeze(OFF_subject_specific_spectral...
            {subject_number, 1}.state(offstate).coh(factor_number,:,:)) ; 
        data_off = data_off(:);
        data_off = data_off(find(~isnan(data_off)));
        OFF_coh = [data_off;OFF_coh];
        
        %Pwr
        psd_off = diag(squeeze(OFF_subject_specific_spectral...
            {subject_number, 1}.state(offstate).psd(factor_number,:,:))) ;
        
        if (length(rows) == length(cols))
            
           if sum(rows-cols) == 0
               psd_off = psd_off([rows]);
           else
               psd_off = psd_off([rows,cols]);
           end
        end
        
        OFF_psd = [psd_off;OFF_psd];
    end
     mean_off_coh(offstate) = mean(OFF_coh);
     err_off_coh(offstate) = std( OFF_coh ) / sqrt( length( OFF_coh ));
     
     mean_off_psd(offstate) = mean(OFF_psd);
     err_off_psd(offstate) = std( OFF_psd ) / sqrt( length( OFF_psd ));
     
     All_off_coh(:,offstate) = OFF_coh;
     All_off_psd(:,offstate) = OFF_psd;
     OFF_coh = [];
     OFF_psd = [];
end


for onstate = 1:1:6
    for subject_number = 1:1:16
    
        data_on = mask.*squeeze(ON_subject_specific_spectral...
            {subject_number, 1}.state(onstate).coh(factor_number,:,:)) ; 
        data_on = data_on(:);
        data_on = data_on(find(~isnan(data_on)));
        ON_coh = [data_on;ON_coh];
       
        %Pwr
        psd_on = diag(squeeze(ON_subject_specific_spectral...
            {subject_number, 1}.state(onstate).psd(factor_number,:,:))) ;
        
        if (length(rows) == length(cols))
            
           if sum(rows-cols) == 0
                psd_on = psd_on([rows]);
           else
               psd_on = psd_on([rows,cols]);
           end
        end
        
        ON_psd = [psd_on;ON_psd];
    end
     mean_on_coh(onstate) = mean(ON_coh);
     err_on_coh(onstate) = std( ON_coh) / sqrt( length(ON_coh ));
    
     mean_on_psd(onstate) = mean(ON_psd);
     err_on_psd(onstate) = std( ON_psd ) / sqrt( length( ON_psd ));
     
     All_on_coh(:,onstate) = ON_coh;
     All_on_psd(:,onstate) = ON_psd;
     ON_coh = [];
     ON_psd = [];
end

% [h_off_lesser_than_on,p_off_lesser_than_on] = ttest2(OFF_coh,ON_coh,'VarType','unequal','Tail','left') %off < on
% [h_off_greater_than_on,p_off_greater_than_on] = ttest2(OFF_coh,ON_coh,'VarType','unequal','Tail','right') %off > on

% COHERENCE

m_off = mean_off_coh(:,[1,2,3]);
error_off = err_off_coh(:,[1,2,3]);

m_on = mean_on_coh(:,[2,4,1]);
error_on = err_on_coh(:,[2,4,1]);

fig_handle = figure(3);
hold on
e1 = errorbar(m_off,error_off,'-s','MarkerSize',10,...
    'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',4);

e2 = errorbar(m_on,error_on,'-s','MarkerSize',10,...
    'MarkerEdgeColor','green','MarkerFaceColor','green','LineWidth',4);

gca
ans.XLim = [0.5 4.5];
ans.XTick = [1,2,3];
ans.XTickLabel = {'HighPwr','Comms','Local'};
[h_off_lesser_than_on,p_off_lesser_than_on] = ttest2(All_off_coh(:,[1,2,3]),All_on_coh(:,[2,4,1]),'VarType','unequal','Tail','left') %off < on
[h_off_greater_than_on,p_off_greater_than_on] = ttest2(All_off_coh(:,[1,2,3]),All_on_coh(:,[2,4,1]),'VarType','unequal','Tail','right') %off > on

txt1 = ['p-value lesser = ' num2str(p_off_lesser_than_on(1))];
txt2 = ['p-value lesser = ' num2str(p_off_lesser_than_on(2))];
txt3 = ['p-value lesser = ' num2str(p_off_lesser_than_on(3))];

txt4 = ['p-value greater = ' num2str(p_off_greater_than_on(1))];
txt5 = ['p-value greater = ' num2str(p_off_greater_than_on(2))];
txt6 = ['p-value greater = ' num2str(p_off_greater_than_on(3))];

text([1,2,3],[m_off(1)+0.25,m_off(2)+0.25,m_off(3)+0.25],{txt1,txt2,txt3},'FontWeight','bold')

text([1,2,3],[m_on(1)+0.25,m_on(2)+0.25,m_on(3)+0.25],{txt4,txt5,txt6},'FontWeight','bold')

% PSD

m_off = mean_off_psd(:,[1,2,3]);
error_off = err_off_psd(:,[1,2,3]);

m_on = mean_on_psd(:,[2,4,1]);
error_on = err_on_psd(:,[2,4,1]);

fig_handle = figure(5);
hold on
e1 = errorbar(m_off,error_off,'-s','MarkerSize',10,...
    'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',4);

e2 = errorbar(m_on,error_on,'-s','MarkerSize',10,...
    'MarkerEdgeColor','green','MarkerFaceColor','green','LineWidth',4);

gca
ans.XLim = [0.5 4.5];

[h_off_lesser_than_on,p_off_lesser_than_on] = ttest2(All_off_psd(:,[1,2,3]),All_on_psd(:,[2,4,1]),'VarType','unequal','Tail','left') %off < on
[h_off_greater_than_on,p_off_greater_than_on] = ttest2(All_off_psd(:,[1,2,3]),All_on_psd(:,[2,4,1]),'VarType','unequal','Tail','right') %off > on

txt1 = ['p-value lesser = ' num2str(p_off_lesser_than_on(1))];
txt2 = ['p-value lesser = ' num2str(p_off_lesser_than_on(2))];
txt3 = ['p-value lesser = ' num2str(p_off_lesser_than_on(3))];

txt4 = ['p-value greater = ' num2str(p_off_greater_than_on(1))];
txt5 = ['p-value greater = ' num2str(p_off_greater_than_on(2))];
txt6 = ['p-value greater = ' num2str(p_off_greater_than_on(3))];

text([1,2,3],[m_off(1)+0.25,m_off(2)+0.25,m_off(3)+0.25],{txt1,txt2,txt3},'FontWeight','bold')

text([1,2,3],[m_on(1)+0.25,m_on(2)+0.25,m_on(3)+0.25],{txt4,txt5,txt6},'FontWeight','bold')


%% SPECTRAL PERMUTATION TESTING
%% Load the datasets
% These are either the factored ones or the binned ones
% It is advisable to use the factored ones
load('C:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis\OFF\Whole_brain_stn_lfp_medication_OFF_06_Jan_2020_18_39_55_HMM_model_pca_NO_MAR_Motor_cortex_LFP_all_embed_lags\fitmt_subj_fact_4b_filtered1to45.mat')
OFF_subject_specific_spectral = fitmt_subj_fact_4b;
clear fitmt_subj_fact_4b
% ON Dataset
load('C:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis\ON\Whole_brain_stn_lfp_medication_ON_06_Jan_2020_18_45_18_HMM_model_pca_NO_MAR_Motor_cortex_LFP_all_embed_lags\fitmt_subj_fact_4b_filtered1to45_type1.mat')
ON_subject_specific_spectral = fitmt_subj_fact_4b;
clear fitmt_subj_fact_4b

tests_off = specttest(OFF_subject_specific_spectral,5000,0,0);
tests_on = specttest(ON_subject_specific_spectral,5000,0,0);
significant_off_005 = spectsignificance(tests_off,0.05);
significant_on_005 = spectsignificance(tests_on,0.05);

%% Visualise
load('D:\Abhinav_Sharma\RS_peri_MEGLFP\brainstorm_db\MEG_LFP_peri\anat\@default_subject\tess_cortex_pial_low')
atlas_name = 'Mindboggle';
atlas_number = 6;
if ~strcmpi(atlas_name,Atlas(atlas_number).Name)
   error('Check for valid atlas number') 
end
ROIs = [3,4,5,6,11,12,15,16,19,20,21,22,23,24,25,26,27,28,29,30,...
    33,34,35,36,37,38,41,42,45,46,47,48,51,52,53,54,55,56,57,58,59,60];
myLabel = cell(length(ROIs));
for i =1:1:length(ROIs)
  myLabel{i} = Atlas(atlas_number).Scouts(ROIs(i)).Label;
end

total_parcellations = 6;
% Original ROI indices
contacts = [1,2,3,4,5,6];
frontal = [7,8,11,12,21,22,23,24,25,26,33,34,35,36];
medial_PFC = [1,2,15,16];
temporal = [17,18,39,40];
sensory_motor = [27,28,29,30];
parietal = [5,6,19,20,31,32,37,38,41,42];
visual = [3,4,9,10,13,14];
re_ROIs = [frontal,medial_PFC,temporal,sensory_motor,parietal,visual];
parc_lengths = [length(contacts),length(frontal),length(medial_PFC),length(temporal)...
    ,length(sensory_motor),length(parietal),length(visual)];
% Original indices of the Atlas Scouts
re_ROIs_2 = ROIs(re_ROIs);
% Create custom node labels
myLabel_2 = cell(length(re_ROIs) + 6);
schemaball_size = 10*ones(48,1);
for i = 1:6
      myLabel_2{i} = ['contact' num2str(i)];
end
for i =7:1:(length(ROIs)+6)
  myLabel_2{i} = Atlas(atlas_number).Scouts(re_ROIs_2(i-6)).Label;
end

myLabel_3 = cell(length(re_ROIs));
% For cortex only plots
for i =1:1:(length(ROIs))
  myLabel_3{i} = Atlas(atlas_number).Scouts(re_ROIs_2(i)).Label;
end

% Indices in the 48x48 matrix
re_ROIs_48_index = re_ROIs + 6;

% Node colors for schemaball
% 48 X 3
colors = lines(total_parcellations + 1);
figure(1000)
x = 1:(total_parcellations + 1);
y = repmat(5,[(total_parcellations + 1),1]);
scatter(1:(total_parcellations + 1),repmat(5,[(total_parcellations + 1),1]),...
   'CData',colors,'SizeData',repmat(2000,[(total_parcellations + 1),1]),'MarkerFaceColor'...
   ,'flat')
xlim([-1,8])
parc_names = {'LFPs','frnt','mPFC','tmp','sm','par','vis'}';
t = text(x,y,parc_names);

colr_mat_2 = zeros(48,3);
begin_index = 1;
for tp = 1:(total_parcellations + 1)
    colr_mat_2(begin_index:begin_index+parc_lengths(tp)-1,:) = repmat(colors(tp,:),[length(begin_index:begin_index+parc_lengths(tp)-1),1]);
    begin_index = begin_index + parc_lengths(tp);
end
colr_mat_2_orig = colr_mat_2;

myLabel_2_schemaball = myLabel_2(:,1);

graph_schemaball = [];
colr_mat_2_schm = colr_mat_2_orig;
K = 6;
total_band_num = 3;

for k = 1:1:K
    
  for band = 1:1:total_band_num 
    
    graph_schemaball = double( squeeze(significant_on001.lower_corr.state(k).coh(band,:,:)));
    schemaball(graph_schemaball,myLabel_2_schemaball,[],colr_mat_2_schm,schemaball_size,colr_mat_2_orig,...
                ['State ' num2str(k) ' band ' num2str(band)]);
            
    
  end
end
