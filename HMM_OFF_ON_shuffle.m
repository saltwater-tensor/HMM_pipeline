% For this shuffle analysis we are going to use the dataset created for the
% HMM analysis combined
% 
% cd('Q:\HMM_analysis_combined')
% load('T')
% load('original_dataset')
% 
% %% Split along the half
% % First half is treated as OFF data and second half is treated as ON
% % Sign flip correction is applied separately on this dataset.
% 
% OFF_data = data(1:17,:);
% T_OFF = T(1,1:17);
% ON_data = data(18:34,:);
% T_ON = T(1,18:34);
% cd('Q:\HMM_analysis_shuffled\OFF')
% data = OFF_data;
% T = T_OFF;
% save('original_dataset','data','-v7.3')
% save('T','T')
% clear data T OFF_data T_OFF
% cd('Q:\HMM_analysis_shuffled\ON')
% data = ON_data;
% T = T_ON;
% save('original_dataset','data','-v7.3')
% save('T','T')

%% Go to anyone directroy and run the whole pipeline
% cd('Q:\HMM_analysis_shuffled\OFF')
% load('T')
% load('original_dataset')
% sign_flip_correction_HMM
% clear

%%

sampling_freq = 250;

load('T')
load('dataset_sign_flip_corrected')

med_state = 'ON_shuffle';
num_pcs = 96; 


run_combined = 0;
 
% For pca
pca_done = 1;

% For MAR
run_mar = 0;

sign_flag = 1;
states = 6; %input('Enter the number of states');
options.K = states ;%input('enter desired num of states');
options.Fs = sampling_freq; %input('enter sampling freq in Hz');
if run_mar
    run_mar_tag = 'MAR';
    options.order = 11; % (sampling_freq/oreder = smallest_detectable_frequency) set to detect ateast beta
else
    run_mar_tag = 'NO_MAR';
    options.order = 0; % (sampling_freq/oreder = smallest_detectable_frequency) set to detect ateast beta
end

options.timelag = 2;
options.orderoffset = 1;
options.covtype = 'full';
options.zeromean = 1;
options.DirichletDiag = 10; %default
if pca_done
    options.pca = num_pcs; 
end
options.standardise = 1;
%options.standardise_pc = 1;
options.onpower = 0;
options.filter = [];
options.detrend = 0;
options.downsample = 0;

if run_combined
   options.S = S; %not set will use the default value 
end

options.symmetricprior = 1; %For lfp signals driving cortex
options.dropstates = 1;
options.tol = 1e-5;
options.cyc = 1000;
options.inittype = 'hmmmar';
options.initrep = 5;
options.initcyc = 100;
options.repetitions = 1;
options.updateGamma = 1;
options.decodeGamma = 1;
options.useParallel = 0;
% options.DirStats = OutputPath;

options.embeddedlags = -7:7; %[-10 8 6 4 2 0 2 4 6 8 10];
options.embedd_lag_tags = 'embed_lags';
% Stochastic gradient paramters
options.BIGNbatch = 8; 
options.BIGNinitbatch = options.BIGNbatch ; options.BIGtol = 1e-7; options.BIGcyc = 100; % or 
options.BIGundertol_tostop = 5; options.BIGdelay = 1; 
options.BIGforgetrate = 0.7; options.BIGbase_weights = 0.9;


dataset_name = ['Whole_brain_lfp_stn_medication' med_state];
run_hmm_pipeline(dataset_name,data,T,num_pcs,sampling_freq,options)
save('dataset_name','dataset_name');

clear

%% STEP 4: Calculate HMM basics

load('T')
load('dataset_sign_flip_corrected')
load('sampling_freq')
load('med_state')
load('dataset_name')

%% Load HMM model and options
load('MODEL_NAME')
try
    load(MODEL_NAME)
catch
    string_split = strsplit(MODEL_NAME,'\');
    MODEL_NAME = [MODEL_NAME '\' string_split{end}];
end

load('OPTIONS_NAME')
load(OPTIONS_NAME)


%% Temporal aspects
[vpath] = HMM_model.vpath;
save('viterbi_path','vpath');
LifeTimes = getStateLifeTimes (HMM_model.vpath,T,options,5);
save('LifeTimes','LifeTimes');
Intervals = getStateIntervalTimes (HMM_model.vpath,T,options,5);
save('Intervals','Intervals');
FO = getFractionalOccupancy (HMM_model.Gamma,T,options);
save('FO','FO');
maxFO = getMaxFractionalOccupancy(HMM_model.Gamma,T,options);
save('maxFO','maxFO');
switchingRate =  getSwitchingRate(HMM_model.Gamma,T,options);
save('switchingRate','switchingRate');

%% Distance between states
total_states = HMM_model.hmm.K;
C_all = zeros(total_states,total_states);
for k = 1:1:total_states
    % Takes the low dimensional pcs_dims by pc_dims matrix 
    % projects it back to the high dimensional num_lags X contacts space
    % using PC loading matrix A 
    % Not totally accurate but it makes sense 
    % 90% of variance is or should be captured by pcs atleast
    C1 = getAutoCovMat(HMM_model.hmm,k);
    C_all_covs(:,:,k) = C1;
    for k2 = 1:1:total_states
          C2 = getAutoCovMat(HMM_model.hmm,k2);
          [E1] = eig(C1,C2);
          logE1 = log(E1);
          sqr = logE1 .* logE1;
          sumsqr = sum(sqr);
          sqrtsumsqr = sqrt(sumsqr);
          abssq = abs(sqrtsumsqr);
          dist1 = abssq;
          if k == k2
              C_all(k,k2) = 0;
          else
              C_all(k,k2) = dist1;
          end 
    end   
end

figure(26)
heatmap(C_all,'Colormap',jet)
savefig('state_distances')
close all
save(['covariance_matrices_' med_state],'C_all_covs')

%% STEP 5:  MULTITAPER SPECTRAL ANALYSIS

% Load HMM model and options
load('MODEL_NAME')

%%


load('sampling_freq')
options_mt = struct('Fs',sampling_freq); % Sampling rate - for the 25subj it is 300
options_mt.fpass = [1 45];  % band of frequency you're interested in
options_mt.tapers = [4 7]; % taper specification - leave it with default values
options_mt.p = 0; %0.01; % interval of confidence  
options_mt.win = 2 * sampling_freq; % multitaper window
options_mt.to_do = [1 0]; % turn off pdc
options_mt.order = 0;
options_mt.embeddedlags = -7:7;


% Hz = 250;
% options_mt = struct('Fs',Hz); 
% % options_mt.fpass = [1 45];  % band of frequency you're interested in
% options_mt.tapers = [4 7]; % taper specification - leave it with default values
% options_mt.p = 0; %0.01; % interval of confidence  
% % options_mt.win = 2* Hz; % multitaper window
% options_mt.to_do = [1 0]; % turn off pdc
% options_mt.order = 0;
% options_mt.embeddedlags = -7:7;



% Group level spectra
% load('T')
% load('dataset_sign_flip_corrected')
T = cell2mat(T);
T = T';
data = cell2mat(data);
fitmt_group = hmmspectramt(data,T,HMM_model.Gamma,options_mt);
clear data T
save('fitmt_group_filtered1to45','fitmt_group','-v7.3')

%%
%Subject specific spectra
load('T')
load('dataset_sign_flip_corrected')
N = size(data,1);
% per subject
fitmt_subj = cell(N,1);
d = length(options_mt.embeddedlags) - 1;
acc = 0 ;
N = size(data,1);
Gamma = HMM_model.Gamma;
for n=1:N
    X = data{n};
    gamma = Gamma(acc + (1:(sum(T{n})-length(T{n})*d)),:);
    acc = acc + size(gamma,1);
    fitmt_subj{n} = hmmspectramt(X,T{n},gamma,options_mt);
    fitmt_subj{n}.state = rmfield(fitmt_subj{n}.state,'ipsd');
    fitmt_subj{n}.state = rmfield(fitmt_subj{n}.state,'pcoh');
    fitmt_subj{n}.state = rmfield(fitmt_subj{n}.state,'phase');
    disp(['Subject ' num2str(n)])
end
save('fitmt_subj_filtered1to45','fitmt_subj' ,'-v7.3')

%% Four band analysis
% options_fact = struct();
% options_fact.Ncomp = 4; 
% options_fact.Base = 'coh';
% % We choose to do it in the assymetric way for nnmf factor calculation
% % Also we wont do the robust automated fitting because we want to see what 
% % types of factors do we get 
% 
% % SUBSTEP 1: Get back the matrix created by combining all individual subject
% % level multitaper results
% % No NNMF is performed here
% chan = [];
% chan1 = [1:6] ;%for stn contacts
% chan2 = [7:48];%the remaining cortical locations

%% 

%Load the subject level dataset
% load('fitmt_subj_filtered1to45')
% supply_profiles = [];
% [X] = ...
%     spectdecompose_custom_updated(fitmt_subj,options_fact,chan,chan1,chan2,supply_profiles);
% % save('matrix_to_factor_filtered1to45','X','-v7.3')
% 
% %% NNMF with plot, % Perform decomposition here
% % Once you are satisfied with the profiles you will proceed to substep 2
% figure(9022)
% [supply_profiles,H] = nnmf(X,4,'algorithm','als');
% 
% %%
% load('fitmt_group_filtered1to45')
% 
% subplot1 = @(m,n,p) subtightplot (m, n, p, [0.04 0.5], [0.01 0.03], [0.02 0.02]);
% for pl = 1:1:size(supply_profiles,2)
%    subplot1(size(supply_profiles,2),1,pl)
%    plot(fitmt_group.state(1).f,supply_profiles(:,pl),'Linewidth',3)
% end
% 
% Hz = fitmt_group.state(1).f;
% 
% %%
% savestep = input('Enter row numbers as vector of the supply profiles that you wish to save as factors, REMOVE NOISE FACTOR');
% supply_profiles = supply_profiles(:,savestep);
% save('supply_profiles_factors','supply_profiles')
% 
% %% SUBSTEP 2: Now project the data on a group level and subject level
% % calculated using the factors calculated at a group level above
% load('fitmt_subj_filtered1to45')
% [X,fitmt_subj_fact_4b,fitmt_group_fact_4b,sp_profiles_4b] = ...
%     spectdecompose_custom_updated(fitmt_subj,options_fact,chan,chan1,chan2,supply_profiles);


