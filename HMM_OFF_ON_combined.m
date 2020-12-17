
% Load off data non sign flip corrected

OFF_data = load('C:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis\OFF\original_dataset');
OFF_data.data(8) = [];
OFF_T = load('C:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis\OFF\T');
OFF_T.T(8) = [];

S032_data_OFF = load('C:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis_2\OFF\S032_OFF');
S032_T_OFF = load('C:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis_2\OFF\S032_OFF_T');


OFF_data.data(17) = S032_data_OFF.data(1);
OFF_T.T(17) = S032_T_OFF.T(1);



ON_data = load('C:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis\ON\original_dataset');
S032_data_ON = load('C:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis_2\ON\S032_ON');

ON_T = load('C:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis\ON\T');
S032_T_ON = load('C:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis_2\ON\S032_ON_T');

ON_data.data(17) = S032_data_ON.data(1);
ON_T.T(17) = S032_T_ON.T(1);

% Combining random chunks from OFF and ON

OFF = 1:17;
ON = 18:34;

%% Creating combined dataset
% Random seed
rng(1371)
subject_ids = datasample([OFF ON],34);
combined = OFF_data.data;
combined = [combined ON_data.data];
T_comb = OFF_T.T;
T_comb = [T_comb ON_T.T];
for sub = 1:1:length(subject_ids)
    
    if subject_ids(sub) <= 17
        original_dataset_combined(sub,1) = combined(subject_ids(sub),1);
        combined_T(1,sub) = T_comb(subject_ids(sub));
    end
    
    
    if subject_ids(sub) > 17
        original_dataset_combined(sub,1) = combined(subject_ids(sub)-17,2);
        combined_T(1,sub) = T_comb(subject_ids(sub));
    end
    
end

% Time
clearvars -except original_dataset_combined combined_T

%% Sign flip correction for the whole dataset


cd('Q:\HMM_analysis_combined')
load('T')
load('original_dataset')

sign_flip_correction_HMM

%% Run HMM on this dataset with exactly the same options but increase the
% number of states to 12

sampling_freq = 250;

load('T')
load('dataset_sign_flip_corrected')

med_state = 'OFF_ON_combined';
num_pcs = 96; 


run_combined = 0;
 
% For pca
pca_done = 1;

% For MAR
run_mar = 0;

sign_flag = 1;
states = 12; %input('Enter the number of states');
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

