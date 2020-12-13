function [] = run_hmm_pipeline(dataset_name,data,T,num_pcs,sampling_freq,options)
tic
%%
disp('Running HMM...')
pause(2)
disp('Please create a separate folder for the data you running it on')
disp('This will create a model and produce gamma output')
pause(2)
disp('To avoid confusion run it within a folder for the specific data set...')


%%
if ~isempty(options)
    K_states = options.K;
else
    K_states = input('SUPPLY THE APPROPRIATE NUMBER OF STATES!!');
end
for states = K_states

    all = 1;  

        filename = dataset_name;
    OutputPath = [pwd filesep];
%% Load the specific file name
% Works with a single filename
% data = load(strcat(database,filename));

%% For running only MEG dataset
% data = data.MEG_region_signals;
% run_combined = 0;
% run_tag = 'MEG'

%% For running combined lfp MEG data
% S = data.S;
% data = data.total_signals;
% run_combined = 1;
% run_tag = 'TOTAL';

%% For running LFP only
% data = data.electrodes;
% run_combined = 0;
% run_tag = 'LFP';

%% Run wih all the data 
% data = data.total_Data;
% clear total_Data
 run_combined = 0;
 run_tag = 'Motor_cortex_LFP_all';
% T = cellfun(@length,data,'UniformOutput',0);
% all = 1;

%% For pca
pca_done = 1;

%% For MAR
run_mar = 0;
if run_mar
    run_mar_tag = 'MAR';
    options.order = 11; % (sampling_freq/oreder = smallest_detectable_frequency) set to detect ateast beta
else
    run_mar_tag = 'NO_MAR';
    options.order = 0; % (sampling_freq/oreder = smallest_detectable_frequency) set to detect ateast beta
end
%% Sign flip correction
% disp('Running sign flip correction on data...')
% disp('You can pause and stop execution if you dont want to flip')
% pause(5)
% options.maxlag = 4;
sign_flag = 1;
% [flips,scorepath] = findflip(data,T,options);
% data = flipdata(data,T,flips);
% save('data_sign_corrected','data','-v7.3')
% clear options

%% Hmm options
if isempty(options)
    
    options.K = K_states ;%input('enter desired num of states');
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
    
    
end


%% 

% Saving options
clck = char(datetime);
clck = strrep(clck,':','_');
clck = strrep(clck,' ','_');
clck = strrep(clck,'-','_');
save(strcat(OutputPath,filename,'_',clck,'_','options_file'),'options')
OPTIONS_NAME = strcat(OutputPath,filename,'_',clck,'_','options_file');
save('OPTIONS_NAME','OPTIONS_NAME');
% data to be fed to hmm is TIME X VARIABLES (tXN)
if all
    [hmm, Gamma, Xi, vpath, GammaInit, residuals, fehist, feterms, rho] = hmmmar (data,T,options);
end

if ~all
    T = size(data,2);
    [hmm, Gamma, Xi, vpath, GammaInit, residuals, fehist, feterms, rho] = hmmmar (data',T,options);
end

fe = hmmfe(data,T,hmm,Gamma,Xi);
HMM_model.hmm = hmm;
HMM_model.Gamma = Gamma;
HMM_model.Xi = Xi;
HMM_model.vpath = vpath;
HMM_model.GammaInit = GammaInit;
HMM_model.residuals = residuals;
HMM_model.fehist = fehist;
HMM_model.feterms = feterms;
HMM_model.rho = rho;
HMM_model.sign_corrected = sign_flag;
HMM_model.dataset_name = dataset_name;
HMM_model.dataset_location = pwd;
HMM_model.free_energy = fe;

if ~pca_done
    save(strcat(OutputPath,filename,'_',clck,'_','HMM_model','_',run_mar_tag,'_',run_tag,'_',options.embedd_lag_tags),...
    'HMM_model','-v7.3')
else
    save(strcat(OutputPath,filename,'_',clck,'_','HMM_model_pca','_',run_mar_tag,'_',run_tag,'_',options.embedd_lag_tags),...
        'HMM_model','-v7.3')
    
    MODEL_NAME = strcat(OutputPath,filename,'_',clck,'_','HMM_model_pca','_',run_mar_tag,'_',run_tag,'_',options.embedd_lag_tags);
    save('MODEL_NAME','MODEL_NAME')
    
end

end_time = toc;

save(strcat(OutputPath,filename,'_',clck,'_','total_run_time'),'end_time')

end

end

