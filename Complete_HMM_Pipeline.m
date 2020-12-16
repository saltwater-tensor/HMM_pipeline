%% COMPLETE HMM PIPELINE

% The preprocessinng ,cleaning and projection to source space of the data
% happens separately in brainstorm and those steps are not included here.
% But since projection to default anatomy takes place while creating the
% dataset brainstorm with nogui is invoked
% NOTE: Please take care that within brainstorm the correct database is
% loaded otherwise projection to default anatomy will not work.

% This complete pipeline is a call to scripts and functions please check
% the required list of scripts before running
% This would also require manual check of certain things within those
% scripts

work_dir = input('enter complete path for the place where all the model analysis will be saved');
cd(work_dir)

%% STEP 1: DATASET CREATION
% NOTE : PLEASE CHECK THE ATLAS THAT YOU ARE USING AND THE LIST OF ROIs
% THAT YOU INTEND TO USE. THIS IS CRITICAL
% SCRIPT INVOKED: MEG_lfp_HMM_dataset_creator_ROI.m

atlas_check = input('Did you check the ROIs list line 40 and atlas name line 335 ;1 for yes');
atlas_name_check = input('ENTER the atlas name');
if ~atlas_check
    error('Atlas not checked');
end

MEG_lfp_HMM_dataset_creator_ROI
save('original_dataset','DATA','-v7.3')
save('T','T')
sampling_freq = downsample_rate;
atlas_used = atlas_name;
save('sampling_freq','sampling_freq');
save('atlas_used','atlas_used');
save('original_dataset_subject_level','data_subject','-v7.3')
save('med_state','med_state');
clear

load('original_dataset')
data = DATA;
clear DATA
save('original_dataset','data','-v7.3');
clear

%% STEP 2: SIGN FLIP CORRECTION
load('T')
load('original_dataset')

sign_flip_correction_HMM
clear

%% STEP 3: Run HMM
% 
load('T')
load('dataset_sign_flip_corrected')
load('sampling_freq')
load('med_state')
num_pcs = 96; 

%% Hmm options

run_combined = 0;
 
% For pca
pca_done = 1;

% For MAR
run_mar = 0;

sign_flag = 1;
states = input('Enter the number of states');
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


% warm start
options.hmm = HMM_model.hmm;
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
[vpath] = hmmdecode(data,T,HMM_model.hmm,1);
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
save('fitmt_group_filtered13to30','fitmt_group','-v7.3')

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
save('fitmt_subj_filtered13to30','fitmt_subj' ,'-v7.3')


%% STEP 6: NNMF BASED SPECTRAL ANALYSIS

% This is the most variable and user defined part
% Here you have multiple options to conduct the nnmf analysis
% You choose to focus on specific channels or the entire matrix for
% factorisation
% Also you can perform a robust automated nnmf factorisation or perform
% custom repeats to get factors that you think are desirable and useful
% We will use a custom function here that is very different from the one
% that is supplied with the HMM toolbox 

%% Four band analysis
options_fact = struct();
options_fact.Ncomp = 4; 
options_fact.Base = 'coh';
% We choose to do it in the assymetric way for nnmf factor calculation
% Also we wont do the robust automated fitting because we want to see what 
% types of factors do we get 

%% SUBSTEP 1: Get back the matrix created by combining all individual subject
% level multitaper results
% No NNMF is performed here
chan = [];
chan1 = [1:6] ;%for stn contacts
chan2 = [7:48];%the remaining cortical locations

%%
supply_profiles = [];
[X] = ...
    spectdecompose_custom_updated(fitmt_subj,options_fact,chan,chan1,chan2,supply_profiles);
save('matrix_to_factor_filtered1to45','X','-v7.3')

%% NNMF with plot, % Perform decomposition here
% Once you are satisfied with the profiles you will proceed to substep 2
figure(9022)
[supply_profiles,H] = nnmf(X,4,'algorithm','als');
subplot1 = @(m,n,p) subtightplot (m, n, p, [0.04 0.5], [0.01 0.03], [0.02 0.02]);
for pl = 1:1:size(supply_profiles,2)
   subplot1(size(supply_profiles,2),1,pl)
   plot(fitmt_group.state(1).f,supply_profiles(:,pl),'Linewidth',3)
end
Hz = fitmt_group.state(1).f;

%%
savestep = input('Enter row numbers as vector of the supply profiles that you wish to save as factors, REMOVE NOISE FACTOR');
supply_profiles = supply_profiles(savestep,:);
save('supply_profiles_4b','supply_profiles')

%% SUBSTEP 2: Now project the data on a group level and subject level
% calculated using the factors calculated at a group level above

[X,fitmt_subj_fact_4b,fitmt_group_fact_4b,sp_profiles_4b] = ...
    spectdecompose_custom_updated(fitmt_subj,options_fact,chan,chan1,chan2,supply_profiles);

%% clear supply_profiles
save('fitmt_subj_fact_4b_filtered1to45_type4','fitmt_subj_fact_4b')
save('fitmt_group_fact_4b_filtered1to45_type4','fitmt_group_fact_4b')

%% SUBSTEP 3: Repeat for broadband by changing number of factors to 2
% 
figure(9022)
[supply_profiles,H] = nnmf(X,2,'algorithm','als');
subplot1 = @(m,n,p) subtightplot (m, n, p, [0.04 0.5], [0.01 0.03], [0.02 0.02]);
for pl = 1:1:size(supply_profiles,2)
   subplot(size(supply_profiles,2),1,pl)
   plot(fitmt_group.state(1).f,supply_profiles(:,pl),'Linewidth',3)
end
%%
savestep = input('Enter row numbers as vector of the supply profiles that you wish to save as factors, REMOVE NOISE FACTOR');
supply_profiles = supply_profiles(savestep,:);
%%
save('supply_profiles_wb','supply_profiles')

%% SUBSTEP 4 : Project to 2b
[X,fitmt_subj_fact_2b,fitmt_group_fact_2b,sp_profiles_2b] = ...
    spectdecompose_custom_updated(fitmt_subj,options,chan,chan1,chan2,supply_profiles);

save('fitmt_subj_fact_2b','fitmt_subj_fact_2b')
save('fitmt_group_fact_2b','fitmt_group_fact_2b')
clear supply_profiles


%% -------------------------------- VISUALISATION -------------------------%%




%% STEP 7: CREATING FIGURES
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

%% Node colors for schemaball
% 48 X 3
colors = lines(total_parcellations + 1);
figure(1000)
x = 1:(total_parcellations + 1);
y = repmat(5,[(total_parcellations + 1),1]);
scatter(1:(total_parcellations + 1),repmat(5,[(total_parcellations + 1),1]),'s',...
   'CData',colors,'SizeData',repmat(200,[(total_parcellations + 1),1]),'MarkerFaceColor'...
   ,'flat','MarkerEdgeColor','none')
xlim([-1,8])
set(gca,'visible','off')
parc_names = {'LFPs','frnt','mPFC','tmp','sm','par','vis'}';
% t = text(x,y,parc_names);

colr_mat_2 = zeros(48,3);
begin_index = 1;
for tp = 1:(total_parcellations + 1)
    colr_mat_2(begin_index:begin_index+parc_lengths(tp)-1,:) = repmat(colors(tp,:),[length(begin_index:begin_index+parc_lengths(tp)-1),1]);
    begin_index = begin_index + parc_lengths(tp);
end
colr_mat_2_orig = colr_mat_2;


%% RING PLOTS
% Selecting significant connections
%  Uncorrected test
P_VALUE = input('enter your desired p-value for testing');

K = 6;
ndim = 48;
total_band_num = 3;
cortex = [7;48];
lfp_pos = [1;6];
mask_connections = 0;
within_state = 0;
threshold = 1;
fitmt_group_fact = fitmt_group_fact_4b;
% Mean activation/coherence across all states
if total_band_num == 1
    Band_tag = 'Wideband';
    M = zeros(ndim);
    M_psd = zeros(ndim);
    for k = 1:K
        M = M + squeeze(abs(fitmt_group_fact_4b.state(k).coh(1,:,:))) / K;
        M_psd = M_psd + squeeze(abs(fitmt_group_fact.state(k).psd(1,:,:))) / K;        
    end
end 

% if (total_band_num ~= (size(fitmt_group_fact.state(1).coh,1) - 1))
%     error('The number of bands and the fitmt group loaded does not match')
% end

for k = 1:1:K
        
  for band = 1:1:total_band_num 
   
        % WITHIN A BAND BUT FOR ALL STATES FOR THIS BAND
        if total_band_num > 1
            Band_tag = 'Fourband';
            M = zeros(ndim);
            M_psd = zeros(ndim);
            % Mean activation/coherence across all states
            if ~within_state
            for k1 = 1:K
                M = M + squeeze(abs(fitmt_group_fact.state(k1).coh(band,:,:))) / K;
                M_psd = M_psd + squeeze(abs(fitmt_group_fact.state(k1).psd(band,:,:))) / K;
            end
            end
            if within_state
            
            M = zeros(ndim);
            M_psd = zeros(ndim);
            % Within this state(k) taking mean across all the bands.
            % We will subtract this mean from the connections across all states
            % and plot the connections
            for band1 = 1:total_band_num
                M = M + squeeze(abs(fitmt_group_fact.state(k).coh(band1,:,:))) / total_band_num;
                M_psd = M_psd + squeeze(abs(fitmt_group_fact.state(k).psd(band1,:,:))) / total_band_num;
            end
            end
        end    
        
        coherence_value = squeeze(abs(fitmt_group_fact.state(k).coh(band,:,:))); 
        % PSD controls the size of the schemaball nodes
        % Subtracting the mean PSD calculated across states or within
        % states for the psd values for the current state and the band
        psd_amp_centered = diag(squeeze(fitmt_group_fact.state(k).psd(band,:,:))) - diag(M_psd);
        schemaball_size = 100*psd_amp_centered;
        [r_schm,c_schm] = find(schemaball_size <= 0);
        colr_mat_2_schm = colr_mat_2_orig;
        colr_mat_2_schm(r_schm,:) = repmat([0 0 0],length(find(schemaball_size <= 0)),1);
        schemaball_size(schemaball_size <= 0) = 30;
        schemaball_size(schemaball_size<10) = 30;
        
        % Coherence values for connectivity plot
        graph = coherence_value;        
        %Subtract the mean across states from the current state
        %coherence matrix
        graph = (graph - M);  
        if mask_connections
            mask = NaN(ndim,ndim);
            mask(lfp_pos(1):lfp_pos(2),cortex(1):cortex(2)) = ...
            ones(length(lfp_pos(1):lfp_pos(2)),...
            length(cortex(1):cortex(2)));
            mask(cortex(1):cortex(2),lfp_pos(1):lfp_pos(2)) = ...
            ones(length(cortex(1):cortex(2)),length(lfp_pos(1):lfp_pos(2))...
            );
            graph = graph.*mask;
        end
            tmp = squash(tril(graph));
            inds2 = find(tmp>1e-10);
            %             inds2 = find(abs(tmp) >1e-10);
            data = tmp(inds2);
        if ~exist('mask','var')
            mask = ones(ndim,ndim);
        end
        if (~isempty(data) && (length(find(mask == 1)) > 36 ) && threshold)
            
            S2 = [];
            S2.data = data;
            S2.do_fischer_xform = false;
            S2.do_plots = 0;
            S2.pvalue_th = P_VALUE/length(S2.data); %(Division corrects for multiple corrections)
            graph_ggm = teh_graph_gmm_fit(S2); 
            th = graph_ggm.normalised_th;
            graph = graph_ggm.data';
            graph(graph<th) = NaN;
            graphmat = zeros(ndim, ndim);
            graphmat(inds2) = graph;
            graph = graphmat;
            
        else
        if threshold
            graphmat = NaN(ndim, ndim);
            graph(graph < 0 | graph == 0) = 0;
            graph(21,22) = 0.05;
        end
            graph(graph < 0 | graph == 0) = 0;
        end
            graph(isnan(graph)) = 0;
            graph_schemaball = graph;
            graph_schemaball = normalize(graph_schemaball,1,'range',[0 1]);
            graph = tril(graph);
            % Rearranging graph for region based clustering
         if ndim == 48
            graph = graph([1:6,re_ROIs_48_index],[1:6,re_ROIs_48_index]);
            graph_schemaball = normalize(graph,1,'range',[0 1]);
%             C = figure((k*100 + band));
            set(gcf,'NumberTitle','off') %don't show the figure number
            set(gcf,'Name',['State ' num2str(k) ' band ' num2str(band)]) %select the name you want
%             circularGraph((graph),'Label',myLabel_2,'ColorMap',colr_mat_2);
            myLabel_2_schemaball = myLabel_2(:,1);
            
            %Numbering based labels
            for LBL = 1:1:48
                myLabel_2_schemaball{LBL,1} = num2str(LBL);
            end
            
%             savefig(C,[ Band_tag '_State ' num2str(k) ' band ' num2str(band)])
%             saveas(C,[Band_tag '_State ' num2str(k) ' band ' num2str(band)],'png');
%             close(C)
%             clear C
            C = figure((k*100 + band));
            set(gcf,'NumberTitle','off') %don't show the figure number
            set(gcf,'Name',['State ' num2str(k) ' band ' num2str(band)]) %select the name you want
            schemaball(graph_schemaball,myLabel_2_schemaball,[],colr_mat_2_schm,schemaball_size,colr_mat_2_orig,...
                ['State ' num2str(k) ' band ' num2str(band)]);
            clear graph graphmat th grap_ggm
%             close(C)
%             clear C
        elseif ndim == 42
            indices_gr = re_ROIs_48_index - 6;
            graph = graph(indices_gr,indices_gr);
            C = figure((k*100 + band));
            set(gcf,'NumberTitle','off') %don't show the figure number
            set(gcf,'Name',['State ' num2str(k) ' band ' num2str(band)]) %select the name you want
            circularGraph(graph,'Label',myLabel_3,'ColorMap',colr_mat_2(7:48,:));
            savefig(C,[Band_tag '_State ' num2str(k) ' band ' num2str(band)])
            saveas(C,[Band_tag '_State ' num2str(k) ' band ' num2str(band)],'png');
            clear graph graphmat th grap_ggm 
            close(C)
            clear C
        end
            

  end
end

%% TEMPORAL FEATURES

variable = Intervals;
for ii = 1:1:size(variable,2)
   
        state_int{1,ii} = cell2mat(cellfun(@transpose,variable(:,ii),'UniformOutput',false));
        
end

state_int = cellfun(@(x) x*0.001,state_int,'UniformOutput',false);
figure(100)
violin(state_int)

% box plots
y = num2cell(1:numel(state_int));
x = cellfun(@(x, y) [x(:) y*ones(size(x(:)))], state_int, y, 'UniformOutput', 0); % adding labels to the cells
X = vertcat(x{:});
boxplot(X(:,1), X(:,2))


%% Distance between different states of different HMMs

for c = 1:1:size(C_all_covs_ON,3)
    
    C1 = C_all_covs_ON(:,:,c);
    
    for c1 = 1:1:size(C_all_covs_OFF,3)
        C2 = C_all_covs_OFF(:,:,c1);
        E = eig(C1*C2);
          E1 = eig(C1,C2);
          logE = log(E);
          logE1 = log(E1);
          sqr = logE1 .* logE1;
          sumsqr = sum(sqr);
          sqrtsumsqr = sqrt(sumsqr);
          abssq = abs(sqrtsumsqr);
          dist1 = abssq;
%           dist = sqrt(sum(log
          C_all(c,c1) = dist1;
    end
end

h = heatmap(C_all,'Colormap',jet);
h.XLabel = 'ON1'; %columns
h.YLabel = 'ON2'; %rows ON

 % Get state matches
 
 % This results in columns corresponding to ON states and data entries as
 % their corresponding OFF
%  [assignment,cost] = munkres(C_all); 
 
 % This results in columns corresponding to OFF states and data entries as
 % their corresponding ON
 [assignment,cost] = munkres(C_all');