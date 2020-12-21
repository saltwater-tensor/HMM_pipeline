load('C:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis\OFF\Whole_brain_stn_lfp_medication_OFF_06_Jan_2020_18_39_55_HMM_model_pca_NO_MAR_Motor_cortex_LFP_all_embed_lags\fitmt_subj_fact_4b_filtered1to45.mat')
OFF_subject_specific_spectral = fitmt_subj_fact_4b;
clear fitmt_subj_fact_4b
% ON Dataset
load('C:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis\ON\Whole_brain_stn_lfp_medication_ON_06_Jan_2020_18_45_18_HMM_model_pca_NO_MAR_Motor_cortex_LFP_all_embed_lags\fitmt_subj_fact_4b_filtered1to45_type1.mat')
ON_subject_specific_spectral = fitmt_subj_fact_4b;
clear fitmt_subj_fact_4b

%% OFF 
N = length(OFF_subject_specific_spectral); %Total number of subjects
K = length(OFF_subject_specific_spectral{1}.state); % Total number of states
[Nf,ndim,~] = size(OFF_subject_specific_spectral{1}.state(1).psd); 

% PSD
p = ndim;
X_off = zeros(N,K,p,Nf); 
% Coh
ind2 = triu(true(ndim),1);
p = ndim * (ndim-1) / 2; 

% X is collecting the coherence across all subjects in the OFF condition
% Subject x States x spatial_dimensions x No.of frequency modes
X_off = zeros(N,K,p,Nf); 
for n = 1:N
    for nf = 1:Nf
        for k = 1:K
            d = OFF_subject_specific_spectral{n}.state(k).coh(nf,ind2)';
            X_off(n,k,:,nf) = d(:);
        end
    end
end

%%
N = length(ON_subject_specific_spectral); %Total number of subjects
K = length(ON_subject_specific_spectral{1}.state); % Total number of states
[Nf,ndim,~] = size(ON_subject_specific_spectral{1}.state(1).psd); 

% PSD
p = ndim;
X_on = zeros(N,K,p,Nf); 
% Coh
ind2 = triu(true(ndim),1);
p = ndim * (ndim-1) / 2; 

% X is collecting the coherence across all subjects in the ON condition
% Subject x States x spatial_dimensions x No.of frequency modes
X_on = zeros(N,K,p,Nf); 
for n = 1:N
    for nf = 1:Nf
        for k = 1:K
            d = ON_subject_specific_spectral{n}.state(k).coh(nf,ind2)';
            X_on(n,k,:,nf) = d(:);
        end
    end
end


%% Permutation testing

% Delta band
factor = 3;
k_off = 1;
k_on = 2;

X_off_delta = squeeze(X_off(:,k_off,:,factor));
X_on_delta = squeeze(X_on(:,k_on,:,factor));

for conns = 1:1:size(X_off_delta,2) 
    
    conn_pair_off = X_off_delta(:,conns);
    conn_pair_on = X_on_delta(:,conns);
    
    [permutation_p(conns), observeddifference, effectsize] = permutationTest(conn_pair_off,conn_pair_on,...
          500,'sidedness','smaller' );
      
      [permutation_p_2(conns), observeddifference_2, effectsize_2] = permutationTest(conn_pair_off,conn_pair_on,...
          500,'sidedness','larger' );
    
end


