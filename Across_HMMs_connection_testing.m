% Across HMMs all connections test
% We will test this OFF vs ON HMM only for the comparable states
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


%%
mask = ones(48,48);
% This is for coherence between all possible pairs of contacts 
coh_conns = triu(mask,1);
[coh_rows, coh_cols] = find(coh_conns == 1);
OFF_coh = [];
ON_coh = [];
coh_off_mean = NaN(6,3,48,48);
coh_on_mean = NaN(6,3,48,48);
%Since we are comparing only off and on for the three states under
%consideration we would store values for those three comparisons only
% 1 - HghPwr , 2 - Comms , 3 - Local
test_lesser = NaN(3,3,48,48);
test_greater = NaN(3,3,48,48);

for factor_number = 1:1:3
    for cns = 1:1:length(coh_rows)
        
        r = coh_rows(cns);
        c = coh_cols(cns);
        
        % Collect this connection across off and on conditions
        for offstate = 1:1:6
            for subject_number = 1:1:17
                
                data_off = squeeze(OFF_subject_specific_spectral...
                    {subject_number, 1}.state(offstate).coh(factor_number,r,c)) ; 
                data_off = data_off(:);
%                 data_off = data_off(find(~isnan(data_off)));
                OFF_coh = [data_off;OFF_coh];
            end
            %Store data for this state
             mean_off_coh(offstate) = mean(OFF_coh);
             err_off_coh(offstate) = std( OFF_coh ) / sqrt( length( OFF_coh ));
             All_off_coh(:,offstate) = OFF_coh;
             
             coh_off_mean(offstate,factor_number,r,c) = mean_off_coh(offstate);
             coh_off_raw{offstate,factor_number,r,c} = (OFF_coh);
             OFF_coh = [];
        end
        
        for onstate = 1:1:6
            for subject_number = 1:1:16
                
                data_on = squeeze(ON_subject_specific_spectral...
                    {subject_number, 1}.state(onstate).coh(factor_number,r,c)) ; 
                data_on = data_on(:);
%                 data_on = data_on(find(~isnan(data_on)));
                ON_coh = [data_on;ON_coh];

            end
            %Store data for this state
             mean_on_coh(onstate) = mean(ON_coh);
             err_on_coh(onstate) = std( ON_coh) / sqrt( length(ON_coh ));
             All_on_coh(:,onstate) = ON_coh;
             
             coh_on_mean(onstate,factor_number,r,c) = mean_off_coh(onstate);
             coh_on_raw{onstate,factor_number,r,c} = (ON_coh);
             ON_coh = [];
        end
      
      % Now compare this connection
      off_comp_hghpwr = coh_off_raw{1,factor_number,r,c};
      on_comp_hghpwr = coh_on_raw{2,factor_number,r,c};
      
      [h_off_lesser_hghpwr,p_off_lesser_hghpwr] = ttest2(off_comp_hghpwr,on_comp_hghpwr,'VarType','unequal','Tail','left'); %off < on
      [h_off_greater_hghpwr,p_off_greater_hghpwr] = ttest2(off_comp_hghpwr,on_comp_hghpwr,'VarType','unequal','Tail','right'); %off > on
      
%       [permutation_p, observeddifference, effectsize] = permutationTest(off_comp_hghpwr,on_comp_hghpwr,...
%           500,'sidedness','smaller' );
%       [permutation_p_2, observeddifference, effectsize] = permutationTest(off_comp_hghpwr,on_comp_hghpwr,...
%           500,'sidedness','larger' );
      
      test_lesser(1,factor_number,r,c) = p_off_lesser_hghpwr;
      test_greater(1,factor_number,r,c) = p_off_greater_hghpwr;
      
   
      off_comp_comms = coh_off_raw{2,factor_number,r,c};
      on_comp_comms = coh_on_raw{4,factor_number,r,c};
      
      [h_off_lesser_comms,p_off_lesser_comms] = ttest2(off_comp_comms,on_comp_comms,'VarType','unequal','Tail','left'); %off < on
      [h_off_greater_comms,p_off_greater_comms] = ttest2(off_comp_comms,on_comp_comms,'VarType','unequal','Tail','right'); %off > on
      
      
      
      test_lesser(2,factor_number,r,c) = p_off_lesser_comms;
      test_greater(2,factor_number,r,c) = p_off_greater_comms;
      
      
      
      off_comp_local = coh_off_raw{3,factor_number,r,c};
      on_comp_local = coh_on_raw{1,factor_number,r,c};
      
      [h_off_lesser_local,p_off_lesser_local] = ttest2(off_comp_local,on_comp_local,'VarType','unequal','Tail','left'); %off < on
      [h_off_greater_local,p_off_greater_local] = ttest2(off_comp_local,on_comp_local,'VarType','unequal','Tail','right'); %off > on
      
      test_lesser(3,factor_number,r,c) = p_off_lesser_local;
      test_greater(3,factor_number,r,c) = p_off_greater_local;
      
    end
end

%% Across specific regions clustered by anatomy
% IT IS PREFERRED TO RUN THIS PORTION OF THE ANALYSIS
contacts = [1,2,3,4,5,6];
frontal = [7,8,11,12,21,22,23,24,25,26,33,34,35,36];
medial_PFC = [1,2,15,16];
temporal = [17,18,39,40];
sensory_motor = [27,28,29,30];
parietal = [5,6,19,20,31,32,37,38,41,42];
visual = [3,4,9,10,13,14];

regions = {[1,2,3,4,5,6],[7,8,11,12,21,22,23,24,25,26,33,34,35,36],[1,2,15,16],[17,18,39,40],...
   [27,28,29,30],[5,6,19,20,31,32,37,38,41,42],[3,4,9,10,13,14] }';
mask = ones(48,48);
OFF_coh = [];
ON_coh = [];
coh_off_mean = NaN(6,3,7,7);
coh_on_mean = NaN(6,3,7,7);
%Since we are comparing only off and on for the three states under
%consideration we would store values for those three comparisons only
% 1 - HghPwr , 2 - Comms , 3 - Local
test_lesser = NaN(3,3,7,7);
test_greater = NaN(3,3,7,7);

for factor_number = 1:1:3
    
    for region1 = 1:1:7
        for region2 = 1:1:7
            
            mask = ones(48,48);
            %---------Select specific connections you want to run tests on
            rows = regions{region1,1};
            cols = regions{region2,1}; %1:6;
            mask(rows,cols) = 0;
            
            %Just for safety
            mask(cols,rows) = 0;
            mask = ~mask;
            %------------------
            
            mask = triu(mask,1);
            mask = double(mask);
            mask(find(mask == 0)) = NaN;
            
            for offstate = 1:1:6
                for subject_number = 1:1:17
                    
                    data_off = mask.* squeeze(OFF_subject_specific_spectral...
                        {subject_number, 1}.state(offstate).coh(factor_number,:,:)); 
                    data_off = data_off(:);
                    data_off = data_off(find(~isnan(data_off)));
                    OFF_coh = [data_off;OFF_coh];
                    
                end
                 mean_off_coh(offstate) = mean(OFF_coh);
                 err_off_coh(offstate) = std( OFF_coh ) / sqrt( length( OFF_coh ));
                 
                 
                 coh_off_raw_region{offstate,factor_number,region1,region2} = (OFF_coh);
                 
                 OFF_coh = [];
                 
            end
            
            
            for onstate = 1:1:6
                for subject_number = 1:1:16
                    
                    data_on = mask.*squeeze(ON_subject_specific_spectral...
                        {subject_number, 1}.state(onstate).coh(factor_number,:,:)); 
                    data_on = data_on(:);
                    data_on = data_on(find(~isnan(data_on)));
                    ON_coh = [data_on;ON_coh];
                    
                end
                 mean_on_coh(onstate) = mean(ON_coh);
                 err_on_coh(onstate) = std( ON_coh) / sqrt( length(ON_coh ));
                 
                 
                 coh_on_raw_region{onstate,factor_number,region1,region2} = (ON_coh);
                 
                 ON_coh = [];
            end
        
        % Now we have collected the the coherence values for this pair
        % of regions across all subjects
      off_comp_hghpwr = coh_off_raw_region{1,factor_number,region1,region2};
      on_comp_hghpwr = coh_on_raw_region{2,factor_number,region1,region2};
      
      [h_off_lesser_hghpwr,p_off_lesser_hghpwr] = ttest2(off_comp_hghpwr,on_comp_hghpwr,'VarType','unequal','Tail','left'); %off < on
      [h_off_greater_hghpwr,p_off_greater_hghpwr] = ttest2(off_comp_hghpwr,on_comp_hghpwr,'VarType','unequal','Tail','right'); %off > on
      
%       [permutation_p, observeddifference, effectsize] = permutationTest(off_comp_hghpwr,on_comp_hghpwr,...
%           500,'sidedness','smaller' );
%       [permutation_p_2, observeddifference, effectsize] = permutationTest(off_comp_hghpwr,on_comp_hghpwr,...
%           500,'sidedness','larger' );
      
      test_lesser(1,factor_number,region1,region2) = p_off_lesser_hghpwr;
      test_greater(1,factor_number,region1,region2) = p_off_greater_hghpwr;
      
      off_comp_comms = coh_off_raw_region{2,factor_number,region1,region2};
      on_comp_comms = coh_on_raw_region{4,factor_number,region1,region2};
      
      [h_off_lesser_comms,p_off_lesser_comms] = ttest2(off_comp_comms,on_comp_comms,'VarType','unequal','Tail','left'); %off < on
      [h_off_greater_comms,p_off_greater_comms] = ttest2(off_comp_comms,on_comp_comms,'VarType','unequal','Tail','right'); %off > on
      
      test_lesser(2,factor_number,region1,region2) = p_off_lesser_comms;
      test_greater(2,factor_number,region1,region2) = p_off_greater_comms;
      
      off_comp_local = coh_off_raw_region{3,factor_number,region1,region2};
      on_comp_local = coh_on_raw_region{1,factor_number,region1,region2};
      
      [h_off_lesser_local,p_off_lesser_local] = ttest2(off_comp_local,on_comp_local,'VarType','unequal','Tail','left'); %off < on
      [h_off_greater_local,p_off_greater_local] = ttest2(off_comp_local,on_comp_local,'VarType','unequal','Tail','right'); %off > on
      
      test_lesser(3,factor_number,region1,region2) = p_off_lesser_local;
      test_greater(3,factor_number,region1,region2) = p_off_greater_local;
      
        end
    end


end







%% Figures 

% All greater tests are plotted


p_thresh = 0.01;
% Beta band
F = figure(1);
D = squeeze(test_greater(1,3,:,:));
D = reshape(mafdr(D(:),'BHFDR',true),[7,7]);
%D(find(D > p_thresh)) = NaN;
h = heatmap(D);
h.XLabel = 'Regions'; %columns
h.YLabel = 'Regions'; %rows
title_str = 'High Power beta band p OFFgreaterON';
h.Title = 'High Power beta band p OFF > ON';
h.XDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
h.YDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
saveas(F,title_str);
saveas(F,title_str,'png');
close(F)
clear F



F =figure(2)
D = squeeze(test_greater(2,3,:,:));
D = reshape(mafdr(D(:),'BHFDR',true),[7,7]);
%D(find(D > p_thresh)) = NaN;
h = heatmap(D);
h.XLabel = 'Regions'; %columns
h.YLabel = 'Regions'; %rows
title_str ='Comms beta band p OFFgreaterON';
h.Title = 'Comms beta band p OFF > ON';
h.XDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
h.YDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
saveas(F,title_str);
saveas(F,title_str,'png');
close(F)
clear F



F =figure(3)
D = squeeze(test_greater(3,3,:,:));
D = reshape(mafdr(D(:),'BHFDR',true),[7,7]);
%D(find(D > p_thresh)) = NaN;
h = heatmap(D);
h.XLabel = 'Regions'; %columns
h.YLabel = 'Regions'; %rows
title_str ='Local beta band p OFFgreaterON';
h.Title = 'Local beta band p OFF > ON';
h.XDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
h.YDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
saveas(F,title_str);
saveas(F,title_str,'png');
close(F)
clear F



% Alpha Band
F =figure(4)
D = squeeze(test_greater(1,2,:,:));
D = reshape(mafdr(D(:),'BHFDR',true),[7,7]);
%D(find(D > p_thresh)) = NaN;
h = heatmap(D);
h.XLabel = 'Regions'; %columns
h.YLabel = 'Regions'; %rows
title_str ='High Power alpha band p OFFgreaterON';
h.Title = 'High Power alpha band p OFF > ON';
h.XDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
h.YDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
saveas(F,title_str);
saveas(F,title_str,'png');
close(F)
clear F



F =figure(5)
D = squeeze(test_greater(2,2,:,:));
D = reshape(mafdr(D(:),'BHFDR',true),[7,7]);
%D(find(D > p_thresh)) = NaN;
h = heatmap(D);
h.XLabel = 'Regions'; %columns
h.YLabel = 'Regions'; %rows
title_str = 'Comms alpha band p OFFgreaterON';
h.Title = 'Comms alpha band p OFF > ON';
h.XDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
h.YDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
saveas(F,title_str);
saveas(F,title_str,'png');
close(F)
clear F



F =figure(6)
D = squeeze(test_greater(3,2,:,:));
D = reshape(mafdr(D(:),'BHFDR',true),[7,7]);
%D(find(D > p_thresh)) = NaN;
h = heatmap(D);
h.XLabel = 'Regions'; %columns
h.YLabel = 'Regions'; %rows
title_str = 'Local alpha band p OFFgreaterON';
h.Title = 'Local alpha band p OFF > ON';
h.XDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
h.YDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
saveas(F,title_str);
saveas(F,title_str,'png');
close(F)
clear F



% Delta/theta Band
F =figure(7)
D = squeeze(test_greater(1,1,:,:));
D = reshape(mafdr(D(:),'BHFDR',true),[7,7]);
%D(find(D > p_thresh)) = NaN;
h = heatmap(D);
h.XLabel = 'Regions'; %columns
h.YLabel = 'Regions'; %rows
title_str = 'High Power deltatheta band p OFFgreaterON';
h.Title = 'High Power delta/theta band p OFF > ON';
h.XDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
h.YDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
saveas(F,title_str);
saveas(F,title_str,'png');
close(F)
clear F



F =figure(8)
D = squeeze(test_greater(2,1,:,:));
D = reshape(mafdr(D(:),'BHFDR',true),[7,7]);
%D(find(D > p_thresh)) = NaN;
h = heatmap(D);
h.XLabel = 'Regions'; %columns
h.YLabel = 'Regions'; %rows
title_str = 'Comms deltatheta band p OFFgreaterON';
h.Title = 'Comms deltatheta band p OFF > ON';
h.XDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
h.YDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
saveas(F,title_str);
saveas(F,title_str,'png');
close(F)
clear F



F =figure(9)
D = squeeze(test_greater(3,1,:,:));
D = reshape(mafdr(D(:),'BHFDR',true),[7,7]);
%SD(find(D > p_thresh)) = NaN;
h = heatmap(D);
h.XLabel = 'Regions'; %columns
h.YLabel = 'Regions'; %rows
title_str = 'Local deltatheta band p OFFgreaterON';
h.Title = 'Local delta/theta band p OFF > ON';
h.XDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
h.YDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
saveas(F,title_str);
saveas(F,title_str,'png');
close(F)
clear F


%% Figures 

% All lesser tests are plotted 


p_thresh = 0.01;
% Beta band
F = figure(1);
D = squeeze(test_lesser(1,3,:,:));
D = reshape(mafdr(D(:),'BHFDR',true),[7,7]);
%D(find(D > p_thresh)) = NaN;
h = heatmap(D);
h.XLabel = 'Regions'; %columns
h.YLabel = 'Regions'; %rows
title_str = 'High Power beta band p OFFlesserON';
h.Title = 'High Power beta band p OFF > ON';
h.XDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
h.YDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
saveas(F,title_str);
saveas(F,title_str,'png');
close(F)
clear F



F =figure(2)
D = squeeze(test_lesser(2,3,:,:));
D = reshape(mafdr(D(:),'BHFDR',true),[7,7]);
%D(find(D > p_thresh)) = NaN;
h = heatmap(D);
h.XLabel = 'Regions'; %columns
h.YLabel = 'Regions'; %rows
title_str ='Comms beta band p OFFlesserON';
h.Title = 'Comms beta band p OFF > ON';
h.XDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
h.YDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
saveas(F,title_str);
saveas(F,title_str,'png');
close(F)
clear F



F =figure(3)
D = squeeze(test_lesser(3,3,:,:));
D = reshape(mafdr(D(:),'BHFDR',true),[7,7]);
%D(find(D > p_thresh)) = NaN;
h = heatmap(D);
h.XLabel = 'Regions'; %columns
h.YLabel = 'Regions'; %rows
title_str ='Local beta band p OFFlesserON';
h.Title = 'Local beta band p OFF > ON';
h.XDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
h.YDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
saveas(F,title_str);
saveas(F,title_str,'png');
close(F)
clear F



% Alpha Band
F =figure(4)
D = squeeze(test_lesser(1,2,:,:));
D = reshape(mafdr(D(:),'BHFDR',true),[7,7]);
%D(find(D > p_thresh)) = NaN;
h = heatmap(D);
h.XLabel = 'Regions'; %columns
h.YLabel = 'Regions'; %rows
title_str ='High Power alpha band p OFFlesserON';
h.Title = 'High Power alpha band p OFF > ON';
h.XDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
h.YDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
saveas(F,title_str);
saveas(F,title_str,'png');
close(F)
clear F



F =figure(5)
D = squeeze(test_lesser(2,2,:,:));
D = reshape(mafdr(D(:),'BHFDR',true),[7,7]);
%D(find(D > p_thresh)) = NaN;
h = heatmap(D);
h.XLabel = 'Regions'; %columns
h.YLabel = 'Regions'; %rows
title_str = 'Comms alpha band p OFFlesserON';
h.Title = 'Comms alpha band p OFF > ON';
h.XDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
h.YDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
saveas(F,title_str);
saveas(F,title_str,'png');
close(F)
clear F



F =figure(6)
D = squeeze(test_lesser(3,2,:,:));
D = reshape(mafdr(D(:),'BHFDR',true),[7,7]);
%D(find(D > p_thresh)) = NaN;
h = heatmap(D);
h.XLabel = 'Regions'; %columns
h.YLabel = 'Regions'; %rows
title_str = 'Local alpha band p OFFlesserON';
h.Title = 'Local alpha band p OFF > ON';
h.XDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
h.YDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
saveas(F,title_str);
saveas(F,title_str,'png');
close(F)
clear F



% Delta/theta Band
F =figure(7)
D = squeeze(test_lesser(1,1,:,:));
D = reshape(mafdr(D(:),'BHFDR',true),[7,7]);
%D(find(D > p_thresh)) = NaN;
h = heatmap(D);
h.XLabel = 'Regions'; %columns
h.YLabel = 'Regions'; %rows
title_str = 'High Power deltatheta band p OFFlesserON';
h.Title = 'High Power delta/theta band p OFF > ON';
h.XDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
h.YDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
saveas(F,title_str);
saveas(F,title_str,'png');
close(F)
clear F



F =figure(8)
D = squeeze(test_lesser(2,1,:,:));
D = reshape(mafdr(D(:),'BHFDR',true),[7,7]);
%D(find(D > p_thresh)) = NaN;
h = heatmap(D);
h.XLabel = 'Regions'; %columns
h.YLabel = 'Regions'; %rows
title_str = 'Comms deltatheta band p OFFlesserON';
h.Title = 'Comms deltatheta band p OFF > ON';
h.XDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
h.YDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
saveas(F,title_str);
saveas(F,title_str,'png');
close(F)
clear F



F =figure(9)
D = squeeze(test_lesser(3,1,:,:));
D = reshape(mafdr(D(:),'BHFDR',true),[7,7]);
%SD(find(D > p_thresh)) = NaN;
h = heatmap(D);
h.XLabel = 'Regions'; %columns
h.YLabel = 'Regions'; %rows
title_str = 'Local deltatheta band p OFFlesserON';
h.Title = 'Local delta/theta band p OFF > ON';
h.XDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
h.YDisplayLabels = {'contacts','frontal','medial PFC','temporal','sensory motor','parietal','visual'};
saveas(F,title_str);
saveas(F,title_str,'png');
close(F)
clear F

