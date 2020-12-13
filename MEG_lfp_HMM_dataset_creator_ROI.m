%% ROI based exctraction procedure

% The atlas is specified in the section ATLAS BASED SIGNAL EXTRACTION
% This script exctracts a single principal component for each roi for each
% subject.

% Then concatenates the lfp data from each subject to the matrix 

% Finally the numerical value of the lfp data is adjusted to prevent bias

% No PCA is performed on the lfp data. 

%% Running brainstorm with no gui
% This is essential
brainstorm nogui

%% Local computer path names
%   database ='D:\Abhinav_Sharma\RS_peri_MEGLFP\brainstorm_db\MEG_LFP_peri\';
  database = 'Q:\RS_MEG_LFP_peri\brainstorm_db\RS_PD_peri\';

%% Subject type

%If you want to loop through all the subjects set autogenerate to 1
autogenerate = 1;

%Subject medication condition ON or OFF
med_state = input('Enter subject medication state ON or OFF');

%% PERFORM PCA 
% LFP only
lfp_only = 1;
% MEG only

% Whole data

%% STN only data

right_contacts = [2,4,7];
left_contacts =  [10,12,15];% 3 left and 3 right %14 was also tried

% THESE ROIS are for the MINDBOGGLE ATLAS
ROIs = [3,4,5,6,11,12,15,16,19,20,21,22,23,24,25,26,27,28,29,30,...
    33,34,35,36,37,38,41,42,45,46,47,48,51,52,53,54,55,56,57,58,59,60];

if ~exist('atlas_name_check','var')
    atlas_name_check = input('Enter atlas name for the above ROIs');
end


%% Atlas based processing
atlas_process = 1;

multi_site_dataset = 0;

%% Source space projection
use_source_space = 1;
% Project to default anatomy
project_to_default_anat = 1;
destSurfFile = '@default_subject/tess_cortex_pial_low.mat';


%% Clip data
%Timing considerations
% dwnsmpl = input('Do you want to downsample the signal?');
% if dwnsmpl
%     int_factor = input('By what factor? 2, 4 etc the new rate will originalrate/factor');
% end
% reducetime = input('Do you want to reduce signal time?');
% if reducetime
%     red_time = input('Enter reduced time in milliseconds 120,000 ms is 2 minutes');
% end

% In seconds
red_time_in_seconds = 300; %600 seconds is 10 minutes
reducetime = 0;


%% Autogenerating subject list
if autogenerate

parent_dir = fullfile(database,'data');
parent = dir(fullfile(database,'data'));
sub_num = 0;
run_once = 0;
subject = [];

for pr = 1:1:length(parent)
    
    parent(pr).name
    if (~isempty(regexpi(parent(pr).name,'S0','once'))) && (~isempty(regexpi(parent(pr).name,'peri','once')))
        subject_dir = dir(fullfile(parent_dir,filesep,parent(pr).name));
        sub_num = sub_num + 1;

        run_count = 1;
       
        for sbj = 1:length(subject_dir)
            
           if ((~isempty(regexpi(subject_dir(sbj).name,'rest','once')))...
                   && (isempty(regexpi(subject_dir(sbj).name,'@','once')))...
                   && (~isempty(regexp(subject_dir(sbj).name,med_state,'once'))))
                
               subject{sub_num,1} = parent(pr).name;
               subject{sub_num,2}{1,run_count} = subject_dir(sbj).name;
               run_count = run_count + 1;
               
               run_once = 1;     
           end
  
        end
        if size(subject,1) ~= sub_num && run_once
            sub_num = sub_num -1;
        end
    end
   
end

end


%% Data processing
total_instance = 1;
total_Data = {};
% left or right
side = input('L , R or B');

for isubject= 1:size(subject,1) %loop over subjects
    Bad_tmp=[];
    Data=[];
    total_data = [];
    total_data_cortical = [];
    subject_data{isubject,1} = subject{isubject,1};
    Data_cortex_lfp_total = [];
    combined_across_runs = [];
    
    for irun = 1:length(subject{isubject,2})
        lfp_save = 0;
        sRate=[];
        subject_region_signals_per_run = [];
        a=dir(fullfile([database 'data' filesep subject{isubject,1} filesep  subject{isubject,2}{irun}  filesep]));
        
        %% FINDING IMAGING KERNEL
        % As the analysis we will do is at the source level, the core MEG
        % cortex we will need the source MEG kernel
        ImagingKernel = [];
        
        for k=1:length(a)
            
            %% For MNE
%             if ~isempty(strfind(a(k).name, 'results_MN_MEG_GRAD_MEG_MAG_KERNEL_19'))
%                 Kernel = load ([database 'data' filesep subject{isubject,1} filesep subject{isubject,2}{irun} filesep a(k).name]);
%                 ImagingKernel = Kernel.ImagingKernel;
%                 surface_file = Kernel.SurfaceFile;
%                     
%             end
            
            %% For LCMV
            a(k).name
            if ~isempty(strfind(a(k).name, 'results_PNAI_MEG_GRAD_MEG_MAG_KERNEL_19'))
                Kernel = load ([database 'data' filesep subject{isubject,1} filesep subject{isubject,2}{irun} filesep a(k).name]);
                ImagingKernel = Kernel.ImagingKernel;
                surface_file = Kernel.SurfaceFile;
                break    
            end
            
        end
        
        %% FINDING RAW RECORDING IMPORTED AFTER CORRECTION INTO BRAINSTORM
        iCond = [];
        
        for i=1:length(a)
            if ~isempty(strfind(a(i).name, 'block'))
                iCond(end+1) = i;
            end
        end
        
        if size(iCond)>1 
            disp('There were more than one datafile')
            break
        end
        
        if isempty(iCond)
            disp('There was no datafile')
            break
        end
        
        data=  load ([database 'data' filesep  subject{isubject,1} filesep subject{isubject,2}{irun} filesep a(iCond(1)).name]);
        lfp_data = data.F(310:325,:);
        recordings = data.F(Kernel.GoodChannel(1:end),:);
        data.F = recordings;
        clear recordings
        
        if ~isempty(data)
            sRate = round(abs(1 / (data.Time(2) - data.Time(1)))); %new sampling rate This is the sampling rate that we obtain after BRAINSTORM preprocessing.
        end
        
%% MANUALLY REMOVING BAD DATA
      
      % All type of events have to be corrected together otherwise after
      % correction the time stamps would be changed
      event_type = 'BAD/BAD_LFP';
      Data = clean_bad(event_type,sRate,data);
      data.F = Data;
      recordings = data.F;
      lfp_Data = clean_bad(event_type,sRate,lfp_data,data);
      clear lfp_data
      clear Data

%% TIME REDUCTION
if reducetime

     red_time = sRate*red_time_in_seconds;
     time_for_run = size(recordings,2)./sRate;
     time_for_run = time_for_run/60;
      if time_for_run <= red_time_in_seconds/60 %Not processing trials which are not 10
%              mins long
          continue
      end

     time_ok = 1;
     if reducetime 
         if red_time > size(recordings,2)
             disp('Cannot reduce time, already below specified amount')
             time_ok = 0;
         else
             disp(['Reducing time to ' num2str(red_time_in_seconds/60) ' minutes '])
             recordings = recordings(:,1:red_time);
             data.Time = data.Time(1:red_time);
             lfp_Data = lfp_Data(:,1:red_time);
         end
     end
end


%%
     data.F = recordings;
     clear recordings
     
     
     
      %% RESAMPLING AND SOURCE LEVEL PROJECTION
      
      % Resampling
      % Forward model is unaffected by data and the inverse model relies on
      % data covariance for computation. Because the computation is linear
      % and treats every source as independent sampling frequency would
      % have a limited effect on the gain matrix or Kernel computed. Hence
      % we can resample at the sensor level and then project to the source
      % level.
      % Data = The new cleaned variable
      % lfp_Data = The new lfp signal
      % Resample both to 1000 Hz.
      
      Time = data.Time;
      downsample_rate = 250;
      
       if use_source_space
          [data.F, data.Time] = process_resample('Compute', data.F,data.Time,downsample_rate,'resample-cascade');
      end
      
      [lfp_Data, lfp_Data_time] = process_resample('Compute', lfp_Data,Time,downsample_rate,'resample-cascade');
      
      % Projecting MEG data to source level, lfp data remains as it is
%       clear Data
      if use_source_space
           try
              Data = ImagingKernel * data.F;
              ImagingKernel = [];
              if project_to_default_anat


                  % Calculate interpolation matrix
                  % Interpolation is calculated here
                  srcSurfFile = surface_file;
                  isInteractive = 0;
                  isStopWarped = [];
                  try
                  [Wmat, sSrcSubj, sDestSubj, srcSurfMat, destSurfMat, isStopWarped] = ...
                        tess_interp_tess2tess( srcSurfFile, destSurfFile, isInteractive, isStopWarped );
                  catch
                      error('BRAINSTORM SHOULD BE RUNNING IN THE BACKGROUND, DO NOT USE CLEAR ALL')
                  end
                  % Apply interpolation to the Data/ImageGridAmp matrix
                  % A local copy of the interpolation function exists in the
                  % hmm folder
                  disp('DEFAULT PROJECTION')
                  Data = muliplyInterp_local...
                    (Wmat, double(Data), Kernel.nComponents);



               end
           catch
               continue
           end
           
           
      else
          disp('Not using source space only using LFP data')
      end
      
      
      %% LFP Dataset
      if ~multi_site_dataset
      
      right_electrodes = lfp_Data(right_contacts,:);
      left_electrodes = lfp_Data(left_contacts,:);
      lfp_Data = [right_electrodes;left_electrodes];
      
      switch side
          case 'L'
              total_Data{total_instance,1} = left_electrodes';
              total_instance = total_instance + 1;
              
              total_data = vertcat(total_data,left_electrodes');
              T{1,isubject}(1,irun) = size(left_electrodes,2); 
              subject_data{isubject,2}{1,irun} = left_electrodes';
          case 'R'
              total_Data{total_instance,1} = right_electrodes';
              total_instance = total_instance + 1;
              
              total_data = vertcat(total_data,right_electrodes');
              T{1,isubject}(1,irun) = size(right_electrodes,2); 
              subject_data{isubject,2}{1,irun} = right_electrodes';
              
         case 'B'
              %Instance where you want to run the HMM on the combined left
              %and right electrode data
              total_Data{total_instance,1} = lfp_Data';
              total_instance = total_instance + 1;
              
              total_data = vertcat(total_data,lfp_Data');
              T{1,isubject}(1,irun) = size(lfp_Data,2); 
              subject_data{isubject,2}{1,irun} = lfp_Data';
              
      end
      
      end
     
      
      %% ATLAS BASED SIGNAL EXTRACTION
      if atlas_process && ~multi_site_dataset
          
          atlas_name = 'Mindboggle';
          if ~strcmpi(atlas_name,atlas_name_check)
              error('YOU HAVE NOT CHECKED THE ATLAS CHECK ROIs LINE 40 above too')
          end
          if project_to_default_anat
               surface_file = destSurfFile;
               surface_file_data = load([database 'anat' filesep surface_file]);
          end
          brk = 0;
          for atname = 1:1:length(surface_file_data.Atlas)
             if strcmpi(surface_file_data.Atlas(atname).Name,atlas_name)
                 brk = 1;
                 Scouts = surface_file_data.Atlas(atname).Scouts;
                 PCA_scout = [];
                 for sc = 1:1:length(ROIs) %length(Scouts)
                     
                     V = Scouts(ROIs(sc)).Vertices;
                     dt = Data(V,:);
                     
                     %%  PCA
    
                    pca_opts.stack_data='yes';
                    pca_opts.storage='full';
                    pca_opts.precision='double';
                    pca_opts.eig_solver='selective';
                    try
                        [V, Lambda] = icatb_calculate_pca(dt, min(25,size(dt,1)-1),...
                            'type', 'standard', 'whiten', 0, 'verbose', 1, ...
                            'preproc_type', 'remove mean per timepoint', ...
                        'pca_options', pca_opts, 'varToLoad',...
                        [], 'dims', size(dt));
                    catch
                        [V, Lambda] = icatb_calculate_pca(dt, 20,...
                            'type', 'standard', 'whiten', 0, 'verbose', 1, ...
                            'preproc_type', 'remove mean per timepoint', ...
                        'pca_options', pca_opts, 'varToLoad',...
                        [], 'dims', size(dt));
                        
                    end

                    
                    V = V';
                    V_last = V(end,:); 
                    PCA_scout = [PCA_scout;V_last];
                 end
                 
                 
             end 
             if brk==1
               break
             end
          end
          
          icasig = icatb_remove_mean(PCA_scout)';
          
          mean_PCA_scout = mean(PCA_scout,2);
          variance_PCA_scout = std(PCA_scout,[],2);
          
          mean_PCA_scout_1 = repmat(mean_PCA_scout,[1,size(PCA_scout,2)]);
          variance_PCA_scout_1 = repmat(variance_PCA_scout,[1,size(PCA_scout,2)]);
          
          PCA_scout_mean_sub = PCA_scout - mean_PCA_scout_1;
          PCA_scout_var_norm = PCA_scout_mean_sub ./ variance_PCA_scout_1;
          
          
          [L]= symmetric_orthogonalise_custom(PCA_scout_var_norm);
          
          L = L - repmat(mean(L,2),[1,size(L,2)]);
          L = L./repmat(std(L,[],2),[1,size(L,2)]);
           
          
           MEG_region_signals = L;
           combined_MEG_LFP_per_run = [lfp_Data',MEG_region_signals'];
           
           data_subject{1,isubject}{1,irun} = combined_MEG_LFP_per_run;
           T{1,isubject}(1,irun) = size(lfp_Data,2); 
           
           
           
           
           
      end 
     combined_across_runs = [combined_across_runs;combined_MEG_LFP_per_run];
      
      
    end
    DATA{isubject,1} = combined_across_runs;
    
    
end

%% Rescaling lFP DATA by subtracting mean and dividing by variance
% This rescales the data to the cortical PC scale so that HMM is not biased
% by numerical values of correlation between lfp signals
for d = 1:1:length(DATA)
    
    D = DATA{d,1};
    D_lfp = D(:,1:6);
    L = D_lfp';
    L = L - repmat(mean(L,2),[1,size(L,2)]);
    L = L./repmat(std(L,[],2),[1,size(L,2)]);
    
    D_lfp = L';
    D(:,1:6) = D_lfp;
    
    DATA{d,1} = D;
    
    
end












