database ='D:\Abhinav_Sharma\RS_peri_MEGLFP\brainstorm_db\MEG_LFP_peri\';
ROIs = [3,4,5,6,11,12,15,16,19,20,21,22,23,24,25,26,27,28,29,30,...
    33,34,35,36,37,38,41,42,45,46,47,48,51,52,53,54,55,56,57,58,59,60];
destSurfFile = '@default_subject/tess_cortex_pial_low.mat';
surface_file = destSurfFile;
surface_file_data = load([database 'anat' filesep surface_file]);
Atlas_name = surface_file_data.Atlas(6).Name
Scouts = surface_file_data.Atlas(6).Scouts;

%% List files for the ROIs

directory = 'Q:\RS_MEG_LFP_peri\default_projected_data';
files = dir(directory);

for lsc = 3:1:length(files)
    
    file_name = files(lsc).name;
    file_name_num = strsplit(file_name,'_');
    file_name_num = strsplit(file_name_num{4},'.');
    ROI_num = str2num(file_name_num{1});
    
    
    custom_Scouts(ROI_num).Vertices = Scouts(ROI_num).Vertices;
    custom_Scouts(ROI_num).Label = Scouts(ROI_num).Label;
    
    load(file_name)
    
    sz = cellfun('size',DATA_scout,1);
    sz = num2cell(sz);
    [A_pc,B_pc,e_pc] = highdim_pca(DATA_scout,sz,1);
    
    custom_Scouts(ROI_num).spatialPC = A_pc;
    custom_Scouts(ROI_num).var_explained = e_pc(1);
    
    clear A_pc B_pc e_pc ROI_num sz file_name file_name_num

end

%% 

load('custom_Scouts')

Amp_file = zeros(15002,2);

    for lsc = 1:1:length(custom_Scouts)
     
     if ~isempty(custom_Scouts(lsc).spatialPC)
         
         AA = custom_Scouts(lsc).spatialPC - min(custom_Scouts(lsc).spatialPC);
         custom_Scouts(lsc).spatialPC = AA./(max(custom_Scouts(lsc).spatialPC)-min(custom_Scouts(lsc).spatialPC));
         Amp_file(custom_Scouts(lsc).Vertices) = custom_Scouts(lsc).spatialPC;
         Amp_file(custom_Scouts(lsc).Vertices,2) = r(lsc)*custom_Scouts(lsc).spatialPC; 
         
     end   
        
    end
    kernelMat.ImageGridAmp = Amp_file;
    kernelMat.Comment      = ['trial_new'];
    kernelMat.sRate      = 1;
    kernelMat.ImageGridTime         = 1:2;
    kernelMat.DataFile = [];
    kernelMat.Time         = 1:2;
    kernelMat.SurfaceFile = ['@default_subject\tess_cortex_pial_low.mat'];
    kernelMat.HeadModelType = 'surface';
    % Output filename
    OPTIONS.ResultFile = fullfile([database 'data\Group_analysis\'], ['results_Group_trial_new'] );
    % Save file
    save(OPTIONS.ResultFile, '-struct', 'kernelMat');
    