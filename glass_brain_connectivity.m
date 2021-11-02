database ='D:\Abhinav_Sharma\RS_peri_MEGLFP\brainstorm_db\MEG_LFP_peri\';
ROIs_new = [3,4,5,6,11,12,15,16,19,20,21,22,23,24,25,26,27,28,29,30,...
    33,34,35,36,37,38,41,42,45,46,47,48,51,52,53,54,55,56,57,58,59,60];
destSurfFile = '@default_subject/tess_cortex_pial_low.mat';
surface_file = destSurfFile;
surface_file_data = load([database 'anat' filesep surface_file]);
Atlas_name = surface_file_data.Atlas(6).Name
Scouts = surface_file_data.Atlas(6).Scouts;
load('D:\Abhinav_Sharma\RS_peri_MEGLFP\brainstorm_db\MEG_LFP_peri\anat\@default_subject\tess_cortex_pial_low.mat')

a= surface_file_data;
load('Q:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis\OFF\Whole_brain_stn_lfp_medication_OFF_06_Jan_2020_18_39_55_HMM_model_pca_NO_MAR_Motor_cortex_LFP_all_embed_lags\supply_profiles_4b_type4_filtered1to45\keyfunc_dot.mat')
% close all
%% Load a connectivity matrix
% clear ROIs
% load('re_ROIs_2')
r = [];
c = [];
adj_mat = [];
X = [];
Y = [];
Z = [];
connsx = [];
connsy = [];
connsz = [];
unique_label_points = [];
text_coordinate = [];
coordinate_names = {};

ROIs_new = re_ROIs_2;
% ROIs_new = re_ROIs_48_index_complete;
% adj_mat = graph_schemaball(4:45,4:45);
adj_mat = graph_schemaball;

[r,c] = find(adj_mat > 0 );

if ~isempty(r)

for r1 = 1:1:length(r)
    
    point1 = r(r1);
    ROI_actual_1 = ROIs_new(point1);
    ROI_actual_1_seed_vertex = Scouts(ROI_actual_1).Seed;
    point2 = c(r1);
    ROI_actual_2 = ROIs_new(point2);
    ROI_actual_2_seed_vertex = Scouts(ROI_actual_2).Seed;
    
    X = [a.Vertices(ROI_actual_1_seed_vertex,1);...
        a.Vertices(ROI_actual_2_seed_vertex,1)];
    
    Y = [a.Vertices(ROI_actual_1_seed_vertex,2);...
        a.Vertices(ROI_actual_2_seed_vertex,2)];
    
    Z = [a.Vertices(ROI_actual_1_seed_vertex,3);...
        a.Vertices(ROI_actual_2_seed_vertex,3)];
    
    connsx{r1,1} = X;
    connsy{r1,1} = Y;
    connsz{r1,1} = Z;
    
end
unique_label_points = [r;c];
unique_label_points = unique(unique_label_points);

for u = 1:1:length(unique_label_points)
    
    point = unique_label_points(u);
    ROI_actual = ROIs_new(point);
    ROI_actual_seed_vertex = Scouts(ROI_actual).Seed;
    ROI_actual_seed_label = Scouts(ROI_actual).Label;
    
    text_coordinate(u,1) = a.Vertices(ROI_actual_seed_vertex,1);
    text_coordinate(u,2) = a.Vertices(ROI_actual_seed_vertex,2);
    text_coordinate(u,3) = a.Vertices(ROI_actual_seed_vertex,3);
    coordinate_names{u,1} = ROI_actual_seed_label;
    
end

end

%% open brainstorm
brainstorm_open =0 ;
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
if length(FigList > 0)
  for iFig = 1:length(FigList)
      FigHandle = FigList(iFig);
      FigName   = get(FigHandle, 'Name');
  if (strcmpi(FigName,'Brainstorm'))
      brainstorm_open = 1;
  end
  end
    
end 
if ~brainstorm_open 
     brainstorm
end
clear FigList

% open figure on top of which you want to plot the connections
% ensure that the 3d brain is selected
if ~isempty(r)
    
    view_surface_data([], 'Group_analysis\Multilevel_group_ICA\results_Group_trial.mat');
      
      view_surface_matrix(Vertices, Faces)
    waitfor(msgbox('Ensure that the 3d brain is open and selected click on the figure bar and then ok'));
    gca;
    hold on
    for rr = 1:1:length(connsx)

         X = connsx{rr};
         Y = connsy{rr};
         Z = connsz{rr};
         plot3(X,Y,Z,'Marker','o','MarkerSize',12,'MarkerFaceColor','r',...
             'LineWidth',5,'Color','y','MarkerEdgeColor','none')
    end
clear X Y Z connsx connsy connsz

    Display text
    for ut = 1:1:length(coordinate_names)

        text(text_coordinate(ut,1),text_coordinate(ut,2),text_coordinate(ut,3),...
            coordinate_names{ut},'Color','white','FontSize',10,'FontWeight','normal')


    end
    waitfor(msgbox('Set transparency to 60%'));
    figure_3d('FigureKeyPressedCallback',   gcf, keyEvent)
    
    figure_3d_save(['State ' num2str(k)],[' band ' num2str(band) ' '])
waitfor(msgbox('Save all figures'));
close all
end

%% 
function [] = figure_3d_save( state,band )
% FolderName = 'C:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis\ON\Whole_brain_stn_lfp_medication_ON_06_Jan_2020_18_45_18_HMM_model_pca_NO_MAR_Motor_cortex_LFP_all_embed_lags\supply_profiles_4b_filtered1to45_type1A\Glassbrain\';
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
FolderName = 'Q:\RS_MEG_LFP_peri\Mindboggle_analysis_for_Elife\Panel_B_0.05\Final_figures';
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');
  FigName = strrep(FigName,'/','_');
  FigName = strrep(FigName,': ','_');
%   FigHandle.Renderer = 'painters';
  if (~strcmpi(FigName,'Brainstorm'))
    FigName = [state band FigName num2str(iFig)];
    saveas(FigHandle, [FolderName FigName '.svg']);
  end
end

end