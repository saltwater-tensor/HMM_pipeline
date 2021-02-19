function rs_state_power_topo(timepoints, nstates, normalize, varargin)
%plots normalized spectral power for each state, specifically the spectra
%for power averaged over sources and the topograhpies for power averaged
%over frequencies within the delta, alpha and beta band, respectively
%inputs:
%timepoints: string; epoch for which states were estimated; can be 'baseline' or
%'spiral_slow'
%nstates: integer scalar; number of states

if nargin > 3
    freq_boundaries = varargin{1};
    band_labels = varargin{2};
    spec_measures = varargin{3};
    col_lim = varargin{4};
    filt = varargin{5};
else
    %freq band definition
    freq_boundaries = [0 3;4 7; 8 12; 13 35; 60 90];
    band_labels = {'delta','theta','alpha','beta','gamma'};
    spec_measures = {'pow','coh'};
    col_lim = [];
    filt = 'unfilt';
end

addpath(genpath('/net/citta/storage/home/jan/matlab_toolboxes/ScreenCapture'));

%subject info
load /data/project/megdbs/Jan/projects/rest_stim/Info/Info
subjects= fieldnames(Info);

%state spectra
if strcmp(filt,'filt')
load(['/data/user/jan/rest_stim/state_data2/state_spectra/',timepoints,'_',num2str(nstates),'.mat']);
elseif strcmp(filt,'unfilt')
load(['/data/user/jan/rest_stim/state_data2/state_spectra/',timepoints,'_',num2str(nstates),'_unfilt.mat']);
end

%parcellation
load /data/user/jan/rest_stim/headmodels/parcellation.mat

%template grid
load /data/user/jan/rest_stim/headmodels/template_grid_2052_sources.mat

%set to true if you want to relate to mean over all bl, spont_pre and
%spont_post states
rel_to_all_states = false;


%normalize by power averaged across states
freqs = state_spec.state(1).f;
pow_all_states = nan(size(state_spec.state(1).psd,2),length(freqs),nstates);
coh_all_states = nan(size(state_spec.state(1).psd,2),size(state_spec.state(1).psd,2),length(freqs),nstates);
for st = 1:nstates
    for f = 1:length(freqs)
        pow_this_freq= squeeze(state_spec.state(st).psd(f,:,:));
        coh_this_freq = squeeze(state_spec.state(st).coh(f,:,:));
        pow_all_states(:,f,st) = abs(diag(pow_this_freq));
        coh_all_states(:,:,f,st) = coh_this_freq;
    end
end
n = squeeze(nanmedian(pow_all_states,3));
n = repmat(n,[1 1 nstates]);

n2 = squeeze(nanmedian(coh_all_states,4));
n2 = repmat(n2,[1 1 1 nstates]);

if normalize
pow_all_states = (pow_all_states-n)./n;
coh_all_states = (coh_all_states-n2)./n2;
suffix = '_norm';
meas_y_labels = {'norm. state power [a.u.]','norm. state coherence [a.u.]'};
else
    suffix = '_unnorm';
    meas_y_labels = {'power [(T/cm)^2]','coherence [a.u.]'};
end

suffix = [suffix,'_',filt];

%replace pow and coh by spectra relative to mean over all bl, spont_pre and
%spont_post states
if rel_to_all_states
[pow_all_states,coh_all_states,freqs] = rs_state_spec_norm_over_all_states(timepoints);
suffix = '_grand_norm';
end

%average coherence over all sources
coh_all_states = squeeze(nanmean(coh_all_states,1));

%some settings for plots
pics(1).view = [0 90];
pics(1).brain = 'interp_b';
pics(1).light = [90 90];

pics(2).view = [90 0];
pics(2).brain = 'interp_b';
pics(2).light = [90 30];

pics(3).view = [-90 0];
pics(3).brain = 'interp_b';
pics(3).light = [-90 30];

pics(4).view = [0 0];
pics(4).brain = 'interp_b';
pics(4).light = [0 10];
% 
% pics(5).view =[90,0];
% pics(5).brain = 'interp_l';
% pics(5).light =[90,30];
% 
% pics(6).view =[-90,0];
% pics(6).brain = 'interp_r';
% pics(6).light =[-90,30];

pics_per_state_and_band = length(pics);
screen_width = 1684;
screen_height = 974;
pic_width = floor(screen_width/nstates);
pic_height = floor(screen_height/pics_per_state_and_band);
lefts = 1:pic_width:screen_width ;
heights = screen_height-pic_height:-pic_height:1;

col_lims.pow.spont_pre = [-0.1 1; -0.2 0.6; -0.1 0.2; -0.2 0.3];
col_lims.pow.baseline = [-0.35 0.3; -0.3 0.45; -0.3 0.6;-0.15 0.2];
col_lims.pow.spont_post= [-0.4 0.1;-0.2 0.4;0 1;-0.4 0.8];

for sp = 1:numel(spec_measures)
    meas = spec_measures{sp};
    meas_all_states = eval([meas,'_all_states']);
    
    figure_dir = ['/data/user/jan/rest_stim/figures2/for_pub/',meas,'/',timepoints,'_',num2str(nstates),'/'];
    if ~exist(figure_dir,'dir')
        mkdir(figure_dir)
    end
    
    colmax = max(max(max(meas_all_states)));
    colmin = min(min(min(meas_all_states)));
    collim = 0.8 * max([abs(colmax),abs(colmin)]);
    
    for fr = 1:size(freq_boundaries,1)

        bound_start = freq_boundaries(fr,1);
        start_ind = closest_index(bound_start,freqs);
        bound_end = freq_boundaries(fr,2);
        end_ind = closest_index(bound_end,freqs);
        
        %color lim common to states
        colmax = max(max(mean(meas_all_states(:,start_ind:end_ind,:),2)));
        colmin = min(min(mean(meas_all_states(:,start_ind:end_ind,:),2)));
        collim_common = max([abs(colmax),abs(colmin)]);
        
        for st = 1:nstates
            
            l = lefts(st);
            dummy_source = template_grid;
            dummy_source.avg.pow = nan(size(template_grid.pos,1),1);
            
            %average power within band
            mean_bandpow = squeeze(mean(meas_all_states(:,start_ind:end_ind,st),2));
            
            %write from parcellated source to original source space
            index_into_parcel = parcel.afni;
            has_parcel = parcel.afni~=0;
            dummy_source.avg.pow(has_parcel) = mean_bandpow(index_into_parcel(has_parcel));
            
            %get rid of sources that do not belong to any parcel
            del_ind = isnan(dummy_source.avg.pow);
            dummy_source.avg.pow(del_ind) = [];
            dummy_source.pos(del_ind,:) = [];
            
            %interpolate power on different brain surfaces
            cfg = [];
            cfg.parameter = 'pow';
            surf = ft_read_headshape('surface_pial_left.mat');
%             interp_l = ft_sourceinterpolate(cfg, dummy_source, surf);
%             
%             surf = ft_read_headshape('surface_pial_right.mat');
%             interp_r = ft_sourceinterpolate(cfg, dummy_source, surf);
            
            surf = ft_read_headshape('surface_pial_both.mat');
            interp_b = ft_sourceinterpolate(cfg, dummy_source, surf);
            
            %color lim for each state
            colmax = max(mean(meas_all_states(:,start_ind:end_ind,st),2));
            colmin = min(mean(meas_all_states(:,start_ind:end_ind,st),2));
            collim_state = max([abs(colmax),abs(colmin)]);
            
            %and plot
            cfg = [];
            cfg.method = 'surface';
            cfg.funcolormap = 'jet';
            cfg.funparameter = 'pow';
            cfg.anaparameter = 'anatomy';
            cfg.projmethod = 'project';
            cfg.colorbar = 'no';
            cfg.maskparameter = 'pow';
            if ~normalize
                cfg.funcolorlim = [colmin colmax];
            elseif ~isempty(col_lim)
                cfg.funcolorlim = col_lim;
            else
                cfg.funcolorlim = 'zeromax';
            end
            %cfg.funcolorlim = [-collim_state collim_state];
            
            for v = 1:length(pics)
                
                h = heights(v);
                
                ft_sourceplot(cfg,eval(pics(v).brain));
                set(gcf,'OuterPosition',[l h pic_width pic_height],'MenuBar','none','ToolBar','none','DockControls','off');
                view(pics(v).view)
                lightangle(pics(v).light(1),pics(v).light(2))
                if v==2
                    colorbar('Fontsize',18,'Location','WestOutside')
                end
                
            end
        end
        
        %save a screenshot
        drawnow;
        fig_name = [figure_dir,band_labels{fr},suffix,'.png'];
        pause(180)
        screencapture('handle',0,'position',[1 1 screen_width screen_height],'target',fig_name);
        close all
        
    end
    
    %plot power spectra
    meas_max= max(max(mean(meas_all_states,1)));
    meas_min = min(min(mean(meas_all_states,1)));
    f1 = closest_index(55,freqs);
    f2 = closest_index(65,freqs);
    
    figure;
    for st = 1:nstates
        subplot(1,nstates,st)
        state_powspec = squeeze(mean(meas_all_states(:,:,st),1));
        if ~any([isempty(f1),isempty(f2)])
        plot(freqs(1:f1),state_powspec(1:f1),'k')
        hold on
        plot(freqs(f2:end),state_powspec(f2:end),'k')
        hold off
        else
             plot(freqs,state_powspec,'k')
        end
        xlabel('frequency [Hz]');
        ylim([meas_min,meas_max]);
        if st==1
        ylabel(meas_y_labels{sp});
        end
        %axis('tight');
        set(gca,'FontSize',16);
        %title(['state ',num2str(st)])
        set(gcf,'Position',[79         520        1783         384]);
    end
    %print(gcf,'-r1200','-dpng',[figure_dir,meas,'_spectra',suffix,'_highres.png']);
    saveas(gcf,[figure_dir,meas,'_spectra',suffix,'._test.fig']);
    saveas(gcf,[figure_dir,meas,'_spectra',suffix,'._test.png']);
    
end
