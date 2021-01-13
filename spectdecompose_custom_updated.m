function [X,sp_fit,sp_fit_group,sp_profiles] = spectdecompose_custom_updated(sp_fit,options,chan,chan1,chan2,supply_profiles)


%% INPUTS
% chan             Vector specifying channels.If supplied a chan X chan submatrix is selected from coherence or
%                  the psd channels accordingly are selected. This is the symmetric version
%                  used for factorisation.
% chan1            Vector specifying channels.This option has to be used with chan2. This is used for assymetric selection
%                  of a subamtrix where chan1 selects the row numbers and chan2 selects the
%                  columns.
% chan2            Vector specifying channels.Selectrs matrix columns, see chan1 for explanation.
% supply_profiles  The decomposed profiles obtained in some previous
%                  attempt or from anywhere else to which you would like to project your
%                  data.
% options          options supplied during spectral decomposition using
%                  multitaper approach

%% OUTPUTS
% X               If you only need back the entire data matrix on the group level
%                 assimilated for factorisation then you can just get back X.

%% CHECKS AND STRUCTURE CONFIGURATION

% Number of input and output check
nargoutchk(1,4)
narginchk(0,7)


if ~iscell(sp_fit)
    error('Variable fit needs to be a cell, with one estimation per subject')
end

if nargin < 2, options = struct(); end
if ~isfield(options,'Niterations'), options.Niterations = 10; end
if ~isfield(options,'Ncomp'), options.Ncomp = 4; end
if ~isfield(options,'Method'), options.Method = 'NNMF'; end
if ~isfield(options,'Base'), options.Base = 'coh'; 
else
    if strcmpi(options.Base,'psd') 
        error('Current version of this program cannot function with psd as base') 
    end
end
if ~isfield(options,'plot'), options.plot = 0; end

% Total number of channels
ndim = size(sp_fit{1}.state(1).psd,2); 

if ~isempty(chan) && ~isempty(chan1)
    error('Both symmetric and assymetric channel/submatrix selection is not possible, either make chan as [] or chan1 and chan2 as []')
end
if isempty(chan), chan = 1:ndim; 
else
    ndim = length(chan);
end

if ~isempty(chan1)
    ndim1 = length(chan1); 
    ndim2 = length(chan2);
    [p,q] = meshgrid(chan1,chan2);
    pairs = [p(:) q(:)];
    ndim_3 = size(pairs,1);
    assym_flag = 1;
else
    chan1 = chan;
    chan2 = chan;
    ndim1 = length(chan1); 
    ndim2 = length(chan2);
    assym_flag = 0;
end



Nf = size(sp_fit{1}.state(1).psd,1); % no. frequencies
N = length(sp_fit); % no. subjects
K = length(sp_fit{1}.state); % no. states
ind_offdiag = triu(true(ndim),1)==1;
freq_bins = length(sp_fit{1,1}.state(1).f);
ndim4 = ndim*(ndim-1)/2; 
if ~isempty(supply_profiles)
    assert(size(supply_profiles,1) == freq_bins)
try
    assert(size(supply_profiles,2) == options.Ncomp)
catch
    options.Ncomp = size(supply_profiles,2);
end
end

% check that all estimations have the same number of freq bins
for n = 1:N
    if isempty(sp_fit{n})
        error('One or more elements of sp_fit are empty - remove them first')
    end
    if size(sp_fit{n}.state(1).psd,1) ~= Nf
        error(['It is necessary for the spectral estimation of all subjects ' ...
            'to have the same number of frequency bins. ' ...
            'In the case of the multitaper, this is can be done by setting ' ...
            'the options tapers and win to a fixed value.'])
    end
end

%% SELECTION OF RELEVANT OPERATIONS

% Compute factors or not
if ~isempty(supply_profiles)
    stop_factorisation = 1;
    disp('Factors supplied,  NNMF or PCA factorisation will NOT take place')
else
     disp('Factors not supplied, NNMF or PCA factorisation will take place')
end

% Projections
project_group_level = input('Do you want to project at the group level');
project_subject_level = input('Do you want to project at the subject level');


%% GROUP LEVEL MATRIX CONSTRUCTION
% put coh and psd in temporary arrays
coh_comps = zeros(N,K,Nf,ndim1,ndim2);

% for psd we will always keep all the channels
psd_comps = zeros(N,K,Nf,ndim);

% accumulating psd coherence values state wise and subject wise
for n = 1:N
    for k = 1:K
        % For a given subject select the proper channels
        psd = sp_fit{n}.state(k).psd(:,chan,chan);
        coh = sp_fit{n}.state(k).coh(:,chan1,chan2);
        % The selected indexes will now begin from 1
        
        for j = 1:1:length(chan)
            psd_comps(n,k,:,j) = psd(:,j,j);
%             for l = 1:1:length(chan2)
%                 coh_comps(n,k,:,j,l) = coh(:,j,l);
%             end
        end
        
        % The selected indexes will now begin from 1
        for j = 1:1:length(chan1)
            
            for l = 1:1:length(chan2)
                coh_comps(n,k,:,j,l) = coh(:,j,l);
            end
        end
        
    end
end

% Build the matrix that is to be factorised

% PSD
Xpsd = zeros(Nf,ndim*K);
for k=1:K
    ind = (1:ndim) + (k-1)*ndim;
    Xpsd(:,ind)= squeeze(mean(abs(psd_comps(:,k,:,:)),1));
end

% COH
if assym_flag
    Xcoh = zeros(Nf,K*ndim_3);
    for k = 1:K
        ind = (1:ndim_3) + (k-1)*ndim_3;
        % Taking mean across subjects
        ck = squeeze(mean(abs(coh_comps(:,k,:,:,:)),1));
    %     Xcoh(:,ind) = ck(:,ind_offdiag);
        ck = reshape(ck,[freq_bins,ndim1*ndim2]);
        Xcoh(:,ind) = ck;
    end
    else
        % This calculates total number of pairs
        % I know we could have replaced this with ndim_3 but for safety sake we
        % will keep the original part of the code
        ndim4 = ndim*(ndim-1)/2; 
        Xcoh = zeros(Nf,K*ndim4);
    for k = 1:K
        ind = (1:ndim4) + (k-1)*ndim4;
        ck = squeeze(mean(abs(coh_comps(:,k,:,:,:)),1));
        Xcoh(:,ind) = ck(:,ind_offdiag);
    end
end
% choose which matrix to factorise
if strcmpi(options.Base,'psd')
    X = Xpsd;
else
    X = Xcoh;
end
% Returns back X and ends the function
if nargout == 1
    return
end
%% GROUP LEVEL FACTORISATION 
if ~stop_factorisation
    
%% ROBUST NNMF FITTING
% This Gaussian function and fit has nothing to do with the description of
% mixture of two gaussians fit to find the connections to be visualised
% This gaussian is used to measure the goodness of fit of the found NNMF
% solution
% Specify fit function, a unimodal gaussian
gauss_func = @(x,f) f.a1.*exp(-((x-f.b1)/f.c1).^2);
% Default fit options
options_fit = fitoptions('gauss1');
% constrain lower and upper bounds
options_fit.Lower = [0,1,0];
options_fit.Upper = [Inf,size(X,1),size(X,1)];


%% 
if ~isfield(options,'sp_profiles') || isempty(options.sp_profiles)
    % Doing the decomposition
    if strcmpi(options.Method,'NNMF')
        bestval = Inf; fitval = zeros(options.Niterations,1); 
        for ii = 1:options.Niterations
            try
                [A,~] = nnmf(X,options.Ncomp,'replicates',500,'algorithm','als');
            catch
                error('nnmf not found - perhaps the Matlab version is too old')
            end 
            % Iterate over the four modes found above 
            for jj = 1:size(A,2)
                f = fit( linspace(1,size(A,1),size(A,1))',A(:,jj), 'gauss1',options_fit);
                residuals = A(:,jj) - gauss_func(1:size(A,1),f)';
                fitval(ii) = fitval(ii) + sum( residuals.^2 ) / size(A,2);
            end
            if fitval(ii) < bestval
                bestval = fitval(ii);
                bestfit = A; 
            end
        end
        %sp_profiles = pinv(X') * b'; % you lose non-negativity by doing this
        sp_profiles = bestfit;
        [~,ind] = max(sp_profiles,[],1);
        [~,neworder_auto] = sort(ind);
        sp_profiles = sp_profiles(:,neworder_auto);
    else
        try
            m = mean(X);
            X = X - repmat(m,Nf,1);
            [~,A] = pca(X,'NumComponents',options.Ncomp);
        catch
            error('Error running pca - maybe not matlab''s own?')
        end
        sp_profiles = A;
    end
else
    sp_profiles = options.sp_profiles;
end

% plot if required
if options.plot
    figure;
    if options.Ncomp > 4
        j1 = ceil(options.Ncomp/2); j2 = 2;
    else
        j1 = options.Ncomp; j2 = 1;
    end
    for j = 1:options.Ncomp
        subplot(j1,j2,j)
        plot(sp_profiles(:,j),'LineWidth',2.5)
    end
end


else
   
    sp_profiles = supply_profiles;
    
end


%% GROUP LEVEL PROJECTION

% Factorisation takes place for the given channels 
% But projection will happen for all the data
if project_group_level
    
    sp_fit_group = struct();
    sp_fit_group.state = struct();
    
    %% Reconstructing matrices for whole data projection
    clear coh_comps psd_comps Xpsd Xcoh
    coh_comps = zeros(N,K,Nf,ndim,ndim);
    psd_comps = zeros(N,K,Nf,ndim);
for n = 1:N
    for k = 1:K
        % For a given subject select the proper channels
        psd = sp_fit{n}.state(k).psd(:,chan,chan);
        coh = sp_fit{n}.state(k).coh(:,chan,chan);
        % The selected indexes will now begin from 1
        % PSD is self it is not cross spectral across 
        for j = 1:1:length(chan)
            psd_comps(n,k,:,j) = psd(:,j,j);
            for l = 1:1:length(chan)
                coh_comps(n,k,:,j,l) = coh(:,j,l);
            end
        end
    end
end

% Build the matrix that is to be factorised
Xpsd = zeros(Nf,ndim*K);
for k=1:K
    ind = (1:ndim) + (k-1)*ndim;
    Xpsd(:,ind)= squeeze(mean(abs(psd_comps(:,k,:,:)),1));
end
Xcoh = zeros(Nf,K*ndim4);
for k = 1:K
    ind = (1:ndim4) + (k-1)*ndim4;
    ck = squeeze(mean(abs(coh_comps(:,k,:,:,:)),1));
    Xcoh(:,ind) = ck(:,ind_offdiag);
end
    
    
    
    
    
    %% 
    if strcmpi(options.Method,'NNMF')
       psd = (Xpsd' * sp_profiles); % ndim by components
       coh = (Xcoh' * sp_profiles);
    else
        Xpsd = Xpsd - repmat(mean(Xpsd),Nf,1);
        Xcoh = Xcoh - repmat(mean(Xcoh),Nf,1);
        psd = (Xpsd' * sp_profiles); % ndim by componentss
        coh = (Xcoh' * sp_profiles);
    end
    
    for k = 1:K
        sp_fit_group.state(k).psd = zeros(options.Ncomp,ndim,ndim);
        sp_fit_group.state(k).coh = ones(options.Ncomp,ndim,ndim);
        ind = (1:ndim) + (k-1)*ndim;
        for i = 1:options.Ncomp
            sp_fit_group.state(k).psd(i,:,:) = diag(psd(ind,i));
        end
        ind = (1:ndim4) + (k-1)*ndim4; % This should exist if the code is not wrong
        for i = 1:options.Ncomp
            graphmat = zeros(ndim);
            graphmat(ind_offdiag) = coh(ind,i);
            graphmat=(graphmat+graphmat') + eye(ndim);
            sp_fit_group.state(k).coh(i,:,:) = graphmat;
        end
    end
        
    
   

%% SUBJECT LEVEL
if N>=1
    for n = 1:N
         %% CREATE MATRICES FOR FACTORISATION
        % put coh and psd in temporary arrays
        coh_comps = zeros(K,Nf,ndim1,ndim2);

        % for psd we will always keep the all the channels
        psd_comps = zeros(K,Nf,ndim);

        % accumulating psd coherence values state wise and subject wise
        
            for k = 1:K
                % For a given subject select the proper channels
                psd = sp_fit{n}.state(k).psd(:,chan,chan);
                coh = sp_fit{n}.state(k).coh(:,chan1,chan2);
                % The selected indexes will now begin from 1

                for j = 1:1:length(chan)
                    psd_comps(k,:,j) = psd(:,j,j);
        %             for l = 1:1:length(chan2)
        %                 coh_comps(n,k,:,j,l) = coh(:,j,l);
        %             end
                end

                % The selected indexes will now begin from 1
                for j = 1:1:length(chan1)

                    for l = 1:1:length(chan2)
                        coh_comps(k,:,j,l) = coh(:,j,l);
                    end
                end

            end
               
        % prepare matrix
        Xpsd = zeros(Nf,ndim*K);
        for k=1:K
            ind = (1:ndim) + (k-1)*ndim;
            Xpsd(:,ind)= squeeze(abs(psd_comps(k,:,:)));
        end        
        
        if assym_flag
            Xcoh = zeros(Nf,K*ndim_3);
            for k = 1:K
                ind = (1:ndim_3) + (k-1)*ndim_3;
                % Taking mean across subjects
                ck = squeeze((abs(coh_comps(k,:,:,:))));
            %     Xcoh(:,ind) = ck(:,ind_offdiag);
                ck = reshape(ck,[freq_bins,ndim1*ndim2]);
                Xcoh(:,ind) = ck;
            end
        else
                % This calculates total number of pairs
                % I know we could have replaced this with ndim_3 but for safety sake we
                % will keep the original part of the code
                ndim4 = ndim*(ndim-1)/2; 
                Xcoh = zeros(Nf,K*ndim4);
                for k = 1:K
                    ind = (1:ndim4) + (k-1)*ndim4;
                    ck = squeeze(mean(abs(coh_comps(k,:,:,:)),1));
                    Xcoh(:,ind) = ck(:,ind_offdiag);
                end
        end        
       
        % NNMF / PCA
        if strcmpi(options.Method,'NNMF')
            opt = statset('maxiter',1);
            if strcmpi(options.Base,'coh')
                 [a,b] = nnmf(Xcoh,options.Ncomp,'algorithm','als',...
                'w0',sp_profiles,'Options',opt);
                 sub_profile = a;
            else
                [a,b] = nnmf(Xpsd,options.Ncomp,'algorithm','als',...
                    'w0',sp_profiles,'Options',opt);
                sub_profile = a;               
            end


        else
            Xpsd = Xpsd - repmat(mean(Xpsd),Nf,1);
            Xcoh = Xcoh - repmat(mean(Xcoh),Nf,1);
            psd = (Xpsd' * sp_profiles);
            coh = (Xcoh' * sp_profiles);
        end        
        
        
        %% PROJECTION
    coh_comps = zeros(K,Nf,ndim,ndim);
    psd_comps = zeros(K,Nf,ndim);
    for k = 1:K
        % For a given subject select the proper channels
        psd = sp_fit{n}.state(k).psd(:,chan,chan);
        coh = sp_fit{n}.state(k).coh(:,chan,chan);
        % The selected indexes will now begin from 1
        % PSD is self it is not cross spectral across 
        for j = 1:1:length(chan)
            psd_comps(k,:,j) = psd(:,j,j);
            for l = 1:1:length(chan)
                coh_comps(k,:,j,l) = coh(:,j,l);
            end
        end
    end
    
    
    Xpsd = zeros(Nf,ndim*K);
    for k=1:K
        ind = (1:ndim) + (k-1)*ndim;
        Xpsd(:,ind)= squeeze((abs(psd_comps(k,:,:))));
    end
    Xcoh = zeros(Nf,K*ndim4);
    for k = 1:K
        ind = (1:ndim4) + (k-1)*ndim4;
        ck = squeeze((abs(coh_comps(k,:,:,:))));
        Xcoh(:,ind) = ck(:,ind_offdiag);
    end
    
    if strcmpi(options.Method,'NNMF')
       psd = (Xpsd' * a); % ndim by components
       coh = (Xcoh' * a);
    else
        Xpsd = Xpsd - repmat(mean(Xpsd),Nf,1);
        Xcoh = Xcoh - repmat(mean(Xcoh),Nf,1);
        psd = (Xpsd' * a); % ndim by componentss
        coh = (Xcoh' * a);
    end
    
    
        % Reshape stuff
            for k = 1:K
                sp_fit{n}.state(k).psd = zeros(options.Ncomp,ndim,ndim);
                sp_fit{n}.state(k).coh = ones(options.Ncomp,ndim,ndim);
                ind = (1:ndim) + (k-1)*ndim;
                for i = 1:options.Ncomp
                    sp_fit{n}.state(k).psd(i,:,:) = diag(psd(ind,i));
                end
                ind = (1:ndim4) + (k-1)*ndim4;
                for i = 1:options.Ncomp
                    graphmat = zeros(ndim);
                    graphmat(ind_offdiag) = coh(ind,i);
                    graphmat=(graphmat+graphmat') + eye(ndim);
                    sp_fit{n}.state(k).coh(i,:,:) = graphmat;
                end
            end            
 
    end
    
else
    sp_fit = sp_fit_group;    
end   
    
end

end