%% Main lag script used to compute time delay matrix and lag projection
%
%
addpath /data/nil-bluearc/GMT/Evan/Scripts/lag-code-master/
lag_lim = 4;    % lag limit (in seconds)
lags = -3:3;    % range of TR shifts; max(lags) = round(lag_lim/tr + 1)
trs = [2.2 1.355];
xdist = 30;

subnames = {'SIC01','SIC02','SIC03','ME01','ME02','ME03','ME04'};%

networks_toinclude = [9 1.5 17 10 11]; %CON SCAN Foot Hand Face



brainstructures_touse_groups = {[1 2 9 10],[1 2 10 11]};%{[1 2 9 10 7 8 13:18],[1 2 10 11 8 9 16:21]};%{[1 2 9 10 15:18],[1 2 10 11 18:21]};%
brainstructure_labels = {'Cx','Cx','Cb','Cb'};%{'Cx','Cx','Cb','Cb','Str','Str','Str','Str','Str','Str','Th','Th'}; %'GP','GP',
brainstructure_labels_collapsed = {'Cx','Cb'};%{'Cx','Cb','Str','Th'};%

%%


for subnum = 1:length(subnames)
    
    subname = subnames{subnum};
    disp(subname)
    
    if strcmp(subname(1:3),'SIC')
       
        basedir = '/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/';
        tr = trs(1);
        brainstructures_touse = brainstructures_touse_groups{1};
        if strcmp(subname,'SIC01')
            
            tmasks = smartload([basedir subname '/onoff_tmask.mat']);%
            
            scanlist = textread([basedir subname '/' subname '_cast_onoff.txt'],'%s');%
            
            dmatname =  ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/SIC01/bold1_222/cifti_distances/SIC01_distmat_surf_geodesic_vol_euclidean_xhem1000_uint8.mat'];%
            
            fslrdir = ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/' subname '/7112b_fs_LR_old/fsaverage_LR32k/'];
            
            aparc_L_file = [fslrdir 'MSC02.L.aparc.32k_fs_LR.label.gii'];
            aparc_R_file = [fslrdir 'MSC02.R.aparc.32k_fs_LR.label.gii'];
            
             lsurffile = [fslrdir 'MSC02.L.midthickness.32k_fs_LR.surf.gii'];
            rsurffile = [fslrdir 'MSC02.R.midthickness.32k_fs_LR.surf.gii'];
            
        else
            
            tmasks = smartload([basedir subname '/tmask.mat']);%
            
            scanlist = textread([basedir subname '/cast_scans.txt'],'%s');%
            
            fslrdir = ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/' subname '/7112b_fs_LR/fsaverage_LR32k/'];
            
            aparc_L_file = [fslrdir subname '.L.aparc.32k_fs_LR.label.gii'];
            aparc_R_file = [fslrdir subname '.R.aparc.32k_fs_LR.label.gii'];
            
            lsurffile = [fslrdir subname '.L.midthickness.32k_fs_LR.surf.gii'];
            rsurffile = [fslrdir subname '.R.midthickness.32k_fs_LR.surf.gii'];
            
            
            if strcmp(subname,'SIC02')
                dmatname =  ['/data/nil-bluearc/GMT/Laumann/MSC/MSM_nativeresampled2_TYNDC/MSC06/fsaverage_LR32k/MSC06_distmat_surf_geodesic_vol_euclidean_xhem1000_uint8.mat'];%['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/SIC01/bold1_222/cifti_distances/SIC01_distmat_surf_geodesic_vol_euclidean_xhem1000_uint8.mat'];%
            else
                dmatname = ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/' subname '/bold1_222/cifti_distances/' subname '_distmat_surf_geodesic_vol_euclidean_xhem1000_uint8.mat'];
                
            end
            
        end
        
        infomapdir = ['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' subname '_precast_infomap_wacky2_subcortreg_ignoreverts/'];
        
        motorspots = ft_read_cifti_mod([infomapdir subname '_rawassn_minsize10_regularized_networksplus_motorrestricted_CS.dscalar.nii']);
        %networks = ft_read_cifti_mod([infomapdir subname '_rawassn_minsize10_regularized_recolored.dscalar.nii']);
    
        
        datafolder = [basedir subname '/bold1_222/'];
        
        scanstouse_inds = 1:12; %pre-cast
        
        scanlist = scanlist(scanstouse_inds);
        
        for scanindnum = 1:length(scanlist)
            scanind = scanstouse_inds(scanindnum);
            scanname = scanlist{scanind};
            ciftifile = dir([datafolder scanname '*surfsmooth2.55_volsmooth2.dtseries.nii']);
            
            data = ft_read_cifti_mod([datafolder ciftifile(1).name]);
            data.data = data.data(:,logical(tmasks(scanind,:)));
            
            if scanindnum==1
                alldata = data;
                
            else
                
                alldata.data = [alldata.data data.data];
                
            end
            
            clear data
        end
        
        data = alldata;
        clear alldata;
        
        sessions = [];
        tmask_concat = [];
        inds_withintmaskedconcat = [];
        prevind = 0;
        for s = 1:numel(scanlist)
            tmask = tmasks(s,:)';
            sessions = [sessions ; repmat(s,size(tmask,1),1)];
            tmask_concat = [tmask_concat ; tmask];
            inds_withintmasked = zeros(length(tmask),1);
            inds_withintmasked(tmask) = [(prevind+1) : (prevind + nnz(tmask))];
            inds_withintmaskedconcat = [inds_withintmaskedconcat ; inds_withintmasked];
            prevind = prevind + nnz(tmask);
        end
        
        tmask_concat = logical(tmask_concat);
        
        
        
    else
        
        
        tr = trs(2);
        brainstructures_touse = brainstructures_touse_groups{2};
        basedir = ['/data/nil-bluearc/GMT/Evan/subjects/' subname '/'];
        infomapdir = [basedir 'infomap/REST_adaptive_moreverts_s1p7_subcortregressed/'];
        dmatname = [basedir '/anat/MNINonLinear/fsaverage_LR32k/distances/normalwall_distmat_surf_geodesic_vol_euclidean_xhemlarge_uint8.mat'];
        cd([infomapdir])
        
        data = ft_read_cifti_mod([basedir '/func/rest/ConcatenatedCiftis/Rest_OCME+MEICA+MGTR_s1.7_MotionCensored+Concatenated.dtseries.nii']);
        
        tmask_concat = smartload([basedir 'func/rest/tmasks/Tmask_' subname '.mat']);
        tmask_concat = logical(tmask_concat);
        sessions = smartload([basedir 'func/rest/tmasks/ScanIdx_' subname '.mat']);
        %sessions = ones(size(sessions));
        
        motorspots = ft_read_cifti_mod([infomapdir subname '_rawassn_minsize10_regularized_CONandmotor_oneID_CS.dtseries.nii']);%[infomapdir subname '_rawassn_minsize10_regularized_networksandmotorspots.dtseries.nii']);
        %networks = ft_read_cifti_mod([basedir 'infomap/REST_adaptive_moreverts_s1p7_subcortregressed/rawassn_minsize10_regularized_recolored.dscalar.nii']);
    
        
        aparc_L_file = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.L.aparc.32k_fs_LR.label.gii'];
        aparc_R_file = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.R.aparc.32k_fs_LR.label.gii'];
        
        lsurffile = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.L.midthickness.32k_fs_LR.surf.gii'];
        rsurffile = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.R.midthickness.32k_fs_LR.surf.gii'];
        
    end
    
    %make them the same color
    motorspots.data(abs(motorspots.data-1.5)<.01) = 1.5;
    motorspots.data(abs(motorspots.data-2.5)<.01) = 1.5;
    motorspots.data(abs(motorspots.data-6.6)<.01) = 1.5;
    
    verts_touse_structure = zeros(size(motorspots.data));
    for s = 1:length(brainstructures_touse)
        verts_touse_structure(motorspots.brainstructure(motorspots.brainstructure>0)==brainstructures_touse(s)) = s;
    end
    
   
    verts_touse_network = zeros(size(motorspots.data));
    for n = 1:length(networks_toinclude)
        verts_touse_network(motorspots.data==networks_toinclude(n)) = networks_toinclude(n);
    end
    
    verts_touse = logical(verts_touse_network .* single(logical(verts_touse_structure)));
    
    verts_touse_structure = verts_touse_structure(verts_touse);
    verts_touse_network = verts_touse_network(verts_touse);
    
    num_nodes = nnz(verts_touse);
    
        
    
    
    %tmask_concat = logical(load(tmaskfile));
    %sessions = what; %un-tmasked
    sessionIDs = unique(sessions);
    inds_withintmaskedconcat = zeros(size(tmask_concat));
    inds_withintmaskedconcat(tmask_concat) = 1:nnz(tmask_concat);
    
    %all_lagmat = nan(num_nodes,num_nodes,length(sessionsIDs));
    
    
    
    %% Setup
    % Set parameters
    
    %mkdir(outdir)
    
    
    % Specify data parameters
    
    
    min_block_durn = (max(lags)+1)*tr;   % min. block duration (in seconds)
    
    %% Loop over scanlist
    % initialize group matrices for running sum
    grp_lags = single(nan(num_nodes));    % peak lags
    grp_ZL = grp_lags;      % zero-lag correlation
    grp_peak = grp_lags;    % peak correlation
    
    grp_lags_nans = single(zeros(num_nodes));
    grp_ZL_nans = grp_lags_nans;
    grp_peak_nans = grp_lags_nans;
    
    for s = 1:length(sessionIDs)
        tic
        disp(['Processing session ' num2str(s)]);
        
        % initialize subject matrices
        subj_lags = single(nan(num_nodes)); % peak lags
        subj_ZL = subj_lags;   % zero-lag correlation
        subj_peak = subj_lags; % peak correlation (correlation at optimal subnetwork_IDslag)
        
        format = logical(tmask_concat(sessions == sessionIDs(s)));
        
        BOLD = zeros(length(format),num_nodes);
        BOLD(format,:) = data.data(verts_touse,sessions(tmask_concat)==s)';
        
        
%         for IDnum = 1:num_nodes
%             
%             these_inds = allclusters(:,IDnum);
%             
%             BOLD(format,IDnum) = mean(data.data(these_inds,sessions(tmask_concat)==s),1)';
%         end
        
        inds_withintmasked = inds_withintmaskedconcat(sessions==sessionIDs(s));
        
        
        good = true(num_nodes,1); % read in spatial mask if desired
        
        % ignore pre-steady-state frames
        format(1:2) = false; % ignore first X frames
        
        FORMAT = create_blocks(format,min_block_durn,tr);
        
        %% Do the lagged correlation/covariance computation of TD matrices
        Cov = zeros([sum(good) sum(good) numel(lags)]);
        nblocks = numel(FORMAT);
        nframes = 0;
        
        % De-mean time series
        run_mean = nanmean(BOLD(format,:),1);
        BOLD = bsxfun(@minus,BOLD,run_mean);
        
        % Loop over blocks of contiguous frames
        for j = 1:numel(FORMAT)
            nframes = nframes + numel(FORMAT{j});
            FHCR = false(1,numel(format));
            FHCR(FORMAT{j}) = true;
            Cov = Cov + lagged_cov(BOLD(FHCR,good),BOLD(FHCR,good),max(lags));
        end
        
        % Normalize pairwise cross-covariance functions based on entire run
        for k = 1:numel(lags)
            Cov(:,:,k) = Cov(:,:,k)/(nframes - abs(lags(k))*nblocks);
        end
        
        % Parabolic interpolation to get peak lag/correlation
        [pl,pc] = parabolic_interp(Cov,tr);
        pl(abs(pl) > lag_lim) = nan; % Exclude long lags (generally occur when CCF is flat)
        
        % Get zero-lag correlation
        temp = Cov(:,:,lags==0);  % zero-lag correlation
        d = zeros(size(temp));
        d(logical(eye(length(temp)))) = sqrt(diag(temp));
        temp = d^(-1)*temp/d;
        temp = atanh(temp); % Fisher z transform
        temp(isnan(pl)) = nan;
        
        % Add to group running sum
        subj_lags(good,good) = pl;
        subj_ZL(good,good) = temp;
        subj_peak(good,good) = pc;
        
        grp_lags = cat(3,grp_lags,subj_lags);
        grp_lags = nansum(grp_lags,3);
        grp_ZL = cat(3,grp_ZL,subj_ZL);
        grp_ZL = nansum(grp_ZL,3);
        grp_peak = cat(3,grp_peak,subj_peak);
        grp_peak = nansum(grp_peak,3);
        
        % running sum of nans
        grp_lags_nans = grp_lags_nans + isnan(subj_lags);
        grp_ZL_nans = grp_ZL_nans + isnan(subj_ZL);
        grp_peak_nans = grp_peak_nans + isnan(subj_peak);
        
        toc
        
    end
    
    % Compute group averages
    grp_lags_mean = grp_lags ./ (length(sessionIDs) - grp_lags_nans);
    grp_peak_mean = grp_peak ./ (length(sessionIDs) - grp_peak_nans);
    grp_ZL_mean = grp_ZL ./ (length(sessionIDs) - grp_ZL_nans);
    grp_ZL_mean = tanh(grp_ZL_mean); % un-fisher z transform
    
    %% Sort group matrices
    
    assns = [1:num_nodes]; % import ROI network assignments for sorting TD matrix
    
    % Sort by matrices by lag
    [M,sorted_inds1] = sort(nanmean(grp_lags_mean));
    assns_sort = assns(sorted_inds1);
    
    grp_lags_temp = grp_lags_mean(sorted_inds1,sorted_inds1);
    grp_peak_temp = grp_peak_mean(sorted_inds1,sorted_inds1);
    grp_ZL_temp = grp_ZL_mean(sorted_inds1,sorted_inds1);
    
    % Sort by network
    [N,sorted_inds2] = sort(assns_sort);
    sorted_inds2 = sorted_inds2(find(N,1):end);
    
    grp_lags_mat = grp_lags_temp(sorted_inds2,sorted_inds2);
    grp_peak_mat = grp_peak_temp(sorted_inds2,sorted_inds2);
    grp_ZL_mat = grp_ZL_temp(sorted_inds2,sorted_inds2);
    
    
    %figure;imagesc(grp_peak_mat,[-.7,.7]);colorbar
    %figure;imagesc(grp_ZL_mat,[-.7,.7]);colorbar
    
    %% Make group-level lag projection maps
    % Unweighted lag projection
    
%     dmat = smartload(dmatname);
%     dmat = dmat(verts_touse,verts_touse);
%     dmat = dmat > xdist;
    %grp_lags_mean(~dmat) = NaN;
    
    grp_lags_proj_unweighted = nanmean(grp_lags_mean);
    
    % Weighted lag projection (inversely weight lags by correlation magnitude
    % to reduce sampling error)
    lag_weights = tan((pi/2)*(1-abs(grp_ZL_mean))).^(-2);    % weighted by 1/f^2(r); f(r) = tan[(pi/2)(1-|r|)]
    lag_weights(logical(eye(size(lag_weights)))) = 0;   % zero-out diagonal weights
    lag_weights(isnan(grp_lags_mean)) = nan;
    grp_lags_mean_wghtd = grp_lags_mean.*lag_weights;
    grp_lags_proj = nansum(grp_lags_mean_wghtd)./nansum(lag_weights);
    
    
    
    
    allstructures = zeros(size(verts_touse_structure));
    for i = 1:length(allstructures)
        matchingstringinds = find(strcmp(brainstructure_labels,brainstructure_labels{verts_touse_structure(i)}));
        allstructures(i) = matchingstringinds(1);
        %allstructures(i) = find(brainstructures_touse==allstructures_temp(i));
    end
    allcolors = verts_touse_network;
    
    %save([infomapdir subname  '_all_lagmat_motorandCON_vertexwise.mat'],'grp_lags_mat','grp_lags_proj','allcolors','allstructures')
    
    lagsmap = motorspots;
    lagsmap.data = zeros(size(motorspots.data,1),2);
    lagsmap.dimord = 'pos_scalar';
    lagsmap.data(verts_touse,1) = grp_lags_proj_unweighted;
    lagsmap.mapname{1} = 'mean lags unweighted';
    lagsmap.data(verts_touse,2) = grp_lags_proj;
    lagsmap.mapname{2} = 'mean lags weighted';
    ft_write_cifti_mod([infomapdir subname  '_lags_motorandCON_vertexwise'],lagsmap)
    
    

    
    
    %%
    
    
    
end


%%
lagsorder_bysub_network_struct = zeros(length(subnames),length(networks_toinclude), length(brainstructure_labels_collapsed));
power_surf_colormap = [1 0 0;0 0 .8;1 1 0;1 .8 .6;0 1 0;1 .6 1;0 .6 .6;0 0 0;.35 0 .65;.2 1 1;1 .5 0;.65 .25 1;0 .25 .6;.6 1 .6;.2 .3 1;1 1 1;0 .4 0; repmat([.25 .25 .25],50,1)];


for subnum = 1:length(subnames)
    
    subname = subnames{subnum};
    disp(subname)
    
    if strcmp(subname(1:3),'SIC')
        
        infomapdir = ['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' subname '_precast_infomap_wacky2_subcortreg_ignoreverts/'];
        
        brainstructures_touse = brainstructures_touse_groups{1};
        
        
    else
        basedir = ['/data/nil-bluearc/GMT/Evan/subjects/' subname '/'];
        infomapdir = [basedir 'infomap/REST_adaptive_moreverts_s1p7_subcortregressed/'];
        brainstructures_touse = brainstructures_touse_groups{2};
        
    end
    
    motorspots = ft_read_cifti_mod([infomapdir '/' subname '_rawassn_minsize10_regularized_CONandmotor_oneID_central_precentral.dtseries.nii']);
    
    lagsmap = ft_read_cifti_mod([infomapdir subname  '_lags_motorandCON_vertexwise.dscalar.nii']);
    
    brainstructurenums = lagsmap.brainstructure(lagsmap.brainstructure>0);
    
    structnums = zeros(size(brainstructurenums));
    
    for s = 1:length(brainstructures_touse)
        
        corresponding_name = brainstructure_labels{s};
        this_structnum = find(strcmp(brainstructure_labels_collapsed,corresponding_name));
        this_structnum = this_structnum(1);
        structnums(brainstructurenums==brainstructures_touse(s)) = this_structnum;
    end
    
    for s = 1:length(brainstructure_labels_collapsed)
        for n = 1:length(networks_toinclude)
            
            lagsorder_bysub_network_struct(subnum,n,s) = mean(lagsmap.data(motorspots.data==networks_toinclude(n) & structnums==s,1));
        end
        
    end
    
end
        
sub_mean_lags = squeeze(nanmean(lagsorder_bysub_network_struct,1));
sub_stder_lags = squeeze(nanstd(lagsorder_bysub_network_struct,[],1)) ./ sqrt(length(subnames));
sub_mean_lags = sub_mean_lags';
sub_stder_lags = sub_stder_lags';
%network_sortorder = [1 2 3 4 5];


%sub_mean_lags = sub_mean_lags(network_sortorder,:)';
%sub_stder_lags = sub_stder_lags(network_sortorder,:)';
%networkIDs_resorted = networks_toinclude(network_sortorder);


%%


network_sortorder_disp = [1 2 5 4 3];
networkIDs_resorted_disp = networks_toinclude(network_sortorder_disp);%networkIDs_resorted(network_sortorder_disp);

figure;
set(gcf,'Position',[273 368 2011 467]);%[813 579 1471 256]);%[813 30 1102 805])
set(gcf,'Color',[1 1 1]);%set(gcf,'Color',[.9 .9 .9])
set(gca,'Color',[1 1 1]);
bar(sub_mean_lags(:,network_sortorder_disp))
h = gca;
for i = 1:length(h.Children)
    thisnumber = networkIDs_resorted_disp(i);
    decimalval = mod(thisnumber,1);
    if decimalval==0
        thiscolor = power_surf_colormap(thisnumber,:);
    else
        thiscolor = sum([power_surf_colormap(floor(thisnumber),:) .* (1-decimalval) ; power_surf_colormap(ceil(thisnumber),:) .* (decimalval)],1);
    end
    h.Children(length(h.Children)-i+1).FaceColor = thiscolor;
end
hold on
h.XColor = [0 0 0];
h.YColor = [0 0 0];
h.YAxis.LineWidth = 1.5;
h.XAxis.LineWidth = 1.5;

nbars = size(sub_mean_lags,2);
xcoords = [];
for i = 1:nbars
    xcoords = [xcoords ; h.Children(i).XEndPoints];
end
xcoords = xcoords';
xcoords = xcoords(:,end:-1:1);
eh = errorbar(xcoords,sub_mean_lags(:,network_sortorder_disp), sub_stder_lags(:,network_sortorder_disp),'k','linestyle','none','LineWidth',4);
for i = 1:length(eh)
    eh(i).CapSize = 30;
end

% for i = 1:length(brainstructure_labels_collapsed)
%     h.XTickLabel{i} = brainstructure_labels_collapsed{i};
% end
% title('Sub Avg')
set(gca,'FontSize',30);
xlim([.6 1.4])
set(gca,'XTick',[])
set(gca,'XTickLabel',[])
set(gca,'XAxisLocation','top');
box off

%%
pvec_cort= [];
disp('Cortex:')
for i = 1:3
[~,P,~,STATS] = ttest(lagsorder_bysub_network_struct(:,2,1),lagsorder_bysub_network_struct(:,i+2,1));
disp(['CMI vs network ' num2str(networks_toinclude(i+2)) ': T=' num2str(STATS.tstat) '; p=' num2str(P)])
pvec_cort(end+1) = P;
end
[~,P,~,STATS] = ttest(lagsorder_bysub_network_struct(:,2,1),mean(lagsorder_bysub_network_struct(:,3:end,1),2));
disp(['CMI vs mean effectors : T=' num2str(STATS.tstat) '; p=' num2str(P)])
[~,P,~,STATS] = ttest(lagsorder_bysub_network_struct(:,2,1),lagsorder_bysub_network_struct(:,1,1));
disp(['CMI vs CON : T=' num2str(STATS.tstat) '; p=' num2str(P)])

for i = 1:3
[~,P,~,STATS] = ttest(lagsorder_bysub_network_struct(:,1,1),lagsorder_bysub_network_struct(:,i+2,1));
disp(['CON vs network ' num2str(networks_toinclude(i+2)) ': T=' num2str(STATS.tstat) '; p=' num2str(P)])
pvec_cort(end+1) = P;
end
[~,P,~,STATS] = ttest(lagsorder_bysub_network_struct(:,1,1),mean(lagsorder_bysub_network_struct(:,3:end,1),2));
disp(['CON vs mean effectors : T=' num2str(STATS.tstat) '; p=' num2str(P)])

pvec_cbllm= [];
disp('Cerebellum:')
for i = 1:3
[~,P,~,STATS] = ttest(lagsorder_bysub_network_struct(:,2,2),lagsorder_bysub_network_struct(:,i+2,2));
disp(['CMI vs network ' num2str(networks_toinclude(i+2)) ': T=' num2str(STATS.tstat) '; p=' num2str(P)])
pvec_cbllm(end+1) = P;
end
[~,P,~,STATS] = ttest(lagsorder_bysub_network_struct(:,2,2),mean(lagsorder_bysub_network_struct(:,3:end,2),2));
disp(['CMI vs mean effectors : T=' num2str(STATS.tstat) '; p=' num2str(P)])
[~,P,~,STATS] = ttest(lagsorder_bysub_network_struct(:,2,2),lagsorder_bysub_network_struct(:,1,2));
disp(['CMI vs CON : T=' num2str(STATS.tstat) '; p=' num2str(P)])

for i = 1:3
[~,P,~,STATS] = ttest(lagsorder_bysub_network_struct(:,1,2),lagsorder_bysub_network_struct(:,i+2,2));
disp(['CON vs network ' num2str(networks_toinclude(i+2)) ': T=' num2str(STATS.tstat) '; p=' num2str(P)])
pvec_cbllm(end+1) = P;
end
[~,P,~,STATS] = ttest(lagsorder_bysub_network_struct(:,1,2),mean(lagsorder_bysub_network_struct(:,3:end,2),2));
disp(['CON vs mean effectors : T=' num2str(STATS.tstat) '; p=' num2str(P)])

