% For each subject, compute FC between the portions of each motor region in BA4
% and nearby vertices in BAs 1 and 3.



subnames = {'SIC01','SIC02','SIC03','ME01','ME02','ME03','ME04'};%


subnetwork_vals = [1.5 17 10 11];% SCAN Foot Hand Face

precentral_structs_ofinterest = [5 6];
postcentral_structs_ofinterest = [4 1];

FCvals = zeros(length(subnetwork_vals),length(subnetwork_vals),length(subnames));

mask_with_aparc = true;
structs_ofinterest_all = {{'precentral','paracentral' 'postcentral'} { 'L_precentral' 'L_paracentral'  'R_precentral' 'R_paracentral' 'L_postcentral'  'R_postcentral'}}; %




for subnum = 1:length(subnames)
    
    subname = subnames{subnum};
    disp(subname)
    
    if strcmp(subname(1:3),'SIC')
        structs_ofinterest = structs_ofinterest_all{1};
        basedir = '/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/';
        if strcmp(subname,'SIC01')
            
            tmasks = smartload([basedir subname '/onoff_tmask.mat']);%
            
            scanlist = textread([basedir subname '/' subname '_cast_onoff.txt'],'%s');%
            
            fslrdir = ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/' subname '/7112b_fs_LR_old/fsaverage_LR32k/'];
            
            aparc_L_file = [fslrdir 'MSC02.L.aparc.32k_fs_LR.label.gii'];
            aparc_R_file = [fslrdir 'MSC02.R.aparc.32k_fs_LR.label.gii'];
            
            lsurffile = [fslrdir 'MSC02.L.midthickness.32k_fs_LR.surf.gii'];
            rsurffile = [fslrdir 'MSC02.R.midthickness.32k_fs_LR.surf.gii'];
            
            dmatname =  ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/SIC01/bold1_222/cifti_distances/SIC01_distmat_surf_geodesic_vol_euclidean_xhem1000_uint8.mat'];%
            
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
        cd(infomapdir)
        motorspots = ft_read_cifti_mod([infomapdir subname '_rawassn_minsize10_regularized_networksplus_motorrestricted_CS.dscalar.nii']);%['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' subname '_precast_infomap_wacky2_subcortreg_ignoreverts/' subname '_spots_effectors_templatematch.dtseries.nii']);%'_spots_effectors_only_cort_templatematch_subcort_wta.dtseries.nii']);%
        motorspots.data(59413:end,:) = 0;
        %networks = ft_read_cifti_mod([infomapdir subname '_rawassn_minsize10_regularized_recolored.dscalar.nii']);
        
        %make them the same color
        motorspots.data(abs(motorspots.data-2.5)<.01) = 1.5;
        motorspots.data(abs(motorspots.data-6.6)<.01) = 1.5;
        
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
        
        structs_ofinterest = structs_ofinterest_all{2};
        basedir = ['/data/nil-bluearc/GMT/Evan/subjects/' subname '/'];
        infomapdir = [basedir 'infomap/REST_adaptive_moreverts_s1p7_subcortregressed/'];
        cd([infomapdir])
        dmatname = [basedir '/anat/MNINonLinear/fsaverage_LR32k/distances/normalwall_distmat_surf_geodesic_vol_euclidean_xhemlarge_uint8.mat'];
        data = ft_read_cifti_mod([basedir '/func/rest/ConcatenatedCiftis/Rest_OCME+MEICA+MGTR_s1.7_MotionCensored+Concatenated.dtseries.nii']);
        
        tmask_concat = smartload([basedir 'func/rest/tmasks/Tmask_' subname '.mat']);
        tmask_concat = logical(tmask_concat);
        sessions = smartload([basedir 'func/rest/tmasks/ScanIdx_' subname '.mat']);
        %sessions = ones(size(sessions));
        
        motorspots = ft_read_cifti_mod([infomapdir subname '_rawassn_minsize10_regularized_CONandmotor_oneID_CS.dtseries.nii']);%[infomapdir subname '_rawassn_minsize10_regularized_networksandmotorspots.dtseries.nii']);
        motorspots.data(59413:end,:) = 0;
        
        motorspots.data(abs(motorspots.data-2.5)<.01) = 1.5;
        motorspots.data(abs(motorspots.data-6.6)<.01) = 1.5;
        
        
        aparc_L_file = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.L.aparc.32k_fs_LR.label.gii'];
        aparc_R_file = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.R.aparc.32k_fs_LR.label.gii'];
        
        lsurffile = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.L.midthickness.32k_fs_LR.surf.gii'];
        rsurffile = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.R.midthickness.32k_fs_LR.surf.gii'];
        
    end
    
    motorspots_valsofinterest = motorspots.data .* (any(repmat(motorspots.data,1,length(subnetwork_vals))==repmat(subnetwork_vals,length(motorspots.data),1),2));
    
    BA = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/CIMT/BAs/' subname '.LR.BA.32k_fs_LR.dtseries.nii']);
    
    BA_pre = zeros(size(BA.data));
    for i = 1:length(precentral_structs_ofinterest)
        BA_pre(BA.data==precentral_structs_ofinterest(i)) = 1;
    end
    
    BA_post = zeros(size(BA.data));
    for i = 1:length(postcentral_structs_ofinterest)
        BA_post(BA.data==postcentral_structs_ofinterest(i)) = 1;
    end
    
    lsurf = gifti(lsurffile);
    rsurf = gifti(rsurffile);
    LRcoords = [lsurf.vertices ; rsurf.vertices];
    LRcoords = LRcoords(motorspots.brainstructure(motorspots.brainstructure<3)>0,:);
    
    
    %%medial SM cortex isn't actually properly divided, so let's hack it.
    SMfootinds = find(motorspots_valsofinterest==17 & (logical(BA_pre) | logical(BA_post)));
    SMfoot_Y_coords = LRcoords(SMfootinds,2);
    medianYcoord = median(unique(SMfoot_Y_coords));
    BA_pre(SMfootinds) = 0;
    BA_pre(SMfootinds(SMfoot_Y_coords>=medianYcoord)) = 1;
    BA_post(SMfootinds) = 0;
    BA_post(SMfootinds(SMfoot_Y_coords<medianYcoord)) = 1;
    
    
    if mask_with_aparc
        aparc_L = gifti(aparc_L_file);
        aparc_L_mask = false(size(aparc_L.cdata));
        aparc_R = gifti(aparc_R_file);
        aparc_R_mask = false(size(aparc_R.cdata));
        
        
        for s=1:length(structs_ofinterest)
            Lind = strcmp(aparc_L.labels.name,structs_ofinterest{s});
            if any(Lind)
                Lkey = aparc_L.labels.key(Lind);
                aparc_L_mask(aparc_L.cdata==Lkey) = true;
            end
            
            Rind = strcmp(aparc_R.labels.name,structs_ofinterest{s});
            if any(Rind)
                Rkey = aparc_R.labels.key(Rind);
                aparc_R_mask(aparc_R.cdata==Rkey) = true;
            end
        end
        aparc_LR_mask = [aparc_L_mask ; aparc_R_mask];
        aparc_LR_mask = aparc_LR_mask(motorspots.brainstructure(motorspots.brainstructure<3)>0);
        aparc_LR_mask((end+1):(nnz(motorspots.brainstructure>0))) = false;
        
        
        BA_pre = BA_pre .* single(aparc_LR_mask);
        BA_post = BA_post .* single(aparc_LR_mask);
        motorspots_valsofinterest = motorspots_valsofinterest .* aparc_LR_mask;
        
    end
    
    BA_pre_vals = motorspots_valsofinterest .* BA_pre;
    BA_pre_inds = find(BA_pre);
    BA_post_inds = find(BA_post);
    
    
    
    distances = smartload(dmatname);
    
    [mins,mini] = min(distances(logical(BA_pre),logical(BA_post)),[],1);
    BA_post_vals = zeros(size(motorspots_valsofinterest));
    
    BA_post_vals(logical(BA_post)) = BA_pre_vals(BA_pre_inds(mini));
    
    
    
    
    pretcs = zeros(size(data.data,2),length(subnetwork_vals));
    posttcs = zeros(size(data.data,2),length(subnetwork_vals));
    for v = 1:length(subnetwork_vals)
        pretcs(:,v) = mean(data.data(abs(BA_pre_vals - subnetwork_vals(v))<.01,:),1);
        posttcs(:,v) = mean(data.data(abs(BA_post_vals - subnetwork_vals(v))<.01,:),1);
        
    end
    
    
    corrmat = FisherTransform(paircorr_mod(pretcs,posttcs));
    
    FCvals(:,:,subnum) = corrmat;
    
    
    
    
end

%%
power_surf_colormap = [1 0 0;0 0 .8;1 1 0;1 .8 .6;0 1 0;1 .6 1;0 .6 .6;0 0 0;.35 0 .65;.2 1 1;1 .5 0;.65 .25 1;0 .25 .6;.6 1 .6;.2 .3 1;1 1 1;0 .4 0; repmat([.25 .25 .25],50,1)];




allcolors = cell(1,length(subnetwork_vals));

for v = 1:length(subnetwork_vals)
    
    decimalval = mod(subnetwork_vals(v),1);
    if decimalval==0
        thiscolor = power_surf_colormap(subnetwork_vals(v),:);
    else
        thiscolor = sum([power_surf_colormap(floor(subnetwork_vals(v)),:) .* (1-decimalval) ; power_surf_colormap(ceil(subnetwork_vals(v)),:) .* (decimalval)],1);
    end
    allcolors{v} = thiscolor;
end

%%
figure;
set(gcf,'Position',[812 165 1027 805])
set(gcf,'Color',[1 1 1]);%set(gcf,'Color',[.9 .9 .9])
set(gca,'Color',[1 1 1]);%set(gca,'Color',[.9 .9 .9])



subsind = [1:7];

diags = zeros(length(subsind),length(subnetwork_vals));
for s = 1:length(subsind)
    diags(s,:) = diag(FCvals(:,:,subsind(s)));
end


categoryIdx = repmat([1:size(diags,2)]',1,length(subsind));
plotSpread(diags','categoryIdx',categoryIdx,'categoryColors',allcolors,'MarkerSize',40,'distributionMarkers','o')

h = gca;
for i = 1:length(h.Children)
    h.Children(i).MarkerFaceColor = h.Children(i).Color;
    h.Children(i).Color = [0 0 0];
    h.Children(i).LineWidth = 2.5;
end

set(gca,'FontSize',30)
set(gca,'XTickLabels',[])
box off
ylim([.4 2])

%%
disp('SCAN vs foot')
[~,P,~,STATS] = ttest(diags(:,1),diags(:,2))

disp('SCAN vs hand')
[~,P,~,STATS] = ttest(diags(:,1),diags(:,3))

disp('SCAN vs mouth')
[~,P,~,STATS] = ttest(diags(:,1),diags(:,4))

disp('foot vs hand')
[~,P,~,STATS] = ttest(diags(:,2),diags(:,3))

disp('foot vs mouth')
[~,P,~,STATS] = ttest(diags(:,2),diags(:,4))

disp('hand vs mouth')
[~,P,~,STATS] = ttest(diags(:,3),diags(:,4))
