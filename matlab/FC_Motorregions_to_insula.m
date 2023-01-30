
subnames = {'SIC01','SIC02','SIC03','ME01','ME02','ME03','ME04'};%


subnetwork_vals = [1.5 17 10 11];%SCAN Foot Hand Face

structs_ofinterest_all = {{'S_circular_insula_sup','G_insular_short'} {'L_S_circular_insula_sup','L_G_insular_short','R_S_circular_insula_sup','R_G_insular_short'}}; % 

FCvals = zeros(length(subnetwork_vals),length(subnames));


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
            
            aparc_L_file = [fslrdir 'MSC02.L.aparc.a2009s.32k_fs_LR.label.gii'];
            aparc_R_file = [fslrdir 'MSC02.R.aparc.a2009s.32k_fs_LR.label.gii'];
            
            lsurffile = [fslrdir 'MSC02.L.midthickness.32k_fs_LR.surf.gii'];
            rsurffile = [fslrdir 'MSC02.R.midthickness.32k_fs_LR.surf.gii'];
            
            dmatname =  ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/SIC01/bold1_222/cifti_distances/SIC01_distmat_surf_geodesic_vol_euclidean_xhem1000_uint8.mat'];%
            
        else
            
            tmasks = smartload([basedir subname '/tmask.mat']);%
            
            scanlist = textread([basedir subname '/cast_scans.txt'],'%s');%
            
            fslrdir = ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/' subname '/7112b_fs_LR/fsaverage_LR32k/'];
            
            aparc_L_file = [fslrdir subname '.L.aparc.a2009s.32k_fs_LR.label.gii'];
            aparc_R_file = [fslrdir subname '.R.aparc.a2009s.32k_fs_LR.label.gii'];
            
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
        
        
        aparc_L_file = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.L.aparc.a2009s.32k_fs_LR.label.gii'];
        aparc_R_file = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.R.aparc.a2009s.32k_fs_LR.label.gii'];
        
        lsurffile = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.L.midthickness.32k_fs_LR.surf.gii'];
        rsurffile = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.R.midthickness.32k_fs_LR.surf.gii'];
        
    end
    
    
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
    aparc_LR_target = [aparc_L_mask ; aparc_R_mask];
    aparc_LR_target = aparc_LR_target(motorspots.brainstructure(motorspots.brainstructure<3)>0);
    aparc_LR_target((end+1):(nnz(motorspots.brainstructure>0))) = false;
    
    
    meantcs = zeros(size(data.data,2),length(subnetwork_vals));
    for v = 1:length(subnetwork_vals)
        meantcs(:,v) = mean(data.data(motorspots.data==subnetwork_vals(v),:),1);
    end
    insulatc = mean(data.data(logical(aparc_LR_target),:),1)';
    
    corrmat = FisherTransform(paircorr_mod(meantcs,insulatc));
    
    FCvals(:,subnum) = corrmat;
    
    
    
    
end

%%




%%
disp('CMI vs foot')
[~,P,~,STATS] = ttest(FCvals(1,:),FCvals(2,:))
disp('CMI vs hand')
[~,P,~,STATS] = ttest(FCvals(1,:),FCvals(3,:))
disp('CMI vs face')
[~,P,~,STATS] = ttest(FCvals(1,:),FCvals(4,:))
