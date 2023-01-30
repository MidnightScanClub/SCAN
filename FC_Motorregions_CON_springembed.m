subnames = {'SIC01','SIC02','SIC03','ME01','ME02','ME03','ME04'};%'SIC02',


xdistance = 30;

IDs_toinclude = [1.5 2.5 6.6 10 11 17 9]; %SCAN Mid SCAN Inf SCAN Sup Hand Face Foot CON

thresholdarray = .45;%[.25 : .05 : .5];


for subnum = 1:length(subnames)
    
    subname = subnames{subnum};
    disp(subname)
    
    if strcmp(subname(1:3),'SIC')
        
        basedir = '/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/';
        if strcmp(subname,'SIC01')
            
            tmasks = smartload([basedir subname '/onoff_tmask.mat']);%
            
            scanlist = textread([basedir subname '/' subname '_cast_onoff.txt'],'%s');%
            
            fslrdir = ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/' subname '/7112b_fs_LR_old/fsaverage_LR32k/'];
            
            aparc_L_file = [fslrdir 'MSC02.L.aparc.32k_fs_LR.label.gii'];
            aparc_R_file = [fslrdir 'MSC02.R.aparc.32k_fs_LR.label.gii'];
            
            atlasT1file = ['/data/nil-bluearc/GMT/Laumann/MSC/MSC02/T1/MSC02_mpr_debias_avgT_111_t88.nii.gz'];
            
            dmatname =  ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/SIC01/bold1_222/cifti_distances/SIC01_distmat_surf_geodesic_vol_euclidean_xhem1000_uint8.mat'];%
            
        else
            
            tmasks = smartload([basedir subname '/tmask.mat']);%
            
            scanlist = textread([basedir subname '/cast_scans.txt'],'%s');%
            
            fslrdir = ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/' subname '/7112b_fs_LR/fsaverage_LR32k/'];
            
            aparc_L_file = [fslrdir subname '.L.aparc.32k_fs_LR.label.gii'];
            aparc_R_file = [fslrdir subname '.R.aparc.32k_fs_LR.label.gii'];
            
            atlasT1file = ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/' subname '/T1/' subname '_mpr_debias_avgT_111_t88.nii.gz'];
            
            if strcmp(subname,'SIC02')
                dmatname =  ['/data/nil-bluearc/GMT/Laumann/MSC/MSM_nativeresampled2_TYNDC/MSC06/fsaverage_LR32k/MSC06_distmat_surf_geodesic_vol_euclidean_xhem1000_uint8.mat'];%['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/SIC01/bold1_222/cifti_distances/SIC01_distmat_surf_geodesic_vol_euclidean_xhem1000_uint8.mat'];%
            else
                dmatname = ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/' subname '/bold1_222/cifti_distances/' subname '_distmat_surf_geodesic_vol_euclidean_xhem1000_uint8.mat'];
            end
            
        end
        
        infomapdir = ['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' subname '_precast_infomap_wacky2_subcortreg_ignoreverts/'];
        cd(infomapdir)
        motorspots = ft_read_cifti_mod([infomapdir subname '_rawassn_minsize10_regularized_networksplus_motorrestricted_CS.dscalar.nii']);%['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' subname '_precast_infomap_wacky2_subcortreg_ignoreverts/' subname '_spots_effectors_templatematch.dtseries.nii']);%'_spots_effectors_only_cort_templatematch_subcort_wta.dtseries.nii']);%
        %motorspots = ft_read_cifti_mod([infomapdir subname '_spots_effectors_only_cort_templatematch_subcort_wta.dtseries.nii']);%['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' subname '_precast_infomap_wacky2_subcortreg_ignoreverts/' subname '_spots_effectors_templatematch.dtseries.nii']);%'_spots_effectors_only_cort_templatematch_subcort_wta.dtseries.nii']);%
        motorspots.data(59413:end,:) = 0;
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
        
        
        basedir = ['/data/nil-bluearc/GMT/Evan/subjects/' subname '/'];
        infomapdir = [basedir 'infomap/REST_adaptive_moreverts_s1p7_subcortregressed/'];
        cd([infomapdir])
        
        data = ft_read_cifti_mod([basedir '/func/rest/ConcatenatedCiftis/Rest_OCME+MEICA+MGTR_s1.7_MotionCensored+Concatenated.dtseries.nii']);
        
        tmask_concat = smartload([basedir 'func/rest/tmasks/Tmask_' subname '.mat']);
        tmask_concat = logical(tmask_concat);
        sessions = smartload([basedir 'func/rest/tmasks/ScanIdx_' subname '.mat']);
        %sessions = ones(size(sessions));
        
        motorspots = ft_read_cifti_mod([infomapdir subname '_rawassn_minsize10_regularized_networksplus_oneID_CS.dscalar.nii']);%[infomapdir subname '_rawassn_minsize10_regularized_networksandmotorspots.dtseries.nii']);
        motorspots.data(59413:end,:) = 0;
        networks = ft_read_cifti_mod([basedir 'infomap/REST_adaptive_moreverts_s1p7_subcortregressed/rawassn_minsize10_regularized_recolored.dscalar.nii']);
        
        
        aparc_L_file = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.L.aparc.32k_fs_LR.label.gii'];
        aparc_R_file = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.R.aparc.32k_fs_LR.label.gii'];
        
        atlasT1file = [basedir '/anat/T1w/T1w_acpc.nii.gz'];
        
        dmatname = [basedir '/anat/MNINonLinear/fsaverage_LR32k/distances/normalwall_distmat_surf_geodesic_vol_euclidean_xhemlarge_uint8.mat'];
        
    end
    
    
    %
    
    
    %
    motorspots.data(59413:end,:) = 0;
    
    %make them the same color
    motorspots.data(abs(motorspots.data-2.5)<.01) = 1.5;
    motorspots.data(abs(motorspots.data-6.6)<.01) = 1.5;
    
    IDs = unique(motorspots.data); IDs(IDs==0) = [];
    for IDnum = length(IDs) : -1 : 1
        if ~any(abs(IDs_toinclude - IDs(IDnum))<.01)
            IDs(IDnum) = [];
        end
    end
    
    allclusters = [];
    allclusterIDs = [];
    for IDnum = 1:length(IDs)
        ID = IDs(IDnum);
        clusters_thisID = cifti_cluster_separateanatstructures(motorspots,ID-.01,ID+.01,20);
        
        allclusters = [allclusters clusters_thisID];
        allclusterIDs = [allclusterIDs repmat(ID,1,size(clusters_thisID,2))];
    end
    
    
    cluster_sizes = sum(allclusters,1);
    
    allclusters = logical(allclusters);
    
    
    clustertcs = zeros(size(data.data,2),size(allclusters,2));
    for c = 1:size(allclusters,2)
        clustertcs(:,c) = mean(data.data(allclusters(:,c),:),1);
    end
    rmat = paircorr_mod(clustertcs);
    
    
    distances = smartload(dmatname);
    centroids = zeros(size(rmat,1),1);
    for c = 1:size(allclusters,2)
        clusterinds = find(allclusters(:,c));
        [~,mini] = min(sum(distances(clusterinds,clusterinds),2));
        centroids(c) = clusterinds(mini);
    end
    dmat = distances(centroids,centroids);
    
    rmat(dmat < xdistance) = 0;
    
    %%
    
    
    for t = 1:length(thresholdarray)
        
        sorted_corrvals = sort(rmat(:),'descend');
        thresh_val = sorted_corrvals(floor(length(sorted_corrvals).*thresholdarray(t)));
        
        graph = rmat >= thresh_val;
        
        Make_gephi_files_func_v3(graph,allclusterIDs',{'NodeID'},[subname '_Motor_CMI_CON_thr' num2str(thresholdarray(t))],[cluster_sizes'],{'Size'})
        
    end
    
end
    
