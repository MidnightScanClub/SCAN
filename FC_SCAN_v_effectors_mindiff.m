
subnames = {'SIC01','SIC02','SIC03','ME01','ME02','ME03','ME04'};%

subnetwork_vals = [1.5 17 10 11];%SCAN Foot Hand Face


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
            
        else
            
            tmasks = smartload([basedir subname '/tmask.mat']);%
            
            scanlist = textread([basedir subname '/cast_scans.txt'],'%s');%
            
            fslrdir = ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/' subname '/7112b_fs_LR/fsaverage_LR32k/'];
            
            aparc_L_file = [fslrdir subname '.L.aparc.32k_fs_LR.label.gii'];
            aparc_R_file = [fslrdir subname '.R.aparc.32k_fs_LR.label.gii'];
            
        end
        
        infomapdir = ['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' subname '_precast_infomap_wacky2_subcortreg_ignoreverts/'];
        cd(infomapdir)
        motorspots = ft_read_cifti_mod([infomapdir subname '_rawassn_minsize10_regularized_networksplus_motorrestricted_CS.dscalar.nii']);%['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' subname '_precast_infomap_wacky2_subcortreg_ignoreverts/' subname '_spots_effectors_templatematch.dtseries.nii']);%'_spots_effectors_only_cort_templatematch_subcort_wta.dtseries.nii']);%
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
        
        
        
        aparc_L_file = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.L.aparc.32k_fs_LR.label.gii'];
        aparc_R_file = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.R.aparc.32k_fs_LR.label.gii'];
        
    end
    
    motorspots.data(abs(motorspots.data-1.5)<.01) = 1.5;
    motorspots.data(abs(motorspots.data-2.5)<.01) = 1.5;
    motorspots.data(abs(motorspots.data-6.6)<.01) = 1.5;
    
    
    
    
        %% spots as separate subnetworks
    
    FCmaps = zeros(size(data.data,1),length(subnetwork_vals));
    for v = 1:length(subnetwork_vals)
        FCmaps(:,v) = FisherTransform(paircorr_mod(data.data',mean(data.data(abs(motorspots.data - subnetwork_vals(v))<.01,:),1)'));
    end
    
    out = data;
    out.data = FCmaps;
    ft_write_cifti_mod([infomapdir subname '_Allspots_and_Effectors_FC'],out)
    
    FCdifferencemaps = repmat(FCmaps(:,1),1,3) - FCmaps(:,2:end);
    out = data;
    out.data = FCdifferencemaps;
    ft_write_cifti_mod([infomapdir subname '_Allspots_v_Effectors'],out)
    
end

%%