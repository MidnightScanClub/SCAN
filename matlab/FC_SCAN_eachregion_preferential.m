
subnames = {'SIC01','SIC02','SIC03','ME01','ME02','ME03','ME04'};%'};%,


subnetwork_vals = [1.5 2.5 6.6]; %Mid Inferior Superior


mask_with_aparc = true;
structs_ofinterest_all = {{'precentral','paracentral' } { 'L_precentral' 'L_paracentral'  'R_precentral' 'R_paracentral' }}; % 'postcentral' 'L_postcentral'  'R_postcentral'



for subnum = [1:length(subnames)]
    
    subname = subnames{subnum};
    disp(subname)
    
    if strcmp(subname(1:3),'SIC')
       
        basedir = '/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/';
        structs_ofinterest = structs_ofinterest_all{1};
        if strcmp(subname,'SIC01')
            
            tmasks = smartload([basedir subname '/onoff_tmask.mat']);%
            
            scanlist = textread([basedir subname '/' subname '_cast_onoff.txt'],'%s');%
            
            fslrdir = ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/' subname '/7112b_fs_LR_old/fsaverage_LR32k/'];
            
            aparc_L_file = [fslrdir 'MSC02.L.aparc.32k_fs_LR.label.gii'];
            aparc_R_file = [fslrdir 'MSC02.R.aparc.32k_fs_LR.label.gii'];
            
            atlasT1file = ['/data/nil-bluearc/GMT/Laumann/MSC/MSC02/T1/MSC02_mpr_debias_avgT_111_t88.nii.gz'];
            
        else
            
            tmasks = smartload([basedir subname '/tmask.mat']);%
            
            scanlist = textread([basedir subname '/cast_scans.txt'],'%s');%
            
            fslrdir = ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/' subname '/7112b_fs_LR/fsaverage_LR32k/'];
            
            aparc_L_file = [fslrdir subname '.L.aparc.32k_fs_LR.label.gii'];
            aparc_R_file = [fslrdir subname '.R.aparc.32k_fs_LR.label.gii'];
            
            atlasT1file = ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/' subname '/T1/' subname '_mpr_debias_avgT_111_t88.nii.gz'];
            
        end
        
        infomapdir = ['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' subname '_precast_infomap_wacky2_subcortreg_ignoreverts/'];
        cd(infomapdir)
        motorspots = ft_read_cifti_mod([infomapdir subname '_rawassn_minsize10_regularized_networksplus_motorrestricted_CS.dscalar.nii']);%['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' subname '_precast_infomap_wacky2_subcortreg_ignoreverts/' subname '_spots_effectors_templatematch.dtseries.nii']);%'_spots_effectors_only_cort_templatematch_subcort_wta.dtseries.nii']);%
        %motorspots = ft_read_cifti_mod([infomapdir subname '_spots_effectors_only_cort_templatematch_subcort_wta.dtseries.nii']);%['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' subname '_precast_infomap_wacky2_subcortreg_ignoreverts/' subname '_spots_effectors_templatematch.dtseries.nii']);%'_spots_effectors_only_cort_templatematch_subcort_wta.dtseries.nii']);%
        motorspots.data(59413:end,:) = 0;
        %networks = ft_read_cifti_mod([infomapdir subname '_rawassn_minsize10_regularized_recolored.dscalar.nii']);
        motorspots_orig = motorspots;
    
        
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
        
        motorspots_orig.data(:) = 0;
    for v = [subnetwork_vals(:)' 10 11 17]
        motorspots_orig.data(abs(motorspots.data-v)<.01) = v;
    end
        
    else
        
        
        structs_ofinterest = structs_ofinterest_all{2};
        basedir = ['/data/nil-bluearc/GMT/Evan/subjects/' subname '/'];
        infomapdir = [basedir 'infomap/REST_adaptive_moreverts_s1p7_subcortregressed/'];
        cd([infomapdir])
        
        data = ft_read_cifti_mod([basedir '/func/rest/ConcatenatedCiftis/Rest_OCME+MEICA+MGTR_s1.7_MotionCensored+Concatenated.dtseries.nii']);
        
        tmask_concat = smartload([basedir 'func/rest/tmasks/Tmask_' subname '.mat']);
        tmask_concat = logical(tmask_concat);
        sessions = smartload([basedir 'func/rest/tmasks/ScanIdx_' subname '.mat']);
        %sessions = ones(size(sessions));
        
        motorspots = ft_read_cifti_mod([infomapdir subname '_rawassn_minsize10_regularized_motorspotsclustered.dtseries.nii']);
        motorspots.data(59413:end,:) = 0;
        motorspots_orig = motorspots;
        
        moremotor = ft_read_cifti_mod([infomapdir subname '_rawassn_minsize10_regularized_CONandmotor_oneID_CS.dtseries.nii']);
        moremotor.data(59413:end,:) = 0;
        motorspots_orig.data(moremotor.data==10) = 10;
        motorspots_orig.data(moremotor.data==11) = 11;
        motorspots_orig.data(moremotor.data==17) = 17;
    
        
        aparc_L_file = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.L.aparc.32k_fs_LR.label.gii'];
        aparc_R_file = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.R.aparc.32k_fs_LR.label.gii'];
        
        atlasT1file = [basedir '/anat/T1w/T1w_acpc.nii.gz'];
        
    end
    
    
    motorspots_mask = motorspots_orig.data;
    
    
    
    
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
        
        
        motorspots.data(motorspots.brainstructure(motorspots.brainstructure>0)<3 & (~aparc_LR_mask)) = 0;
    end
    
    
    
        %% spots as separate subnetworks
    
    FCmaps = zeros(size(data.data,1),length(subnetwork_vals));
    out = data;
    out.dimord = 'pos_scalar';
    out.mapname = {['Spot ' num2str(subnetwork_vals(1))],['Spot ' num2str(subnetwork_vals(2))],['Spot ' num2str(subnetwork_vals(3))]};
    
    for v = 1:length(subnetwork_vals)
        FCmaps(:,v) = FisherTransform(paircorr_mod(data.data',mean(data.data(abs(motorspots.data-subnetwork_vals(v)) < .01,:),1)'));
    end
    FCmaps(isnan(FCmaps)) = 0;
    out.data = FCmaps;
    
    
    
    for v = 1:length(subnetwork_vals)
        otherinds = setdiff([1:length(subnetwork_vals)],v);
        maxotherFC = max(FCmaps(:,otherinds),[],2);
        out.data(:,length(subnetwork_vals)+v) = (FCmaps(:,v) - maxotherFC) .* single(FCmaps(:,v)>0) .* single(motorspots_mask==0);
        out.mapname{length(subnetwork_vals)+v} = ['Spot ' num2str(subnetwork_vals(v)) ' vs others and positive'];
    end
    
    ft_write_cifti_mod([subname 'motor_subnetworks_spotsisolated_Veachother_mindiff'],out)
    
    
    
    
    
    
    
    
end