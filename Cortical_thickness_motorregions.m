%For each subject, compute cortical thickness 


subnames = {'SIC01','SIC02','SIC03','ME01','ME02','ME03','ME04'};%};%,

values_ofinterest = [1.5 17 10 11];%SCAN Foot Hand Face

thicknesses = zeros(length(subnames),length(values_ofinterest));

mask_with_BA = true;
precentral_structs_ofinterest = [5 6];
postcentral_structs_ofinterest = [4 1];

mask_with_aparc = true;
structs_ofinterest_all = {{'precentral','paracentral' 'postcentral'} { 'L_precentral' 'L_paracentral'  'R_precentral' 'R_paracentral' 'L_postcentral'  'R_postcentral'}}; % 


for subnum = 1:length(subnames)
    
    subname = subnames{subnum};
    disp(subname)
    
    if strcmp(subname(1:3),'SIC')
       
        basedir = '/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/';
        if strcmp(subname,'SIC01')
            structs_ofinterest = structs_ofinterest_all{2};
            fslrdir = ['/data/nil-bluearc/GMT/Laumann/MSC/MSM_nativeresampled2_TYNDC/MSC02/fsaverage_LR32k//'];
            
            
            thickness_L_file = [fslrdir 'MSC02.L.thickness.32k_fs_LR.shape.gii'];
            thickness_R_file = [fslrdir 'MSC02.R.thickness.32k_fs_LR.shape.gii'];
            
            aparc_L_file = [fslrdir 'MSC02.L.aparc.32k_fs_LR.label.gii'];
            aparc_R_file = [fslrdir 'MSC02.R.aparc.32k_fs_LR.label.gii'];
            
            lsurffile = [fslrdir 'MSC02.L.midthickness.32k_fs_LR.surf.gii'];
            rsurffile = [fslrdir 'MSC02.R.midthickness.32k_fs_LR.surf.gii'];
            
            
        else
            structs_ofinterest = structs_ofinterest_all{1};
            fslrdir = ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/' subname '/7112b_fs_LR/fsaverage_LR32k/'];
            
            
            thickness_L_file = [fslrdir subname '.L.thickness.32k_fs_LR.shape.gii'];
            thickness_R_file = [fslrdir subname '.R.thickness.32k_fs_LR.shape.gii'];
            
            lsurffile = [fslrdir subname '.L.midthickness.32k_fs_LR.surf.gii'];
            rsurffile = [fslrdir subname '.R.midthickness.32k_fs_LR.surf.gii'];
            
            aparc_L_file = [fslrdir subname '.L.aparc.32k_fs_LR.label.gii'];
            aparc_R_file = [fslrdir subname '.R.aparc.32k_fs_LR.label.gii'];
            
            
            
        end
        
        infomapdir = ['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' subname '_precast_infomap_wacky2_subcortreg_ignoreverts/'];
        
        motorspots = ft_read_cifti_mod([infomapdir subname '_rawassn_minsize10_regularized_networksplus_motorrestricted_CS.dscalar.nii']);%['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' subname '_precast_infomap_wacky2_subcortreg_ignoreverts/' subname '_spots_effectors_templatematch.dtseries.nii']);%'_spots_effectors_only_cort_templatematch_subcort_wta.dtseries.nii']);%
        motorspots.data(59413:end,:) = 0;
        
        
        
        
    else
        
        structs_ofinterest = structs_ofinterest_all{2};
        basedir = ['/data/nil-bluearc/GMT/Evan/subjects/' subname '/'];
        infomapdir = [basedir 'infomap/REST_adaptive_moreverts_s1p7_subcortregressed/'];
        motorspots = ft_read_cifti_mod([infomapdir subname '_rawassn_minsize10_regularized_CONandmotor_oneID_CS.dtseries.nii']);%[infomapdir subname '_rawassn_minsize10_regularized_networksandmotorspots.dtseries.nii']);
        motorspots.data(59413:end,:) = 0;
        
        
        thickness_L_file = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.L.thickness.32k_fs_LR.shape.gii'];
        thickness_R_file = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.R.thickness.32k_fs_LR.shape.gii'];
        
        aparc_L_file = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.L.aparc.32k_fs_LR.label.gii'];
        aparc_R_file = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.R.aparc.32k_fs_LR.label.gii'];
        
        lsurffile = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.L.midthickness.32k_fs_LR.surf.gii'];
        rsurffile = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.R.midthickness.32k_fs_LR.surf.gii'];
        
        
    end
    
    
   
    
    thickness_L = gifti(thickness_L_file);
    thickness_R = gifti(thickness_R_file);
    thickness_LR = [thickness_L.cdata; thickness_R.cdata];
    thickness_LR = thickness_LR(motorspots.brainstructure(1:64984)>0);
    
%     %make them the same color
    motorspots.data(abs(motorspots.data-2.5)<.01) = 1.5;
    motorspots.data(abs(motorspots.data-6.6)<.01) = 1.5;
    
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
    SMfootinds = find(motorspots.data==17 & (logical(BA_pre) | logical(BA_post)));
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
    motorspots.data = motorspots.data .* aparc_LR_mask;
    
    end
    
    motorspots.data = motorspots.data .* BA_pre;
    

        
    
    for v = 1:length(values_ofinterest)
        thicknesses(subnum,v) = mean(thickness_LR(abs(motorspots.data-values_ofinterest(v))<.01));
    end
    
    if subnum==1
        out = motorspots;
        out.data(:) = 0;
        out.dimord = 'pos_scalar';
    end
    out.data(1:59412,subnum) = 0;
    for v = 1:(length(values_ofinterest))
        out.data(abs(motorspots.data-values_ofinterest(v))<.01,subnum) = thickness_LR(abs(motorspots.data-values_ofinterest(v))<.01);
    end
    out.mapname{subnum} = subname;
    
end

ft_write_cifti_mod('Thickness_in_precentral',out)
    

    
%%

power_surf_colormap = [1 0 0;0 0 .8;1 1 0;1 .8 .6;0 1 0;1 .6 1;0 .6 .6;0 0 0;.35 0 .65;.2 1 1;1 .5 0;.65 .25 1;0 .25 .6;.6 1 .6;.2 .3 1;1 1 1;0 .4 0; repmat([.25 .25 .25],50,1)];

allcolors = cell(1,length(values_ofinterest));

for v = 1:length(values_ofinterest)
    
    decimalval = mod(values_ofinterest(v),1);
    if decimalval==0
        thiscolor = power_surf_colormap(values_ofinterest(v),:);
    else
        thiscolor = sum([power_surf_colormap(floor(values_ofinterest(v)),:) .* (1-decimalval) ; power_surf_colormap(ceil(values_ofinterest(v)),:) .* (decimalval)],1);
    end
    allcolors{v} = thiscolor;
end

figure;
set(gcf,'Position',[812 165 1027 805])
set(gcf,'Color',[1 1 1]);%set(gcf,'Color',[.9 .9 .9])
set(gca,'Color',[1 1 1]);%set(gca,'Color',[.9 .9 .9])

categoryIdx = repmat([1:length(values_ofinterest)]',1,length(subnames));

plotSpread(thicknesses','categoryIdx',categoryIdx,'categoryColors',allcolors,'MarkerSize',40,'distributionMarkers','o','binWidth',.3)
h = gca;
for i = 1:length(h.Children)
    h.Children(i).MarkerFaceColor = h.Children(i).Color;
    h.Children(i).Color = [0 0 0];
    h.Children(i).LineWidth = 2.5;
end

set(gca,'FontSize',30)
set(gca,'XTickLabels',[])


disp('SCAN vs foot')
[~,P,~,STATS] = ttest(thicknesses(:,1),thicknesses(:,2))

disp('SCAN vs hand')
[~,P,~,STATS] = ttest(thicknesses(:,1),thicknesses(:,3))

disp('SCAN vs mouth')
[~,P,~,STATS] = ttest(thicknesses(:,1),thicknesses(:,4))
        
