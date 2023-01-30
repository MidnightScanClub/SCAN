% For each subject, compute mean myelin in BA4 for each motor region


subnames = {'SIC01','SIC02','SIC03','ME01','ME02','ME03','ME04'};

values_ofinterest = [17 10 11 1.5]; %Foot Hand Face SCAN

myelindensity = zeros(length(subnames),length(values_ofinterest));

mask_with_BA = true;
precentral_structs_ofinterest = [5 6];
postcentral_structs_ofinterest = [4 1];

myelinsigma = '2.12330450072004760682';

mask_with_aparc = true;
structs_ofinterest_all = {{'precentral','paracentral' 'postcentral'} { 'L_precentral' 'L_paracentral'  'R_precentral' 'R_paracentral' 'L_postcentral'  'R_postcentral'}}; % 

hems = {'L','R'};

medial_masks = {'/data/nil-bluearc/GMT/Evan/Scripts/CoEScripts/Resources/cifti_masks/L.atlasroi.32k_fs_LR.shape.gii','/data/nil-bluearc/GMT/Evan/Scripts/CoEScripts/Resources/cifti_masks/R.atlasroi.32k_fs_LR.shape.gii'};

cd /data/nil-bluearc/GMT/Evan/CIMT/myelin/

for subnum = 1:length(subnames)
    
    subname = subnames{subnum};
    disp(subname)
    
    myelinfile = [subname '.MyelinMap.32k_fs_LR.dtseries.nii'];
        
    
    if strcmp(subname(1:3),'SIC')
       
        basedir = '/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/';
        
        fslrdir_native = [basedir subname '/7112b_fs_LR/Native/'];
        fslrdir_32k = [basedir subname '/7112b_fs_LR/fsaverage_LR32k/'];
        
        ribbonfile = [basedir subname '/7112b_fs_LR/Ribbon/ribbon.nii.gz'];
        
        for h = 1:length(hems)
            surfpials{h} = [fslrdir_native subname '.' hems{h} '.pial.native.surf.gii'];
            surfwhites{h} = [fslrdir_native subname '.' hems{h} '.white.native.surf.gii'];
            surfmids{h} = [fslrdir_native subname '.' hems{h} '.midthickness.native.surf.gii'];
            surfmids32k{h} = [fslrdir_32k subname '.' hems{h} '.midthickness.32k_fs_LR.surf.gii'];
            nativedefspheres{h} = [fslrdir_native subname '.' hems{h} '.sphere.reg.reg_LR.native.surf.gii'];
            outspheres{h} = [fslrdir_32k subname '.' hems{h} '.sphere.32k_fs_LR.surf.gii'];
            thicknesses{h} = [fslrdir_native subname '.' hems{h} '.thickness.native.shape.gii'];
        end
        
        if strcmp(subname,'SIC01')
            structs_ofinterest = structs_ofinterest_all{2};
            
            fslrdir = ['/data/nil-bluearc/GMT/Laumann/MSC/MSM_nativeresampled2_TYNDC/MSC02/fsaverage_LR32k//'];
            
            aparc_L_file = [fslrdir 'MSC02.L.aparc.32k_fs_LR.label.gii'];
            aparc_R_file = [fslrdir 'MSC02.R.aparc.32k_fs_LR.label.gii'];
            
            lsurffile = [fslrdir 'MSC02.L.midthickness.32k_fs_LR.surf.gii'];
            rsurffile = [fslrdir 'MSC02.R.midthickness.32k_fs_LR.surf.gii'];
            
            
            
            
        else
            
            structs_ofinterest = structs_ofinterest_all{1};
            fslrdir = ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/' subname '/7112b_fs_LR/fsaverage_LR32k/'];
            
            
            lsurffile = [fslrdir subname '.L.midthickness.32k_fs_LR.surf.gii'];
            rsurffile = [fslrdir subname '.R.midthickness.32k_fs_LR.surf.gii'];
            
            aparc_L_file = [fslrdir subname '.L.aparc.32k_fs_LR.label.gii'];
            aparc_R_file = [fslrdir subname '.R.aparc.32k_fs_LR.label.gii'];
            
            if strcmp(subname,'SIC02')
                ribbonfile = ['/data/nil-bluearc/GMT/Evan/CIMT/myelin/temp/ribbon.nii.gz'];
            end
            
           
        end
        
        
        
        
        %%
        infomapdir = ['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' subname '_precast_infomap_wacky2_subcortreg_ignoreverts/'];
        
        motorspots = ft_read_cifti_mod([infomapdir subname '_rawassn_minsize10_regularized_networksplus_motorrestricted_CS.dscalar.nii']);%['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' subname '_precast_infomap_wacky2_subcortreg_ignoreverts/' subname '_spots_effectors_templatematch.dtseries.nii']);%'_spots_effectors_only_cort_templatematch_subcort_wta.dtseries.nii']);%
        motorspots.data(59413:end,:) = 0;
    
        
        
        
    else
        
        structs_ofinterest = structs_ofinterest_all{2};
        basedir = ['/data/nil-bluearc/GMT/Evan/subjects/' subname '/'];
        infomapdir = [basedir 'infomap/REST_adaptive_moreverts_s1p7_subcortregressed/'];
        motorspots = ft_read_cifti_mod([infomapdir subname '_rawassn_minsize10_regularized_CONandmotor_oneID_CS.dtseries.nii']);%[infomapdir subname '_rawassn_minsize10_regularized_networksandmotorspots.dtseries.nii']);
        motorspots.data(59413:end,:) = 0;
        
        
        myelinfile = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.MyelinMap.32k_fs_LR.dscalar.nii'];
        
        aparc_L_file = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.L.aparc.32k_fs_LR.label.gii'];
        aparc_R_file = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.R.aparc.32k_fs_LR.label.gii'];
        
        lsurffile = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.L.midthickness.32k_fs_LR.surf.gii'];
        rsurffile = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.R.midthickness.32k_fs_LR.surf.gii'];
        
        surfpials = {[basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.L.pial.32k_fs_LR.surf.gii'],[basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.R.pial.32k_fs_LR.surf.gii']};
        surfwhites = {[basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.L.white.32k_fs_LR.surf.gii'],[basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.R.white.32k_fs_LR.surf.gii']};
        surfmids = {[basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.L.midthickness.32k_fs_LR.surf.gii'],[basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.R.midthickness.32k_fs_LR.surf.gii']};
        
        
    end
    
    
    
    myelin = ft_read_cifti_mod(myelinfile);
    myelin = myelin.data;
    
%     %make them the same color
    motorspots.data(abs(motorspots.data-1.5)<.01) = 1.5;
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
        myelindensity(subnum,v) = mean(myelin(abs(motorspots.data-values_ofinterest(v))<.01));
    end
    
    if subnum==1
        out = motorspots;
        out.data(:) = 0;
        out.dimord = 'pos_scalar';
    end
    out.data(1:59412,subnum) = 0;
    for v = 1:(length(values_ofinterest)-1)
    out.data(abs(motorspots.data-values_ofinterest(v))<.01,subnum) = myelin(abs(motorspots.data-values_ofinterest(v))<.01);
    end
    out.mapname{subnum} = subname;
    
end


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

subsind = [1:7];


categoryIdx = repmat([1:size(myelindensity,2)]',1,length(subsind));
plotSpread(myelindensity(subsind,:)','categoryIdx',categoryIdx,'categoryColors',allcolors,'MarkerSize',40,'distributionMarkers','o')


h = gca;
for i = 1:length(h.Children)
    h.Children(i).MarkerFaceColor = h.Children(i).Color;
    h.Children(i).Color = [0 0 0];
    h.Children(i).LineWidth = 2.5;
end

set(gca,'FontSize',30)
set(gca,'XTickLabels',[])
        
