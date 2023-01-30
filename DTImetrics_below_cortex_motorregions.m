subnames = {'SIC01','SIC02','SIC03'};%};%,

metrics = {'FA','MD'};

network_vals = [10 11 17 1.5];
network_names = {'SMHand','SMMouth','SMFoot','SCAN'};
    

targetstruct_asegvals = [10 49];

brainstem_asegval = 16;

scaledpaththreshpct = .001;

medial_masks = {'/data/nil-bluearc/GMT/Evan/Scripts/CoEScripts/Resources/cifti_masks/L.atlasroi.32k_fs_LR.shape.gii','/data/nil-bluearc/GMT/Evan/Scripts/CoEScripts/Resources/cifti_masks/R.atlasroi.32k_fs_LR.shape.gii'};


metricval_means = zeros(length(subnames),length(network_vals),length(metrics));

for subnum = 1:length(subnames)
    
    subject = subnames{subnum};
    
    
    disp(subject)
    
       
        basedir = '/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/';
        if strcmp(subject,'SIC01')
            fslrdir = ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/' subject '/7112b_fs_LR_old/fsaverage_LR32k/'];
            
            bedpostdir = ['/data/nil-bluearc/GMT/seidern/SIC/SIC01_2019/Concat_Rescan/ModFSL_f2/ModFSL_f2.bedpostX/'];
            
            atlasT1vol = '/data/nil-bluearc/GMT/seidern/SIC/SIC01_2019/Concat_Rescan/Probtrac/7112b/SIC01_mpr_debias_avgT_111_t88.nii.gz';
            
            atlas_to_diffusion_xfm = '/data/nil-bluearc/GMT/seidern/SIC/SIC01_2019/Concat_Rescan/Probtrac/7112b/7112b_to_diff.mat';
            diffusion_to_atlas_xfm = '/data/nil-bluearc/GMT/seidern/SIC/SIC01_2019/Concat_Rescan/Probtrac/7112b/diff_to_7112b.mat';
            
            surfpials = {[fslrdir 'MSC02.L.pial.32k_fs_LR.surf.gii'],[fslrdir 'MSC02.R.pial.32k_fs_LR.surf.gii']};
            surfwhites = {[fslrdir 'MSC02.L.white.32k_fs_LR.surf.gii'],[fslrdir 'MSC02.R.white.32k_fs_LR.surf.gii']};
            surfmids = {[fslrdir 'MSC02.L.midthickness.32k_fs_LR.surf.gii'],[fslrdir 'MSC02.R.midthickness.32k_fs_LR.surf.gii']};
            surfinterior1s = {[fslrdir 'MSC02.L.interior1.32k_fs_LR.surf.gii'],[fslrdir 'MSC02.R.interior1.32k_fs_LR.surf.gii']};
            surfinterior2s = {[fslrdir 'MSC02.L.interior2.32k_fs_LR.surf.gii'],[fslrdir 'MSC02.R.interior2.32k_fs_LR.surf.gii']};
            surfinterior3s = {[fslrdir 'MSC02.L.interior3.32k_fs_LR.surf.gii'],[fslrdir 'MSC02.R.interior3.32k_fs_LR.surf.gii']};
            
            make_interior_surfaces('MSC02',fslrdir,[0.5 1:5])
            
            asegfile = ['/data/nil-bluearc/GMT/Laumann/MSC/fs5.3_native/MSC02/mri/aseg.mgz'];
            anat2atlas_t4file = ['/data/nil-bluearc/GMT/Laumann/MSC/MSC02/T1/MSC02_mpr1T_to_TRIO_Y_NDC_t4'];
            
            for m = 1:length(metrics)
                metricfiles{m} = ['/data/nil-bluearc/GMT/seidern/SIC/SIC01_2019/Concat_Rescan/Probtrac/SIC01_2019_concat_dwi_' metrics{m} '.nii.gz'];
            end
            
        else
            fslrdir = ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/' subject '/7112b_fs_LR/fsaverage_LR32k/'];
            
            bedpostdir = ['/data/nil-bluearc/GMT/seidern/SIC/' subject '/ConcatPre_DWI/ModFSL_f2/ModFSL_f2.bedpostX/'];
            
            atlasT1vol = ['/data/nil-bluearc/GMT/seidern/SIC/' subject '/ConcatPre_DWI/Probtrac/7112b/' subject '_mpr_debias_avgT_111_t88.nii.gz'];
            
            atlas_to_diffusion_xfm = ['/data/nil-bluearc/GMT/seidern/SIC/' subject '/ConcatPre_DWI/Probtrac/7112b/7112b_to_diff.mat'];
            diffusion_to_atlas_xfm = ['/data/nil-bluearc/GMT/seidern/SIC/' subject '/ConcatPre_DWI/Probtrac/7112b/diff_to_7112b.mat'];
    
            surfpials = {[fslrdir subject '.L.pial.32k_fs_LR.surf.gii'],[fslrdir subject '.R.pial.32k_fs_LR.surf.gii']};
            surfwhites = {[fslrdir subject '.L.white.32k_fs_LR.surf.gii'],[fslrdir subject '.R.white.32k_fs_LR.surf.gii']};
            surfmids = {[fslrdir subject '.L.midthickness.32k_fs_LR.surf.gii'],[fslrdir subject '.R.midthickness.32k_fs_LR.surf.gii']};
            surfinterior1s = {[fslrdir subject '.L.interior1.32k_fs_LR.surf.gii'],[fslrdir subject '.R.interior1.32k_fs_LR.surf.gii']};
            surfinterior2s = {[fslrdir subject '.L.interior2.32k_fs_LR.surf.gii'],[fslrdir subject '.R.interior2.32k_fs_LR.surf.gii']};
            surfinterior3s = {[fslrdir subject '.L.interior3.32k_fs_LR.surf.gii'],[fslrdir subject '.R.interior3.32k_fs_LR.surf.gii']};
            
            make_interior_surfaces(subject,fslrdir,[0.5 1:5])
            
            if strcmp(subject,'SIC02')
                asegfile = ['/data/nil-bluearc/GMT/Laumann/MSC/fs5.3_native/MSC06/mri/aseg.mgz'];
                anat2atlas_t4file = ['/data/nil-bluearc/GMT/Laumann/MSC/MSC06/T1/MSC06_mpr1T_to_TRIO_Y_NDC_t4'];
                
                for m = 1:length(metrics)
                    metricfiles{m} = ['/data/nil-bluearc/GMT/seidern/SIC/SIC02/ConcatPre_DWI/Probtrac/SIC02_dti_' metrics{m} '.nii.gz'];
                end
                
            else
                asegfile = '/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/SIC03/fs5.3_native_default/SIC03/mri/aseg.mgz';
                anat2atlas_t4file = '/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/SIC03/T1/SIC03_mpr1T_to_TRIO_Y_NDC_t4';
                
                for m = 1:length(metrics)
                    metricfiles{m} = ['/data/nil-bluearc/GMT/seidern/SIC/SIC03/ConcatPre_DWI/SIC03_dwi_' metrics{m} '.nii.gz'];
                end
                
            end
            
        end
        
    
    
    
    
    outdir = ['/data/nil-bluearc/GMT/Evan/CIMT/tractography/'];
    mkdir(outdir)
    
    infomapdir = ['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' subject '_precast_infomap_wacky2_subcortreg_ignoreverts/'];
    subnetworksfile = [infomapdir subject '_rawassn_minsize10_regularized_networksplus_motorrestricted_CS.dscalar.nii'];
    
    mask_with_BA = true;
    structs_ofinterest = [5 6];%{ 'BA4a'  'BA4p'};
    
    BAfile = ['/data/nil-bluearc/GMT/Evan/CIMT/BAs/' subject '.LR.BA.32k_fs_LR.dtseries.nii'];
    
    %mask_with_aparc = true;
    %structs_ofinterest = {'precentral','paracentral' 'postcentral'};% { 'L_precentral' 'L_paracentral'  'R_precentral' 'R_paracentral' }}; %  'L_postcentral'  'R_postcentral'
%     aparc_L_file = [fslrdir subject '.L.aparc.32k_fs_LR.label.gii'];
%     aparc_R_file = [fslrdir subject '.R.aparc.32k_fs_LR.label.gii'];
    
    
    
    
    
    
    
    
    
    
    
    %% get clusters
    
    cd(outdir)
    
    subnetworks = ft_read_cifti_mod(subnetworksfile);
    subnetworks.data(abs(subnetworks.data-2.5)<.01) = 1.5;
    subnetworks.data(abs(subnetworks.data-6.6)<.01) = 1.5;
    subnetworks.data(abs(subnetworks.data-1.5)<.01) = 1.5;
    
    % restrict to cortex
    subnetworks.data(subnetworks.brainstructure(subnetworks.brainstructure>0)>2) = 0;
    
    
    surfL = gifti(surfmids{1});
    surfR = gifti(surfmids{2});
    
    surf_coordsLR = [surfL.vertices ; surfR.vertices];
    surf_coordsLR = surf_coordsLR(subnetworks.brainstructure(subnetworks.brainstructure<3)>0,:);
    
    if mask_with_BA %mask with BAs
        
        BA = ft_read_cifti_mod(BAfile);
        
        BA_LR_mask = false(size(BA.data));
        
        
        for s=1:length(structs_ofinterest)
            BA_LR_mask(BA.data==structs_ofinterest(s)) = true;
        end
        
        
        subnetworks.data(~BA_LR_mask) = 0;
    end %end mask
    
    for m = 1:length(metricfiles)
        metricfile_atlas = [outdir '/' subject '_' metrics{m} '_7112b.nii.gz'];
        
        %system(['flirt -in ' metricfiles{m} ' -applyxfm -init ' diffusion_to_atlas_xfm ' -ref ' atlasT1vol ' -out ' metricfile_atlas ' -interp spline']);
        
        hems = {'L','R'};
        for h = 1:length(hems)
            
            %system(['wb_command -volume-to-surface-mapping ' metricfile_atlas ' ' surfinterior2s{h} ' ' metricfile_atlas(1:end-7) '_' hems{h} '.func.gii -ribbon-constrained ' surfinterior3s{h} ' ' surfinterior1s{h}]);%surfwhites{h}])
            system(['wb_command -volume-to-surface-mapping ' metricfile_atlas ' ' surfinterior1s{h} ' ' metricfile_atlas(1:end-7) '_' hems{h} '.func.gii -ribbon-constrained ' surfinterior2s{h} ' ' surfwhites{h}])
            %system(['wb_command -volume-to-surface-mapping ' metricfile_atlas ' ' surfmids{h} ' ' metricfile_atlas(1:end-7) '_' hems{h} '.func.gii -ribbon-constrained ' surfwhites{h} ' ' surfpials{h}])
            
        end
        
        
        %system(['wb_command -cifti-create-dense-timeseries ' metricfile_atlas(1:end-7) '_under_withgap.dtseries.nii -left-metric ' metricfile_atlas(1:end-7) '_' hems{1} '.func.gii -roi-left ' medial_masks{1} ' -right-metric ' metricfile_atlas(1:end-7) '_' hems{2} '.func.gii -roi-right ' medial_masks{2}])
        
        %data = ft_read_cifti_mod([metricfile_atlas(1:end-7) '_under_withgap.dtseries.nii']);
        
         system(['wb_command -cifti-create-dense-timeseries ' metricfile_atlas(1:end-7) '_under.dtseries.nii -left-metric ' metricfile_atlas(1:end-7) '_' hems{1} '.func.gii -roi-left ' medial_masks{1} ' -right-metric ' metricfile_atlas(1:end-7) '_' hems{2} '.func.gii -roi-right ' medial_masks{2}])
         
         data = ft_read_cifti_mod([metricfile_atlas(1:end-7) '_under.dtseries.nii']);

%         system(['wb_command -cifti-create-dense-timeseries ' metricfile_atlas(1:end-7) '_incortex.dtseries.nii -left-metric ' metricfile_atlas(1:end-7) '_' hems{1} '.func.gii -roi-left ' medial_masks{1} ' -right-metric ' metricfile_atlas(1:end-7) '_' hems{2} '.func.gii -roi-right ' medial_masks{2}])
%         
%        data = ft_read_cifti_mod([metricfile_atlas(1:end-7) '_incortex.dtseries.nii']);
        
        for v = 1:length(network_vals)
        
            metricval_means(subnum,v,m) = mean(data.data(subnetworks.data==network_vals(v)));
            
        end
        
    end
    
    
    
end
%%

power_surf_colormap = [1 0 0;0 0 .8;1 1 0;1 .8 .6;0 1 0;1 .6 1;0 .6 .6;0 0 0;.35 0 .65;.2 1 1;1 .5 0;.65 .25 1;0 .25 .6;.6 1 .6;.2 .3 1;1 1 1;0 .4 0; repmat([.25 .25 .25],50,1)];
allcolors = cell(1,length(network_vals));

for v = 1:length(network_vals)
    
    decimalval = mod(network_vals(v),1);
    if decimalval==0
        thiscolor = power_surf_colormap(network_vals(v),:);
    else
        thiscolor = sum([power_surf_colormap(floor(network_vals(v)),:) .* (1-decimalval) ; power_surf_colormap(ceil(network_vals(v)),:) .* (decimalval)],1);
    end
    allcolors{v} = thiscolor;
end

for m = 1:length(metricfiles)

figure;
set(gcf,'Position',[812 165 566 805])
set(gcf,'Color',[1 1 1]);%set(gcf,'Color',[.9 .9 .9])
set(gca,'Color',[1 1 1]);%set(gca,'Color',[.9 .9 .9])

thismetricvals = metricval_means(:,:,m);

categoryIdx = repmat([1:size(thismetricvals,2)]',1,length(subnames));
plotSpread(thismetricvals','categoryIdx',categoryIdx,'categoryColors',allcolors,'MarkerSize',40,'distributionMarkers','o')

h = gca;
for i = 1:length(h.Children)
    h.Children(i).MarkerFaceColor = h.Children(i).Color;
    h.Children(i).Color = [0 0 0];
    h.Children(i).LineWidth = 2.5;
end
title(metrics{m})
set(gca,'FontSize',30)
set(gca,'XTickLabels',[])
box off

export_fig(gca,[metrics{m} '_under.png'])

end


