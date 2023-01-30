function Project_cifti_to_rendered_subcortical_struct(filename,structurename_tomap,smoothing_onsurf,aparcfile)
% Project_cifti_to_rendered_subcortical_struct(filename,structurename_tomap,[smoothing_onsurf],[aparcfile])
%
% Makes smooth 3D rendered surfaces out of one or more subcortical
% structures within a cifti file, projects values within that cifti
% file onto the resulting surfaces, and saves those projections as a cifti
% in the space of the rendered structures. Only voxels intersecting the
% rendered surface will be projected, so this method ignores voxel values
% in the structure's interior.  
% 
% For optimal visualization in Workbench, rendered surfaces are set to
% CORTEX_LEFT and CORTEX_RIGHT regardless of their actual structures.
%
%
% Inputs:
%
% filename - The cifti file to use. 
%
% structurename_tomap - A string indicating The subcortical structure
%  within the cifti to render (e.g., 'AMYGDALA', or 'THALAMUS'). 
% Any structure that exists within the cifti can be used, though the
%  script assumes the presence of a LEFT and RIGHT structure and will not
%  map structures that don't follow this convention. 
% Additionally, you may input 'STRIATUM' to combine and jointly render accumbens,
%  caudate, and putamen; 'BASALGANGLIA' or 'BG' to combine and jointly
%  render the striatal structures plus pallidum; or 'MEDIALTEMPORAL' or
%  'MTL' to combine and jointly render amygdala and hippocampus. 
% Finally, you may input other combinations of multiple structures in a
%  cell array (e.g., {'BG','THALAMUS','MTL'}). Each separate
%  cell will be rendered independently and combined into a single surface.
%
% smoothing_onsurf - An optional parameter specifying the smoothing of rendered
%  data on the subcortical surface, in mm sigma. Omit to use no smoothing.
%
% aparcfile - An optional filename specifying the freesurfer-derived aparc
%  divisions (in the same space as the cifti volumetric elements) which
%  will be used to make higher-resolution subcortical structures.
%  
% E. Gordon 4/29/2021

fsldir = '/usr/local/pkg/fsl/';
surficedir = '/data/nil-bluearc/GMT/Evan/Scripts/surfice_atlas-master/';
rorden_SPMscriptsdir = '/data/nil-bluearc/GMT/Evan/Scripts/spmScripts-master/';
spmdir = '/data/nil-bluearc/GMT/Scott/MSC_Subcortical/spm12/';

addpath(surficedir)
addpath(rorden_SPMscriptsdir);
addpath(spmdir)
identmat = [fsldir '/etc/flirtsch/ident.mat'];

warning off

if iscell(structurename_tomap)
    structures = structurename_tomap;
    temp_structurename_tomap = [];
    for i = 1:length(structurename_tomap)
        temp_structurename_tomap = [temp_structurename_tomap structurename_tomap{i} '_'];
    end
    structurename_tomap = temp_structurename_tomap(1:end-1);
    clear temp_structurename_tomap
else
    structures = {structurename_tomap};
    
end

hems = {'LEFT','RIGHT'};



colormap_uniform = repmat([.75 .75 .75],256,1);



dotinds = strfind(filename,'.');
basename = filename(1:(dotinds(end-1)-1));

data = ft_read_cifti_mod(filename);

uniformdata = data;
uniformdata.data = uniformdata.data(:,1);
uniformdata.dimord = 'pos_time';
uniformdata.data(:) = 1;
ft_write_cifti_mod(['temp_' structurename_tomap 'mapping'],uniformdata);

if exist('aparcfile','var') && ~isempty(aparcfile)
    aparc_structs = load_untouch_nii_2D(aparcfile);
    vals = unique(aparc_structs.img);
    for v = 1:length(vals)
        if nnz(aparc_structs.img==vals(v)) > 50000
            aparc_structs.img(aparc_structs.img==vals(v)) = 0;
        end
    end
    
    
    aparc_struct_coords = get_volume_coords(aparcfile);
    aparcstructnums = cell(1,length(hems));
    
    surface_smoothness = 3;
    
else
    surface_smoothness = 6;
end



for h = 1:length(hems)
    
    
    hemcoords_incifti = [];
    heminds_incifti = [];
        
    for s = 1:length(structures)
        
        this_structure = structures{s};
        switch this_structure
            case {'STRIATUM'}
                this_rendered_structure = {'ACCUMBENS','CAUDATE','PUTAMEN'};
            case {'BASALGANGLIA','BG'}
                this_rendered_structure = {'ACCUMBENS','CAUDATE','PUTAMEN','PALLIDUM'};
            case {'MEDIALTEMPORAL','MTL'}
                this_rendered_structure = {'AMYGDALA','HIPPOCAMPUS'};
            otherwise
                this_rendered_structure = {this_structure};
        end
        
        fslstr = ['fslmaths'];
        
        for s2 = 1:length(this_rendered_structure)
            
            structurename = [this_rendered_structure{s2} '_' hems{h}];
            
            if exist('aparcfile','var') && ~isempty(aparcfile)
                
                aparcstructnums{h} = zeros(1,length(structures));
                
                system(['wb_command -cifti-separate temp_' structurename_tomap 'mapping.dtseries.nii COLUMN -volume ' structurename ' temp.nii.gz'])
                ciftistruct = load_untouch_nii_2D('temp.nii.gz');
                ciftistruct_volcoords = get_volume_coords('temp.nii.gz');
                ciftistruct_volcoords = ciftistruct_volcoords(logical(ciftistruct.img),:);
                dmat = pdist2(ciftistruct_volcoords,aparc_struct_coords);
                [~,mini] = min(dmat,[],2);
                closestvals = aparc_structs.img(mini);
                aparcstructnums{h}(s) = mode(closestvals(closestvals>0));
                
                delete('temp.nii.gz')
                
                structinds = [];
                for i = 1:length(aparcstructnums{h})
                    structinds = [structinds ; find(abs(aparc_structs.img-aparcstructnums{h}(i))<.1)];
                end
                
                out = aparc_structs;
                out.img(:) = 0;
                out.img(structinds) = 1;
                save_untouch_nii_2D(out,['temp_' structurename_tomap 'mapping_' structurename '.nii.gz']);
                
            else
                
                
                system(['wb_command -cifti-separate temp_' structurename_tomap 'mapping.dtseries.nii COLUMN -volume ' structurename ' temp_' structurename_tomap 'mapping_' structurename '.nii.gz'])
                
            end
            
            fslstr = [fslstr ' temp_' structurename_tomap 'mapping_' structurename '.nii.gz -add '];
            
            brainstructurenumber = find(strcmp(data.brainstructurelabel,structurename));
            
            structcoords = data.pos(data.brainstructure==brainstructurenumber,:);
            hemcoords_incifti = [hemcoords_incifti ; structcoords];
            
            heminds_incifti = [heminds_incifti ; find(data.brainstructure(data.brainstructure>0)==brainstructurenumber)];
            
            
            
        end
        
        fslstr = [fslstr(1:end-5) ' -bin temp_' this_structure 'mapping_' hems{h}];
        system(fslstr)
        
        if exist('aparcfile','var') && ~isempty(aparcfile)
            movefile(['temp_' this_structure 'mapping_' hems{h} '.nii.gz'],['temp_' this_structure 'mapping_' hems{h} '_upsample.nii.gz'])
        
        else
        
            system(['flirt -in temp_' this_structure 'mapping_' hems{h} ' -ref temp_' this_structure 'mapping_' hems{h} ' -applyisoxfm 0.5 -init ' identmat ' -interp trilinear -out temp_' this_structure 'mapping_' hems{h} '_upsample.nii.gz'])
            system(['fslmaths temp_' this_structure 'mapping_' hems{h} '_upsample -thr 0.5 -bin temp_' this_structure 'mapping_' hems{h} '_upsample'])
        end
        
        voldata = load_untouch_nii_2D(['temp_' this_structure 'mapping_' hems{h} '_upsample.nii.gz']);
        if s==1
            combinedvoldata = voldata;
        end
        combinedvoldata.img(logical(voldata.img)) = s;
        
        
    end
    
    save_untouch_nii_2D(combinedvoldata, ['temp_' structurename_tomap 'mapping_' hems{h} '_upsample.nii.gz'])
    
    nii_nii2atlas(['temp_' structurename_tomap 'mapping_' hems{h} '_upsample.nii.gz'],surface_smoothness,['temp_' structurename_tomap 'mapping_' hems{h} '_upsample'],[colormap_uniform].*255)
    [faces, vertices, ~, ~] = readMz3(['temp_' structurename_tomap 'mapping_' hems{h} '_upsamplemerge.mz3']);
    
    
    surf.vertices = single(vertices);
    surf.faces = single(faces);
    surf.mat = diag(ones(4,1));
    save(gifti(surf),[basename '_' structurename_tomap '_' hems{h} '.surf.gii'])
    system(['wb_command -set-structure ' basename '_' structurename_tomap '_' hems{h} '.surf.gii CORTEX_' hems{h} ' -surface-type ANATOMICAL'])
    
    [~,mini] = pdist2(hemcoords_incifti,vertices,'euclidean','Smallest',1);
    
    metric = data.data(heminds_incifti(mini),:);
    save(gifti(single(metric)),[basename '_' structurename_tomap '_' hems{h} '.func.gii'])
    
    
    if exist('smoothing_onsurf','var') && smoothing_onsurf>0
        
        system(['wb_command -metric-smoothing ' basename '_' structurename_tomap '_' hems{h} '.surf.gii ' basename '_' structurename_tomap '_' hems{h} '.func.gii ' num2str(smoothing_onsurf) ' ' basename '_' structurename_tomap '_' hems{h} '_smooth' num2str(smoothing_onsurf) '.func.gii'])
    end
    
end

 delete(['temp_' structurename_tomap 'mapping*'])
 delete(['temp_' this_structure 'mapping*'])
 
    
