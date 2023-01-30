cd /data/nil-bluearc/GMT/Evan/Atlases/StimSites/Roux2020/

[NUMs,Types,~]=xlsread('Roux_2020_stim_sites.xlsx');

stimcoords = NUMs(:,2:end);

classification_column = 1;


%% surface

surfL = gifti('/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.L.pial.32k_fs_LR.surf.gii');
surfR = gifti('/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.R.pial.32k_fs_LR.surf.gii');

coordsLR = [surfL.vertices ; surfR.vertices];
cifti = ft_read_cifti_mod('/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Cifti_normalwall_vertexindices.dtseries.nii');

coords_cifti = coordsLR(cifti.brainstructure(cifti.brainstructure<3)>0,:);


[distances,closestind] = pdist2(coords_cifti,stimcoords,'euclidean','Smallest',1);

vals = zeros(length(closestind),1);
for i = 1:length(closestind)
    if strcmp(Types{i,classification_column},'Hand')
        vals(i) = 10;
    elseif strcmp(Types{i,classification_column},'Face')
        vals(i) = 11;
    elseif strcmp(Types{i,classification_column},'Leg')
        vals(i) = 17;
    end
end

cifti.data(:) = 0;
cifti.data(closestind) = vals;

ft_write_cifti_mod('Roux_2020_stimsites_pial',cifti)

cifti_to_border_v2('Roux_2020_stimsites_pial.dtseries.nii',0,1,'default')

