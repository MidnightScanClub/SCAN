function Project_cifti_to_cerebellum_flat(orig_filename, anatomical_inciftispace,discrete)
%Project_cifti_to_cerebellum_flat(orig_filename, anatomical_inciftispace,[discrete])


warning off
dotinds = strfind(orig_filename,'.');
basename = [orig_filename(1:(dotinds(end-1)-1)) '_temp4flatmap'];
filename = orig_filename;


%%

MNIatlas = '/data/cn2/data28/fsl_64/fsl_5.0.6/data/standard/MNI152_T1_1mm_brain.nii.gz';


system(['wb_command -cifti-separate ' filename ' COLUMN -volume-all ' basename '.nii.gz'])

if exist('anatomical_inciftispace','var') && ~isempty(anatomical_inciftispace)
    
    system(['niftigz_4dfp -4 ' basename ' ' basename])
    system(['t4img_4dfp none ' basename ' ' basename '_111 -n -O111']);
    system(['niftigz_4dfp -n ' basename '_111 ' basename '_111'])
    
    bet_anatomical = 'temp4flatmap_anatomical_bet.nii.gz';
    orig_anatomical = 'temp4flatmap_anatomical.nii.gz';
    copyfile(anatomical_inciftispace,orig_anatomical)
    
    
    if ~exist([bet_anatomical(1:end-7) '_warpfield_2MNI.nii.gz'])
        system(['bet ' orig_anatomical ' ' bet_anatomical]);
        system(['fnirt --ref=' MNIatlas ' --in=' bet_anatomical ' --iout=' bet_anatomical(1:end-7) '_nlwarpMNI.nii.gz --fout=' bet_anatomical(1:end-7) '_warpfield_2MNI.nii.gz']);
    end
    
    if exist('discrete','var') && logical(discrete)
        system(['applywarp -i ' basename '_111.nii.gz -r ' MNIatlas ' -o ' basename '_MNI.nii.gz -w ' bet_anatomical(1:end-7) '_warpfield_2MNI.nii.gz --interp=nn'])
    else
        system(['applywarp -i ' basename '_111.nii.gz -r ' MNIatlas ' -o ' basename '_MNI.nii.gz -w ' bet_anatomical(1:end-7) '_warpfield_2MNI.nii.gz']);% --interp=nn'])
    end
    
    delete(bet_anatomical);
    
else
    
    disp('No anatomical provided. Assuming cifti data is already in MNI space.')
    movefile([basename '.nii.gz'],[basename '_MNI.nii.gz'])
    
end
gunzip([basename '_MNI.nii.gz'])

%%

addpath /data/nil-bluearc/GMT/Scott/MSC_Subcortical/spm12/
addpath /data/nil-bluearc/GMT/Evan/Scripts/suit/

datasize = get_nifti_datasize([basename '_MNI.nii']);
if datasize(4)>1
    system(['fslsplit ' basename '_MNI.nii ' basename '_MNI -t'])
    for i = 1:datasize(4)
       gunzip([basename '_MNI' sprintf('%04i',(i-1)) '.nii.gz'])
       
       if exist('discrete') && logical(discrete)
           
           Data(:,i) = suit_map2surf([basename '_MNI' sprintf('%04i',(i-1)) '.nii'],'space','FSL','stats',@mode); 
           
       else
       
            Data(:,i) = suit_map2surf([basename '_MNI' sprintf('%04i',(i-1)) '.nii'],'space','FSL','stats',@nanmean); 
       
       end
       
       delete([basename '_MNI' sprintf('%04i',(i-1)) '.nii'])
       delete([basename '_MNI' sprintf('%04i',(i-1)) '.nii.gz'])
    end
else
    if exist('discrete') && logical(discrete)
        Data = suit_map2surf([basename '_MNI.nii'],'space','FSL','stats',@mode);
    else
        Data = suit_map2surf([basename '_MNI.nii'],'space','FSL','stats',@nanmean);
    end
end

save(gifti(single(Data)),[orig_filename(1:(dotinds(end-1)-1)) '_cerebellumflat.func.gii'])
system(['wb_command -set-structure ' orig_filename(1:(dotinds(end-1)-1)) '_cerebellumflat.func.gii CEREBELLUM'])

system(['wb_command -add-to-spec-file ' orig_filename(1:(dotinds(end-1)-1)) '_cerebellumflat.spec CEREBELLUM /data/nil-bluearc/GMT/Evan/Scripts/suit/flatmap/FLAT.surf.gii'])
system(['wb_command -add-to-spec-file ' orig_filename(1:(dotinds(end-1)-1)) '_cerebellumflat.spec CEREBELLUM /data/nil-bluearc/GMT/Evan/Scripts/suit/flatmap/SUIT.shape.gii'])
system(['wb_command -add-to-spec-file ' orig_filename(1:(dotinds(end-1)-1)) '_cerebellumflat.spec CEREBELLUM /data/nil-bluearc/GMT/Evan/Scripts/suit/flatmap/fissures.border'])
system(['wb_command -add-to-spec-file ' orig_filename(1:(dotinds(end-1)-1)) '_cerebellumflat.spec CEREBELLUM ' orig_filename(1:(dotinds(end-1)-1)) '_cerebellumflat.func.gii'])

%%
delete([basename '*'])


