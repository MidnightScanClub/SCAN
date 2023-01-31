function Batch_wb_image_capture_subcort(ciftifile,outname,varargin)
%Batch_wb_image_capture(ciftifile,outname,'param1','value1','param2','value2'...)
%
% Automatically capture images of volumetric cifti data without having to
% use the Workbench GUI. Requires wb_command in your linux path and the
% /data/cn/data1/scripts/CIFTI_RELATED/Resources/ folder (and subfolders)
% in your matlab path.
%
% Mandatory inputs:
%
%  ciftifile - data to capture. Should be a dtseries.nii or dscalar.nii file
%
%  outname - root name of output image(s). Output files will be named
%   [outname]_map[#]_[image_view].png
% 
%
% Optional Parameters:
%
%  'image_view' - Subcortical structures to be shown in captured images.
%   Omit for default view, or specify one of: 
%   'All' (default), 'Striatum_Thalamus', 'Thalamus', 'Hipp_Amyg_Cereb'
%
%  'cifti_inds' - a vector specifying which columns of the input cifti to
%   capture images of. Omit to capture images of all columns.
%
%  'palette' - color palette to use when mapping values in cifti. Omit to
%   use the default 'ROY-BIG-BL' palette (or the palette already saved in
%   the cifti file, if any), or specify one of:   
%   'PSYCH','PSYCH-NO-NONE','ROY-BIG','ROY-BIG-BL','Orange-Yellow','Gray_Interp_Positive','Gray_Interp','clear_brain','videen_style','fidl','raich4_clrmid','raich6_clrmid','HSB8_clrmid','RBGYR20','RBGYR20P','POS_NEG','red-yellow','blue-lightblue','FSL','power_surf','fsl_red','fsl_green','fsl_blue','fsl_yellow','JET256'
%
%  'disp_zeros' - specify whether to display zero values (set to 'true') or
%   not (omit or set to 'false'). 
%
%  'colorscaling' - a 1x4 vector specifying the scaling of the color palette.
%   The vector should be formatted as [POS_MIN POS_MAX NEG_MIN NEG_MAX],
%   where entries represent numeric values in the cifti data. Omit to use
%   auto-scaling. 
%
%  'colorthresholding' - a 1x4 vector specifying how to threshold values in
%   the cifti file. The vector should be formatted as 
%   [MIN_INSIDE_THRESHOLD  MAX_INSIDE_THRESHOLD  MIN_OUTSIDE_THRESHOLD  MAX_OUTSIDE_THRESHOLD].
%   Either the INSIDE or OUTSIDE values should be set to NaN, indicating to
%   use the type of thresholding with real (non-NaN) values. With INSIDE
%   thresholding, only values between MIN_INSIDE_THRESHOLD and
%   MAX_INSIDE_THRESHOLD will be displayed; with OUTSIDE thresholding, only
%   values lower than MIN_OUTSIDE_THRESHOLD and higher than
%   MAX_OUTSIDE_THRESHOLD will be displayed. Omit to display all values.
%
%  'volume' - path to the anatomical volume to use as an underlay. 
%   Omit to use the default 1x1x1 Talaraich volume.
%
%
%E.Gordon 08/17


origcapturefolder = '/data/cn/data1/scripts/CIFTI_RELATED/Resources/batch_image_capture/';

[path,~,~] = fileparts(outname);
if isempty(path)
    path = pwd;
end

capturefolder = [path '/temp_image_capture_files/'];
mkdir(capturefolder)
try copyfile([origcapturefolder '/*.scene'],capturefolder); catch; end


%set defaults
image_view = 'All';
cifti_inds = [];
palette = [];
disp_zeros = 'false';
colorscaling = [];
colorthresholding = [];
volume = '/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/TRIO_KY_NDC_111.nii.gz';


%parse inputs
attributes = {'image_view','cifti_inds','palette','disp_zeros','colorscaling','colorthresholding','volume'};
for i = [1:2:length(varargin)]
    attributeindex = strcmp(varargin{i},attributes);
    if ~any(attributeindex)
        error(['Unrecognized parameter: ' varargin{i}])
    else
        eval([attributes{attributeindex} ' = varargin{i+1};'])
    end
end

%copy files
try copyfile(volume,[capturefolder '/volume.nii.gz']); catch; end

%get indices
data = ft_read_cifti_mod(ciftifile);
if isempty(cifti_inds)
    cifti_inds = [1:size(data.data,2)];
end

out = data;
for ind = cifti_inds(:)'
    if length(cifti_inds)>1
        thismapname = ['_map' num2str(ind)];
    else
        thismapname = [];
    end
    out.data = data.data(:,ind);
    out.mapname = {'Column number'};
    out.dimord = 'scalar_pos';
    ft_write_cifti_mod([capturefolder 'cifti'],out)
    
    if (exist('palette') && ~isempty(palette)) || (exist('colorscaling') && ~isempty(colorscaling)) || (exist('colorthresholding') && ~isempty(colorthresholding)) || (exist('disp_zero') && ~isempty(disp_zero))
        wbstring = ['wb_command -cifti-palette ' capturefolder '/cifti.dscalar.nii'];
        
        if ~isempty(colorscaling)
            wbstring = [wbstring ' MODE_USER_SCALE ' capturefolder '/cifti2.dscalar.nii -pos-user ' num2str(colorscaling(1)) ' ' num2str(colorscaling(2)) ' -neg-user ' num2str(colorscaling(3)) ' ' num2str(colorscaling(4))];
        else
            wbstring = [wbstring ' MODE_AUTO_SCALE_PERCENTAGE ' capturefolder '/cifti2.dscalar.nii'];
        end
        if ~isempty(palette)
            wbstring = [wbstring ' -palette-name ' palette];
        end
        if ~isempty(colorthresholding)
            if all(isnan(colorthresholding([3:4])))
                wbstring = [wbstring ' -thresholding THRESHOLD_TYPE_NORMAL THRESHOLD_TEST_SHOW_INSIDE ' num2str(colorthresholding(1)) ' ' num2str(colorthresholding(2))];
            elseif all(isnan(colorthresholding([1:2])))
                wbstring = [wbstring ' -thresholding THRESHOLD_TYPE_NORMAL THRESHOLD_TEST_SHOW_OUTSIDE ' num2str(colorthresholding(3)) ' ' num2str(colorthresholding(4))];
            else
                error('either inside or outside thresholds must have NaN values')
            end
        end
        
        wbstring = [wbstring ' -disp-zero ' disp_zeros];
        
        [~,~] = system(wbstring);
        try copyfile([capturefolder '/cifti2.dscalar.nii'],[capturefolder '/cifti.dscalar.nii']); catch; end
        delete([capturefolder '/cifti2.dscalar.nii']);
    end
    
    height = 925;
    width = 1401;
    
    [~,~] = system(['wb_command -show-scene ' capturefolder '/Capture_subcort.scene ' image_view ' ' outname thismapname '_' image_view '.png ' num2str(width) ' ' num2str(height)]);
end

rmdir(capturefolder,'s')





