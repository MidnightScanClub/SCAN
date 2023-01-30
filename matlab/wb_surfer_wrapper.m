function wb_surfer_wrapper(ciftidconnfile,borderfile,outname,varargin)
%wb_surfer_wrapper(ciftidconnfile,borderfile,outname,'param1','value1','param2','value2'...)
%
% Automatically create an .mp4 movie of FC values along a border file.
% Requires wb_command and wb_surfer in your linux path and the
% /data/cn/data1/scripts/CIFTI_RELATED/Resources/ folder (and subfolders) 
% in your matlab path.
%
% Mandatory inputs:
%
%  ciftidconnfile - FULL PATH of dconn file
%
%  outname - root name of output image(s). Output files will be named
%   [outname]_[image_view].mp4
% 
%
% Optional Parameters:
%
%  'image_view' - hemisphere (Left, Right) and view (Lateral, Medial,
%   Anterior, Posterior, Dorsal, Ventral) to be shown in captured images.
%   Omit for default view, or specify one of: 
%   'LR_LM' (default), 'LR_L','LR_M','L_LM','R_LM','LR_A','LR_P','LR_D','LR_V','L_L','L_R','R_L','R_M','L_A','L_P','L_D','L_V','R_A','R_P','R_D','R_V'
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
%   auto-scaling. Mutually exclusive with 'colorscaling_pct' and
%   'colorscaling_abspct'.
%
%  'colorscaling_pct' - a 1x4 vector specifying the scaling of the color palette.
%   The vector should be formatted as [POS_MIN POS_MAX NEG_MIN NEG_MAX],
%   where entries represent percentiles of values in the cifti data. Omit to use
%   auto-scaling. Mutually exclusive with 'colorscaling' and
%   'colorscaling_abspct'.
%
%  'colorscaling_abspct' - a 1x2 vector specifying the scaling of the color palette.
%   The vector should be formatted as [MIN MAX], where entries represent
%   percentiles of absolute values in the cifti data. Omit to use 
%   auto-scaling. Mutually exclusive with 'colorscaling' and
%   'colorscaling_pct'.
%
%  'colorthresholding' - a 1x4 vector specifing how to threshold values in
%   the cifti file. The vector should be formatted as 
%   [MIN_INSIDE_THRESHOLD  MAX_INSIDE_THRESHOLD  MIN_OUTSIDE_THRESHOLD  MAX_OUTSIDE_THRESHOLD].
%   Either the INSIDE or OUTSIDE values should be set to NaN, indicating to
%   use the type of thresholding with real (non-NaN) values. With INSIDE
%   thresholding, only values between MIN_INSIDE_THRESHOLD and
%   MAX_INSIDE_THRESHOLD will be displayed; with OUTSIDE thresholding, only
%   values lower than MIN_OUTSIDE_THRESHOLD and higher than
%   MAX_OUTSIDE_THRESHOLD will be displayed. Omit to display all values.
%
%  'Lsurf' - path to the left hemisphere surface to use for display. Omit
%   to use the default Conte69 fs_LR_32k inflated surface. 
%
%  'Rsurf' - path to the right hemisphere surface to use for display. Omit
%   to use the default Conte69 fs_LR_32k inflated surface. 
%
%  'Lborder' - path to the left hemisphere networks border file. 
%
%  'Rborder' - path to the right hemisphere networks border file. 
%
%  'Lsulc' - path to the left hemisphere sulcal depth map to display as an underlay
%   beneath the cifti data. Specify 'none' to use no sulcal depth underlay.
%   Omit to use the default Conte69 fs_LR_32k average sulcal depth map.
%
%  'Rsulc' - path to the right hemisphere sulcal depth map to display as an underlay
%   beneath the cifti data. Specify 'none' to use no sulcal depth underlay.
%   Omit to use the default Conte69 fs_LR_32k average sulcal depth map.
%
%  'volume' - path to the anatomical volume to use as an underlay. 
%   Omit to use the default 1x1x1 Talaraich volume.%
%
%E.Gordon 11/21


origcapturefolder = '/data/nil-bluearc/GMT/Evan/Scripts/wb_surfer_wrapper_files/';

[path,~,~] = fileparts(outname);
if isempty(path)
    path = pwd;
end

capturefolder = [path '/temp_image_capture_files_' num2str(randi(100000)) '/'];
mkdir(capturefolder)

system(['ln -s ' origcapturefolder '/Capture.scene ' capturefolder '/Capture.scene'])

%try copyfile([origcapturefolder '/*'],capturefolder); catch; end


%set defaults
image_view = 'LR_LM';
framerate = 7;
Lsurf = '/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.L.inflated.32k_fs_LR.surf.gii';
Rsurf = '/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.R.inflated.32k_fs_LR.surf.gii';
Lsulc = '/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.L.sulc.32k_fs_LR.shape.gii';
Rsulc = '/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.R.sulc.32k_fs_LR.shape.gii';
volume = '/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/TRIO_KY_NDC_111.nii.gz';

%parse inputs
attributes = {'image_view','cifti_inds','palette','disp_zeros','colorscaling','colorscaling_pct','colorscaling_abspct','colorthresholding','Lsurf','Rsurf','Lborder','Rborder','Lsulc','Rsulc','volume','framerate'};
for i = [1:2:length(varargin)]
    attributeindex = strcmp(varargin{i},attributes);
    if ~any(attributeindex)
        error(['Unrecognized parameter: ' varargin{i}])
    else
        eval([attributes{attributeindex} ' = varargin{i+1};'])
    end
end

%copy files
system(['ln -s ' Lsurf ' ' capturefolder '/L.surf.gii'])
system(['ln -s ' Rsurf ' ' capturefolder '/R.surf.gii'])
system(['ln -s ' volume ' ' capturefolder '/volume.nii.gz'])
if exist('Lborder','var')
    system(['ln -s ' Lborder ' ' capturefolder '/L.border'])
end
if exist('Rborder','var')
    system(['ln -s ' Rborder ' ' capturefolder '/R.border'])
end


if strcmp(Lsulc,'none')
    temp = gifti(Lsurf);
    save(gifti(single(zeros(size(temp.vertices,1),1))),[capturefolder '/L.sulc.shape.gii'])
else
    system(['ln -s ' Lsulc ' ' capturefolder '/L.sulc.shape.gii'])
end

if strcmp(Rsulc,'none')
    temp = gifti(Rsurf);
    save(gifti(single(zeros(size(temp.vertices,1),1))),[capturefolder '/R.sulc.shape.gii'])
else
    system(['ln -s ' Rsulc ' ' capturefolder '/R.sulc.shape.gii'])
end

[ciftipathpath,~,~] = fileparts(ciftidconnfile);
if isempty(ciftipathpath)
    ciftidconnfile = [pwd '/' ciftidconnfile];
end
system(['ln -s ' ciftidconnfile ' ' capturefolder '/Corr.dconn.nii']) %copyfile(ciftifile,[capturefolder '/Corr.dconn.nii'])


if exist('palette','var') || exist('colorscaling','var') || exist('colorscaling_pct','var') || exist('colorscaling_abspct','var') || exist('colorthresholding','var') || exist('disp_zero','var')
    
    capturelines = textread([capturefolder '/Capture.scene'],'%s','delimiter','\n','bufsize',16384000);
    outputfilename = [capturefolder '/Capture_temp.scene'];
    warning off 
    delete([outputfilename]);
    fid = fopen([outputfilename],'at'); %open the output file for writing
    fclose(fid);
    
    for l = 1:length(capturelines)
        linetext = capturelines{l};
        if exist('palette','var') && any(strfind(linetext,'PaletteName'))
            newlinetext = '&lt;PaletteName&gt;videen_style&lt;/PaletteName&gt;';%['&lt;PaletteName&gt;' palette ';/PaletteName&gt;'];
        elseif exist('colorscaling','var') && any(strfind(linetext,'ScaleMode'))
            newlinetext = ['&lt;ScaleMode&gt;MODE_USER_SCALE&lt;/ScaleMode&gt;'];
        elseif exist('colorscaling','var') && any(strfind(linetext,'UserScaleValues'))
            newlinetext = ['&lt;UserScaleValues&gt;' num2str(colorscaling(3)) ' ' num2str(colorscaling(4)) ' ' num2str(colorscaling(1)) ' ' num2str(colorscaling(2)) '&lt;/UserScaleValues&gt;'];
        elseif exist('colorscaling_pct','var') && any(strfind(linetext,'ScaleMode'))
            newlinetext = ['&lt;ScaleMode&gt;MODE_AUTO_SCALE_PERCENTAGE&lt;/ScaleMode&gt;'];
        elseif exist('colorscaling_pct','var') && any(strfind(linetext,'AutoScalePercentageValues'))
            newlinetext = ['&lt;AutoScalePercentageValues&gt;' num2str(colorscaling_pct(3)) ' ' num2str(colorscaling_pct(4)) ' ' num2str(colorscaling_pct(1)) ' ' num2str(colorscaling_pct(2)) '&lt;/AutoScalePercentageValues&gt;'];
        elseif exist('colorscaling_abspct','var') && any(strfind(linetext,'ScaleMode'))
            newlinetext = ['&lt;ScaleMode&gt;MODE_AUTO_SCALE_ABSOLUTE_PERCENTAGE&lt;/ScaleMode&gt;'];
        elseif exist('colorscaling_abspct','var') && any(strfind(linetext,'AutoScaleAbsolutePercentageValues'))
            newlinetext = ['&lt;AutoScaleAbsolutePercentageValues&gt;' num2str(colorscaling_abspct(1)) ' ' num2str(colorscaling_abspct(2)) '&lt;/AutoScaleAbsolutePercentageValues&gt;'];
        elseif exist('colorthresholding','var') && any(strfind(linetext,'ThresholdType'))
            newlinetext = ['&lt;ThresholdType&gt;THRESHOLD_TYPE_NORMAL&lt;/ThresholdType&gt;'];
        elseif exist('colorthresholding','var') && any(strfind(linetext,'ThresholdTest')) && all(isnan(colorthresholding([1:2])))
            newlinetext = ['&lt;ThresholdTest&gt;THRESHOLD_TEST_SHOW_OUTSIDE&lt;/ThresholdTest&gt;'];
        elseif exist('colorthresholding','var') && any(strfind(linetext,'ThresholdTest')) && all(isnan(colorthresholding([3:4])))
            newlinetext = ['&lt;ThresholdTest&gt;THRESHOLD_TEST_SHOW_INSIDE&lt;/ThresholdTest&gt;'];
        elseif exist('colorthresholding','var') && any(strfind(linetext,'ThresholdNormalValues'))
            thresholdvals = colorthresholding(~isnan(colorthresholding));
            newlinetext = ['&lt;ThresholdNormalValues&gt;' num2str(thresholdvals(1)) ' ' num2str(thresholdvals(2)) '&lt;/ThresholdNormalValues&gt;'];
        elseif exist('disp_zeros','var') && any(strfind(linetext,'DisplayZeroData'))    
            newlinetext = ['&lt;DisplayZeroData&gt;' disp_zeros '&lt;/DisplayZeroData&gt;'];
        else
            newlinetext = linetext;
        end
        dlmwrite(outputfilename,newlinetext,'delimiter','','-append')
    end
    delete([capturefolder '/Capture.scene'])
    movefile([capturefolder '/Capture_temp.scene'],[capturefolder '/Capture.scene'])
end


     height = 1488;
     width = 2250;
    if any(strcmp(image_view,{'LR_L','LR_M','L_LM','R_LM'}))
        height = ceil(height ./ 2);
    end

system(['wb_surfer -b ' borderfile ' -r ' num2str(framerate) ' ' capturefolder '/Capture.scene ' image_view ' ' outname '_' image_view '.mp4 -w ' num2str(width) ' ' num2str(height) ' -e libx264'])


rmdir(capturefolder,'s')





