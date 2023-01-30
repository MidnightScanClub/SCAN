

subnames = {'SIC01','SIC02','SIC03','ME01','ME02','ME03','ME04'};%

subnetwork_vals = [1.5 17 10 11]; %SCAN Foot Hand Face

divisions = {'DorsAnt','VentAnt','DorsPost','VentPost'};

local_regression = true;

mask_with_aparc = false;
structs_ofinterest_all = {{'precentral','paracentral' } { 'L_precentral' 'L_paracentral'  'R_precentral' 'R_paracentral' }}; % 'postcentral' 'L_postcentral'  'R_postcentral'

fcvals = zeros(length(subnames),length(subnetwork_vals),length(divisions));

for subnum = 1:length(subnames)
    
    subname = subnames{subnum};
    disp(subname)
    
    if strcmp(subname(1:3),'SIC')
       
        basedir = '/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/';
        structs_ofinterest = structs_ofinterest_all{1};
        if strcmp(subname,'SIC01')
            
            tmasks = smartload([basedir subname '/onoff_tmask.mat']);%
            
            scanlist = textread([basedir subname '/' subname '_cast_onoff.txt'],'%s');%
            
            fslrdir = ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/' subname '/7112b_fs_LR_old/fsaverage_LR32k/'];
            
            dmatname =  ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/SIC01/bold1_222/cifti_distances/SIC01_distmat_surf_geodesic_vol_euclidean_xhem1000_uint8.mat'];%
            
            aparc_L_file = [fslrdir 'MSC02.L.aparc.32k_fs_LR.label.gii'];
            aparc_R_file = [fslrdir 'MSC02.R.aparc.32k_fs_LR.label.gii'];
            
        else
            
            tmasks = smartload([basedir subname '/tmask.mat']);%
            
            scanlist = textread([basedir subname '/cast_scans.txt'],'%s');%
            
            fslrdir = ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/' subname '/7112b_fs_LR/fsaverage_LR32k/'];
            
            aparc_L_file = [fslrdir subname '.L.aparc.32k_fs_LR.label.gii'];
            aparc_R_file = [fslrdir subname '.R.aparc.32k_fs_LR.label.gii'];
            
            if strcmp(subname,'SIC02')
                dmatname =  ['/data/nil-bluearc/GMT/Laumann/MSC/MSM_nativeresampled2_TYNDC/MSC06/fsaverage_LR32k/MSC06_distmat_surf_geodesic_vol_euclidean_xhem1000_uint8.mat'];%['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/SIC01/bold1_222/cifti_distances/SIC01_distmat_surf_geodesic_vol_euclidean_xhem1000_uint8.mat'];%
            else
                dmatname = ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/' subname '/bold1_222/cifti_distances/' subname '_distmat_surf_geodesic_vol_euclidean_xhem1000_uint8.mat'];
                
            end
            
        end
        
        infomapdir = ['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' subname '_precast_infomap_wacky2_subcortreg_ignoreverts/'];
        cd(infomapdir)
        motorspots = ft_read_cifti_mod([infomapdir subname '_rawassn_minsize10_regularized_networksplus_motorrestricted_CS.dscalar.nii']);%['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' subname '_precast_infomap_wacky2_subcortreg_ignoreverts/' subname '_spots_effectors_templatematch.dtseries.nii']);%'_spots_effectors_only_cort_templatematch_subcort_wta.dtseries.nii']);%
        motorspots.data(59413:end,:) = 0;
        
        
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
        
        
        structs_ofinterest = structs_ofinterest_all{2};
        basedir = ['/data/nil-bluearc/GMT/Evan/subjects/' subname '/'];
        infomapdir = [basedir 'infomap/REST_adaptive_moreverts_s1p7_subcortregressed/'];
        dmatname = [basedir '/anat/MNINonLinear/fsaverage_LR32k/distances/normalwall_distmat_surf_geodesic_vol_euclidean_xhemlarge_uint8.mat'];
        cd([infomapdir])
        
        data = ft_read_cifti_mod([basedir '/func/rest/ConcatenatedCiftis/Rest_OCME+MEICA+MGTR_s1.7_MotionCensored+Concatenated.dtseries.nii']);
        
        tmask_concat = smartload([basedir 'func/rest/tmasks/Tmask_' subname '.mat']);
        tmask_concat = logical(tmask_concat);
        sessions = smartload([basedir 'func/rest/tmasks/ScanIdx_' subname '.mat']);
        %sessions = ones(size(sessions));
        
        motorspots = ft_read_cifti_mod([infomapdir subname '_rawassn_minsize10_regularized_CONandmotor_oneID_CS.dtseries.nii']);%[infomapdir subname '_rawassn_minsize10_regularized_networksandmotorspots.dtseries.nii']);
        motorspots.data(59413:end,:) = 0;
        
        aparc_L_file = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.L.aparc.32k_fs_LR.label.gii'];
        aparc_R_file = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.R.aparc.32k_fs_LR.label.gii'];
        
    end
    
    
    if local_regression
        dmat = smartload(dmatname);
        
        
        subcort_vox = find(data.brainstructure(data.brainstructure>0)>2);
        for i = 1:length(subcort_vox)
            neighbor_cort = (dmat(:,subcort_vox(i))<=20) & (data.brainstructure(data.brainstructure>0)<3);
            if any(neighbor_cort)
                meansignal = mean(data.data(neighbor_cort,:),1);
                [~,~,r] = regress(data.data(subcort_vox(i),:)',[meansignal' ones(length(meansignal),1)]);
                data.data(subcort_vox(i),:) = r;
            end
        end
    end
    
    
    
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
    
    
    motorspots.data(abs(motorspots.data-2.5)<.01) = 1.5;
    motorspots.data(abs(motorspots.data-6.6)<.01) = 1.5;
    
   %% spots as separate subnetworks
   
   putamen_inds = [find(data.brainstructure(data.brainstructure>0)==find(strcmp(data.brainstructurelabel,'PUTAMEN_LEFT'))) ; find(data.brainstructure(data.brainstructure>0)==find(strcmp(data.brainstructurelabel,'PUTAMEN_RIGHT')))];
   data_coords = data.pos(data.brainstructure>0,:);
   putamen_coords = data_coords(putamen_inds,:);
        
   putamen_ymiddle = mean([min(putamen_coords(:,2)) max(putamen_coords(:,2))]);
   putamen_zmiddle = mean([min(putamen_coords(:,3)) max(putamen_coords(:,3))]);
   
   putameninds = cell(4,1);
   putameninds{1} = putamen_inds((putamen_coords(:,2)>=putamen_ymiddle) & (putamen_coords(:,3)>=putamen_zmiddle));
   putameninds{2} = putamen_inds((putamen_coords(:,2)>=putamen_ymiddle) & (putamen_coords(:,3)<putamen_zmiddle));
   putameninds{3} = putamen_inds((putamen_coords(:,2)<putamen_ymiddle) & (putamen_coords(:,3)>=putamen_zmiddle));
   putameninds{4} = putamen_inds((putamen_coords(:,2)<putamen_ymiddle) & (putamen_coords(:,3)<putamen_zmiddle));
   
   for d = 1:length(divisions)
       for n = 1:length(subnetwork_vals)
           
           subnetwork_tc = mean(data.data(find(abs(motorspots.data(1:59412)-subnetwork_vals(n))<.01),:),1);
           division_tc = mean(data.data(putameninds{d},:),1);
           
           fcvals(subnum,n,d) = FisherTransform(paircorr_mod(subnetwork_tc',division_tc'));
           
       end
   end
           
   
    
end
%%
power_surf_colormap = [1 0 0;0 0 .8;1 1 0;1 .8 .6;0 1 0;1 .6 1;0 .6 .6;0 0 0;.35 0 .65;.2 1 1;1 .5 0;.65 .25 1;0 .25 .6;.6 1 .6;.2 .3 1;1 1 1;0 .4 0; repmat([.25 .25 .25],50,1)];

allcolors = cell(1,length(subnetwork_vals));

for v = 1:length(subnetwork_vals)
    
    decimalval = mod(subnetwork_vals(v),1);
    if decimalval==0
        thiscolor = power_surf_colormap(subnetwork_vals(v),:);
    else
        thiscolor = sum([power_surf_colormap(floor(subnetwork_vals(v)),:) .* (1-decimalval) ; power_surf_colormap(ceil(subnetwork_vals(v)),:) .* (decimalval)],1);
    end
    allcolors{v} = thiscolor;
end

for d = 1:length(divisions)
%     meanvals = nanmean(fcvals(:,:,d),1);
%     stderrvals = nanstd(fcvals(:,:,d),[],1) ./ sqrt(length(subnames));
%     
%     make_network_bargraph(subnetwork_vals,meanvals,stderrvals,false)
%     title(divisions{d})
    
    
    figure;
set(gcf,'Position',[812 165 1027 805])
set(gcf,'Color',[1 1 1]);%set(gcf,'Color',[.9 .9 .9])
set(gca,'Color',[1 1 1]);%set(gca,'Color',[.9 .9 .9])

division_fcvals = fcvals(:,:,d);

subsind = [1:7];

categoryIdx = repmat([1:size(division_fcvals,2)]',1,length(subsind));
plotSpread(division_fcvals(subsind,:)','categoryIdx',categoryIdx,'categoryColors',allcolors,'MarkerSize',40,'distributionMarkers','o')

h = gca;
for i = 1:length(h.Children)
    h.Children(i).MarkerFaceColor = h.Children(i).Color;
    h.Children(i).Color = [0 0 0];
    h.Children(i).LineWidth = 2.5;
end
title(divisions{d})
set(gca,'FontSize',30)
set(gca,'XTickLabels',[])
box off
 
    
disp([divisions{d} 'CMI vs foot'])
[~,P,~,STATS] = ttest(division_fcvals(:,1),division_fcvals(:,2))

disp([divisions{d} 'CMI vs hand'])
[~,P,~,STATS] = ttest(division_fcvals(:,1),division_fcvals(:,3))

disp([divisions{d} 'CMI vs mouth'])
[~,P,~,STATS] = ttest(division_fcvals(:,1),division_fcvals(:,4))
    
    
end