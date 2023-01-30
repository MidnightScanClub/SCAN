subnames = {'SIC01','SIC02','SIC03','ME01','ME02','ME03','ME04'};%

DISTALvol = '/data/nil-bluearc/GMT/Evan/Atlases/DISTAL/DISTAL_Ewert_2017/all_histo_labels_lr_MNI_0.25mm.nii.gz';
MNIatlas = '/data/cn2/data28/fsl_64/fsl_5.0.6/data/standard/MNI152_T1_1mm_brain.nii.gz';
atlas = load_untouch_nii_2D(DISTALvol);

[atlasvals,atlasnames,~]=xlsread('/data/nil-bluearc/GMT/Evan/Atlases/DISTAL/DISTAL_Ewert_2017/files_names.xlsx');
atlasnames = atlasnames(2:end,6);
unique_atlasnames = unique(atlasnames);
unique_atlasnames(strcmp(unique_atlasnames,'z_none')) = [];


mapvals = zeros(length(subnames),length(unique_atlasnames),3);

for subnum = 1:length(subnames)
    
    subname = subnames{subnum};
    disp(subname)
    
    if strcmp(subname(1:3),'SIC')
       
        basedir = '/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/';
        
        
        infomapdir = ['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' subname '_precast_infomap_wacky2_subcortreg_ignoreverts/'];
        
        
        
    else
        
        
        basedir = ['/data/nil-bluearc/GMT/Evan/subjects/' subname '/'];
        infomapdir = [basedir 'infomap/REST_adaptive_moreverts_s1p7_subcortregressed/'];
        
        
    end
    
    
     filename = [infomapdir subname '_Allspots_and_Effectors_FC.dtseries.nii'];%'_Allspots_v_effectors_mindiff.dtseries.nii'];%'MidCMI_corrmap.dtseries.nii'];
     dotinds = strfind(filename, '.');
     filename_short = filename(1:(dotinds(end-1)-1));
    
     system(['wb_command -cifti-separate ' filename ' COLUMN -volume THALAMUS_LEFT ' filename_short '_LThal.nii.gz -volume THALAMUS_RIGHT ' filename_short '_RThal.nii.gz'])
     
     system(['fslmaths ' filename_short '_RThal.nii.gz -add ' filename_short '_LThal.nii.gz ' filename_short '_LRThal.nii.gz'])
     
     %system(['fslselectvols -i ' filename_short '_LRThal.nii.gz -o Allspots_FC_LRThal.nii.gz --vols=0'])
     
    if strcmp(subname(1:3),'SIC')
        system(['niftigz_4dfp -4 ' filename_short '_LRThal.nii.gz ' filename_short '_LRThal'])
        system(['t4img_4dfp /data/nil-bluearc/raichle/lin64-tools/711-2B_to_MNI152lin_T1_t4 ' filename_short '_LRThal.4dfp.img ' filename_short '_LRThal_MNI_DISTALspace.4dfp.img -O' DISTALvol(1:end-7) '.4dfp.img'])
        
        delete([filename_short '_LRThal_MNI_DISTALspace.nii.gz'])
        system(['niftigz_4dfp -n ' filename_short '_LRThal_MNI_DISTALspace.4dfp.img ' filename_short '_LRThal_MNI_DISTALspace'])
    else
        
        system(['applywarp -i ' filename_short '_LRThal.nii.gz -r ' MNIatlas ' -o ' filename_short '_LRThal_MNI.nii.gz -w temp4flatmap_anatomical_bet_warpfield_2MNI.nii.gz']);
        system(['niftigz_4dfp -4 ' filename_short '_LRThal_MNI.nii.gz ' filename_short '_LRThal_MNI'])
        system(['t4img_4dfp none ' filename_short '_LRThal_MNI.4dfp.img ' filename_short '_LRThal_MNI_DISTALspace.4dfp.img -O' DISTALvol(1:end-7) '.4dfp.img'])
        
        
        %system(['t4img_4dfp none ' filename_short '_LRThal.4dfp.img ' filename_short '_LRThal_MNI_DISTALspace.4dfp.img -O' DISTALvol(1:end-7) '.4dfp.img'])
    
    delete([filename_short '_LRThal_MNI_DISTALspace.nii.gz'])
    system(['niftigz_4dfp -n ' filename_short '_LRThal_MNI_DISTALspace.4dfp.img ' filename_short '_LRThal_MNI_DISTALspace'])
    
    end
    data = load_untouch_nii_2D([filename_short '_LRThal_MNI_DISTALspace.nii.gz']);
    

    

    for n = 1:length(unique_atlasnames)
        
        map = false(size(atlas.img));
        thisstruct_atlasvals = atlasvals(strcmp(atlasnames,unique_atlasnames{n}));
        
        for v = 1:length(thisstruct_atlasvals)
            map(atlas.img==thisstruct_atlasvals(v)) = true;
        end
        
        for m = 1:size(data.img,2)
            mapvals(subnum,n,m) = mean(data.img(map,m));
        end
        
    end
        



end

%% SCAN
ps = zeros(length(unique_atlasnames),1);
ts = ps;
means = ps;
for i = 1:length(unique_atlasnames)
    [H,P,CI,STATS] = ttest(mapvals(:,i,1));
    ps(i) = P;
    ts(i) = STATS.tstat;
    means(i) = mean(mapvals(:,i,1));
end

p_fdr = FDR(ps,.05);
disp(['FDR-corrected p = ' num2str(p_fdr)])

for i = 1:length(unique_atlasnames)
    
    if (ps(i) < p_fdr) && (means(i)>0)
    
    disp([unique_atlasnames{i} ':'])
    disp(['minval: ' num2str(min(mapvals(:,i,1)))])
    disp(['mean: ' num2str(means(i))])
    disp(['ts: ' num2str(ts(i))])
    disp(['ps: ' num2str(ps(i))])%disp(['ps: ' num2str(ps(i,2:end)<p_fdr)])
    disp(' ')
    
    end
end

%%
ps = zeros(length(unique_atlasnames),size(mapvals,3)-1);
ts = ps;
for m = 1:(size(mapvals,3)-1)
    for i = 1:length(unique_atlasnames)
        [H,P,CI,STATS] = ttest(mapvals(:,i,1),mapvals(:,i,m+1),'tail','both');
        ps(i,m) = P;
        ts(i,m) = STATS.tstat;
    end
end

p_fdr = FDR(ps,.05);
disp(['FDR-corrected p = ' num2str(p_fdr)])

for i = 1:length(unique_atlasnames)
    
    if all(ps(i,2:end) < .05) && (mean(mapvals(:,i,1))>0)
    
    disp([unique_atlasnames{i} ':'])
    disp(['ts: ' num2str(ts(i,:))])
    disp(['ps: ' num2str(ps(i,:))])%disp(['ps: ' num2str(ps(i,2:end)<p_fdr)])
    disp(' ')
    
    end
end



%%
ps = zeros(length(unique_atlasnames),1);
ts = ps;
for i = 1:length(unique_atlasnames)
    [H,P,CI,STATS] = ttest(mapvals(:,i,1),mean(mapvals(:,i,[2:end]),3),'tail','both');
    ps(i) = P;
    ts(i) = STATS.tstat;
    
end

%p_fdr = FDR(ps,.05);
%disp(['FDR-corrected p = ' num2str(p_fdr)])
p_bonf = .05 ./ length(unique_atlasnames);
disp(['Bonf-corrected p = ' num2str(p_bonf)])

significant_structures = [];

for i = 1:length(unique_atlasnames)   
    
    if (ps(i)<p_bonf) && (mean(mapvals(:,i,1))>0)
    
    disp([unique_atlasnames{i} ': FC value = ' num2str(mean(mapvals(:,i,1)))])
    disp(['ts: ' num2str(ts(i))])
    disp(['ps: ' num2str(ps(i))])% '; significant and positive FC : ' num2str()])
    
    disp(' ')
    
    significant_structures = [significant_structures i];
    
    end
end


    




%%

power_surf_colormap = [1 0 0;0 0 .8;1 1 0;1 .8 .6;0 1 0;1 .6 1;0 .6 .6;0 0 0;.35 0 .65;.2 1 1;1 .5 0;.65 .25 1;0 .25 .6;.6 1 .6;.2 .3 1;1 1 1;0 .4 0; repmat([.25 .25 .25],50,1)];

subnetwork_vals = [1.5 17 10 11];


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





for s = significant_structures
    
    figure;
    set(gcf,'Position',[812 165 1027 805])
    set(gcf,'Color',[1 1 1]);%set(gcf,'Color',[.9 .9 .9])
    set(gca,'Color',[1 1 1]);%set(gca,'Color',[.9 .9 .9])
    
    
    subsind = [2 1 3:7];
    
    
    
    categoryIdx = repmat(1:4,1,length(subsind));
    plotSpread(squeeze(mapvals(subsind,s,:))','categoryIdx',categoryIdx,'categoryColors',allcolors,'MarkerSize',40,'distributionMarkers','o')
    
    h = gca;
    for i = 1:length(h.Children)
        h.Children(i).MarkerFaceColor = h.Children(i).Color;
        h.Children(i).Color = [0 0 0];
        h.Children(i).LineWidth = 2.5;
    end
    h.Children = h.Children(end:-1:1);
    
    set(gca,'FontSize',30)
    set(gca,'XTickLabels',[])
    box off
    disp(unique_atlasnames{s})
    %title(unique_atlasnames{s})
end