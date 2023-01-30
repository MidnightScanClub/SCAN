
subnames = {'SIC01','SIC02','SIC03','ME01','ME02','ME03','ME04'};%


subnetwork_vals = [1.5 17 10 11];% SCAN Foot Hand Face
othernet_vals = [1 2 3 5 6 7 8 9 12 15 16]; %DMN Visual FPN DAN Premot Lang Sal CON Aud PMN CAN

FCvals_toCON = zeros(length(subnames),length(subnetwork_vals));

FCvals_withothernets = zeros(length(subnetwork_vals),length(othernet_vals),length(subnames));

for subnum = 1:length(subnames)
    
    subname = subnames{subnum};
    disp(subname)
    
    if strcmp(subname(1:3),'SIC')
       
        basedir = '/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/';
        if strcmp(subname,'SIC01')
            
            tmasks = smartload([basedir subname '/onoff_tmask.mat']);%
            
            scanlist = textread([basedir subname '/' subname '_cast_onoff.txt'],'%s');%
            
            fslrdir = ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/' subname '/7112b_fs_LR_old/fsaverage_LR32k/'];
            
            aparc_L_file = [fslrdir 'MSC02.L.aparc.32k_fs_LR.label.gii'];
            aparc_R_file = [fslrdir 'MSC02.R.aparc.32k_fs_LR.label.gii'];
            
        else
            
            tmasks = smartload([basedir subname '/tmask.mat']);%
            
            scanlist = textread([basedir subname '/cast_scans.txt'],'%s');%
            
            fslrdir = ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/' subname '/7112b_fs_LR/fsaverage_LR32k/'];
            
            aparc_L_file = [fslrdir subname '.L.aparc.32k_fs_LR.label.gii'];
            aparc_R_file = [fslrdir subname '.R.aparc.32k_fs_LR.label.gii'];
            
        end
        
        infomapdir = ['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' subname '_precast_infomap_wacky2_subcortreg_ignoreverts/'];
        cd(infomapdir)
        motorspots = ft_read_cifti_mod([infomapdir subname '_rawassn_minsize10_regularized_networksplus_motorrestricted_CS.dscalar.nii']);%['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' subname '_precast_infomap_wacky2_subcortreg_ignoreverts/' subname '_spots_effectors_templatematch.dtseries.nii']);%'_spots_effectors_only_cort_templatematch_subcort_wta.dtseries.nii']);%
        motorspots.data(59413:end,:) = 0;
        %networks = ft_read_cifti_mod([infomapdir subname '_rawassn_minsize10_regularized_recolored.dscalar.nii']);
    

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
        
        
        basedir = ['/data/nil-bluearc/GMT/Evan/subjects/' subname '/'];
        infomapdir = [basedir 'infomap/REST_adaptive_moreverts_s1p7_subcortregressed/'];
        cd([infomapdir])
        
        data = ft_read_cifti_mod([basedir '/func/rest/ConcatenatedCiftis/Rest_OCME+MEICA+MGTR_s1.7_MotionCensored+Concatenated.dtseries.nii']);
        
        tmask_concat = smartload([basedir 'func/rest/tmasks/Tmask_' subname '.mat']);
        tmask_concat = logical(tmask_concat);
        sessions = smartload([basedir 'func/rest/tmasks/ScanIdx_' subname '.mat']);
        
        motorspots = ft_read_cifti_mod([infomapdir subname '_rawassn_minsize10_regularized_networksplus_oneID_CS.dscalar.nii']);%[infomapdir subname '_rawassn_minsize10_regularized_networksandmotorspots.dtseries.nii']);
        motorspots.data(59413:end,:) = 0;
        
        
        
        aparc_L_file = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.L.aparc.32k_fs_LR.label.gii'];
        aparc_R_file = [basedir 'anat/MNINonLinear/fsaverage_LR32k/' subname '.R.aparc.32k_fs_LR.label.gii'];
        
    end
    
    motorspots.data(abs(motorspots.data-1.5)<.01) = 1.5;
    motorspots.data(abs(motorspots.data-2.5)<.01) = 1.5;
    motorspots.data(abs(motorspots.data-6.6)<.01) = 1.5;
    
    
    
    
        %% spots as separate subnetworks
        
    othernettcs = zeros(size(data.data,2),length(othernet_vals));
    
    subnetworktcs = zeros(size(data.data,2),length(subnetwork_vals));
    for v = 1:length(subnetwork_vals)
        subnetworktcs(:,v) = mean(data.data(abs(motorspots.data - subnetwork_vals(v))<.01,:),1);
    end
    
    for v = 1:length(othernet_vals)
        othernettcs(:,v) = mean(data.data(abs(motorspots.data - othernet_vals(v))<.01,:),1);
    end
    
    corrmat = FisherTransform(paircorr_mod(subnetworktcs));
    
    FCvals_toCON(subnum,:) = FisherTransform(paircorr_mod(othernettcs(:,8),subnetworktcs));
    
    FCvals_withothernets(:,:,subnum) = FisherTransform(paircorr_mod(subnetworktcs,othernettcs));
    
    
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

%%
figure;
set(gcf,'Position',[812 165 1027 805])
set(gcf,'Color',[1 1 1]);%set(gcf,'Color',[.9 .9 .9])
set(gca,'Color',[1 1 1]);%set(gca,'Color',[.9 .9 .9])



subsind = [1:7];

%plot(FCvals_toCON(subsind,:)','k','LineWidth',3)

%categoryIdx = repmat([1:length(subsind)]',1,size(FCvals_toCON,2));
    
%plotSpread(FCvals_toCON(1:numsubs,:),'distributionColors',allcolors,'MarkerSize',40,'categoryIdx',categoryIdx,'categoryMarkers',{'o','s','d'})
%plotSpread(FCvals_toCON(subsind,:),'distributionColors',allcolors,'MarkerSize',40,'distributionMarkers','^')

categoryIdx = repmat([1:size(FCvals_toCON,2)]',1,length(subsind));
plotSpread(FCvals_toCON(subsind,:)','categoryIdx',categoryIdx,'categoryColors',allcolors,'MarkerSize',40,'distributionMarkers','o')

h = gca;
for i = 1:length(h.Children)
    h.Children(i).MarkerFaceColor = h.Children(i).Color;
    h.Children(i).Color = [0 0 0];
    h.Children(i).LineWidth = 2.5;
end

set(gca,'FontSize',30)
set(gca,'XTickLabels',[])
box off

%%
disp('CMI vs foot')
[~,P,~,STATS] = ttest(FCvals_toCON(:,1),FCvals_toCON(:,2))

disp('CMI vs hand')
[~,P,~,STATS] = ttest(FCvals_toCON(:,1),FCvals_toCON(:,3))

disp('CMI vs mouth')
[~,P,~,STATS] = ttest(FCvals_toCON(:,1),FCvals_toCON(:,4))

disp('foot vs hand')
[~,P,~,STATS] = ttest(FCvals_toCON(:,2),FCvals_toCON(:,3))

disp('foot vs mouth')
[~,P,~,STATS] = ttest(FCvals_toCON(:,2),FCvals_toCON(:,4))

disp('hand vs mouth')
[~,P,~,STATS] = ttest(FCvals_toCON(:,3),FCvals_toCON(:,4))


%%
diffs = repmat(FCvals_withothernets(1,:,:),[3,1,1]) - FCvals_withothernets(2:4,:,:);
absdiff = abs(diffs);

FCvals_withothernets_CMIveffectors = zeros(size(FCvals_withothernets,2),size(FCvals_withothernets,3));
for s = 1:length(subnames)
    [mins,mini] = min(absdiff(:,:,s),[],1);
    signs_ofmins = zeros(size(mins));
    for i = 1:length(signs_ofmins)
        FCvals_withothernets_CMIveffectors(i,s) = diffs(mini(i),i,s);
    end
 end
FCvals_withothernets_CMIveffectors = FCvals_withothernets_CMIveffectors';
%FCvals_withothernets_CMIveffectors = squeeze(min(absdiff,[],1) .* signdiff)';

%FCvals_withothernets_CMIveffectors = squeeze(min(repmat(FCvals_withothernets(1,:,:),[3,1,1]) - FCvals_withothernets(2:4,:,:),[],1))';

meanFCvals_withothernets = nanmean(FCvals_withothernets_CMIveffectors,1);
stderFCvals_withothernets = nanstd(FCvals_withothernets_CMIveffectors,[],1) ./ sqrt(size(FCvals_withothernets_CMIveffectors,1));

make_network_bargraph(othernet_vals,meanFCvals_withothernets,stderFCvals_withothernets,true)

disp('networks vs zero')
pvec_netVzero = [];
for netstotest = 1:11
    [~,P,~,STATS] = ttest(FCvals_withothernets_CMIveffectors(:,netstotest),0,'tail','right');
    disp(['network ' num2str(othernet_vals(netstotest)) ': T=' num2str(STATS.tstat) '; p=' num2str(P)])
    pvec_netVzero(end+1) = P;
end

disp('CON vs others')
pvec_netVcon = [];
for netstotest = [1:7 9:11]
    [~,P,~,STATS] = ttest(FCvals_withothernets_CMIveffectors(:,8),FCvals_withothernets_CMIveffectors(:,netstotest));
    disp(['CON vs net' num2str(othernet_vals(netstotest)) ': T=' num2str(STATS.tstat) '; p=' num2str(P)])
    pvec_netVcon(end+1) = P;
end
