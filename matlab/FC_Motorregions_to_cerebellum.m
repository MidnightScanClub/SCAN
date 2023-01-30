
subnames = {'SIC01','SIC02','SIC03','ME01','ME02','ME03','ME04'};%};%,

ROIs = gifti('/data/nil-bluearc/GMT/Evan/Atlases/cerebellar_atlases-master/Diedrichsen_2009/atl-Anatom_dseg.label.gii');
names = ROIs.labels.name';
names(end+1) = {'Right_X'};
names = names(2:end);




mindiff_vals = zeros(length(subnames),length(names),4);

motornames = {'SCAN','Foot','Hand','Mouth'};

for subnum = 1:length(subnames)
    
    subname = subnames{subnum};
    disp(subname)
    
    if strcmp(subname(1:3),'SIC')
       
        infomapdir = ['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' subname '_precast_infomap_wacky2_subcortreg_ignoreverts/'];
        
        
        
    else
        
        basedir = ['/data/nil-bluearc/GMT/Evan/subjects/' subname '/'];
        infomapdir = [basedir 'infomap/REST_adaptive_moreverts_s1p7_subcortregressed/'];
        
        
    end
    
    mindiff = gifti([infomapdir subname '_Allspots_and_Effectors_FC_cerebellumflat.func.gii']);
    
    for i = 1:length(names)
        mindiff_vals(subnum,i,:) = nanmean(mindiff.cdata(ROIs.cdata==i,:),1);
        
    end
    
end

ps = [];

for i = 1:length(names)
    
    for j = 2:4
        [H,P,CI,STATS] = ttest(mindiff_vals(:,i,1),mindiff_vals(:,i,j));%,'tail','right');
    
        disp([names{i} ', SCAN vs ' motornames{j} ': t=' num2str(STATS.tstat) '; p=' num2str(P)])
        
        ps(end+1) = P;
    end
    

    
    
end
ps(isnan(ps)) = [];

p_fdr = FDR(ps,.05)
        
