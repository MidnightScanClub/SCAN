subnames = {'SIC01','SIC02'};%


xdistance = 30;

IDs_toinclude = [1.5 17 10 11]; %SCAN Foot Hand Face

networks_tocheck = [9];%CON
networknames = {'CON'};

taskinds_toexclude = [7 14];


for s = 1:length(subnames)
    SICname = subnames{s};

motorsubnetworks = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' SICname '_precast_infomap_wacky2_subcortreg_ignoreverts/' SICname '_rawassn_minsize10_regularized_networksplus_motorrestricted_CS.dscalar.nii']);%(['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' SICname '_precast_infomap_wacky2_subcortreg_ignoreverts/' SICname '_rawassn_minsize10_regularized_networksplus_motor_CS_premot.dscalar.nii']);%['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' SICname '_precast_infomap_wacky2_subcortreg_ignoreverts/' SICname '_rawassn_minsize10_regularized_networksplus_motorrestricted_CS.dscalar.nii']);

data = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/CIMT/motortask/' SICname '_ZstatsMotorContrasts_MEonly_smooth6.0.dscalar.nii']);
data.data(:,taskinds_toexclude) = [];
data.mapname(taskinds_toexclude) = [];

%
motorsubnetworks.data(59413:end,:) = 0;
%make them the same color
motorsubnetworks.data(abs(motorsubnetworks.data-2.5)<.01) = 1.5;
motorsubnetworks.data(abs(motorsubnetworks.data-6.6)<.01) = 1.5;


motorsubnetworksL = motorsubnetworks; motorsubnetworksL.data(29697:end) = 0;
motorsubnetworksR = motorsubnetworks; motorsubnetworksR.data(1:29696) = 0;

    
    allclusters = [];
allclusterIDs = [];


for IDnum = 1:length(IDs_toinclude)
    
    ID = IDs_toinclude(IDnum);
    
        allclusters = [allclusters (abs(motorsubnetworksL.data - ID)<.01)];
        allclusters = [allclusters (abs(motorsubnetworksR.data - ID)<.01)];
        
        allclusterIDs = [allclusterIDs ID ID];
        
       
    
end

for n = 1:length(networks_tocheck)

        allclusters = [allclusters (abs(motorsubnetworks.data - networks_tocheck(n))<.01)];
        
        allclusterIDs = [allclusterIDs networks_tocheck(n) ];
end

allclusters = logical(allclusters);


clustervals = zeros(length(allclusterIDs),size(data.data,2));
for c = 1:length(allclusterIDs)
    clustervals(c,:) = mean(data.data(allclusters(:,c),:),1);
end

%%


%%
for n = 1:length(networks_tocheck)
    
    theseclustervals = zeros(0,size(clustervals,2));
    for IDnum = 1:length(IDs_toinclude)
        theseclustervals = [theseclustervals ; clustervals(allclusterIDs==IDs_toinclude(IDnum),:)];
    end
    theseclustervals = [theseclustervals ; clustervals(allclusterIDs==networks_tocheck(n),:)];
    
    allsubsclustervals(:,:,s) = theseclustervals;

[~,sorti] = sort(theseclustervals(end,:),'ascend');
clustervals_resorted = theseclustervals([end [1:(end-1)]],sorti);

figure;
set(gcf,'Position',[813 30 1102 805])
set(gcf,'Color',[1 1 1]);%set(gcf,'Color',[.9 .9 .9])
set(gca,'Color',[1 1 1]);%set(gca,'Color',[.9 .9 .9])

hold on
imagesc(clustervals_resorted',[-10 10])
colormap jet
title([SICname ': ' networknames{n}])


for i = 1:(size(theseclustervals,1)-1)
    [r,p] = paircorr_mod(theseclustervals(i,:)',theseclustervals(end,:)');
    disp([networknames{n} ' vs network ' num2str(i) ': r=' num2str(r) '; p=' num2str(p)])
end

data.mapname(sorti(end : -1 : 1))'

axis off
end

%data.mapname(sorti(end : -1 : 1))'
end

theseclustervals = mean(allsubsclustervals,3);
[~,sorti] = sort(theseclustervals(end,:),'ascend');
clustervals_resorted = theseclustervals([end [1:(end-1)]],sorti);

figure;
set(gcf,'Position',[813 30 1102 805])
set(gcf,'Color',[1 1 1]);%set(gcf,'Color',[.9 .9 .9])
set(gca,'Color',[1 1 1]);%set(gca,'Color',[.9 .9 .9])

hold on
imagesc(clustervals_resorted',[-7 7])
colormap jet
title(['Mean: ' networknames{n}])

disp([networknames{n} ': ' num2str(mean([paircorr_mod(theseclustervals(1,:)',theseclustervals(end,:)') paircorr_mod(theseclustervals(2,:)',theseclustervals(end,:)')]))])

data.mapname(sorti(end : -1 : 1))'

axis off

