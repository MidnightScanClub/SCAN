subname = 'SIC02';

cd /data/nil-bluearc/GMT/Evan/CIMT/motortask/



winnerfile = [subname '_motor_allwinners_LRcombine_central_precentral.dtseries.nii'];%[subname '_motor_allwinners_LRcombine_CS.dtseries.nii'];%
if strcmp(subname,'SIC01')
    dmatname =  ['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/SIC01/bold1_222/cifti_distances/SIC01_distmat_surf_geodesic_vol_euclidean_xhem1000_uint8.mat'];%
end

if strcmp(subname,'SIC02')
    dmatname =  ['/data/nil-bluearc/GMT/Laumann/MSC/MSM_nativeresampled2_TYNDC/MSC06/fsaverage_LR32k/MSC06_distmat_surf_geodesic_vol_euclidean_xhem1000_uint8.mat'];%['/data/nil-bluearc/GMT/Dillan/preproc_2018-07-03/SIC01/bold1_222/cifti_distances/SIC01_distmat_surf_geodesic_vol_euclidean_xhem1000_uint8.mat'];%
end

hem = 'L';

activationfile = [subname '_ZstatsMotorContrasts_MEonly_smooth6.0.dscalar.nii'];

networks = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' subname '_precast_infomap_wacky2_subcortreg_ignoreverts/' subname '_rawassn_minsize10_regularized_motoronly_oneID_CS.dtseries.nii']);

Face = true;
LE = true;
UE = true;

mapnum_to_color_mapping_orig = [1 2 2 3 3 4 5 6 7 7 8 9 9 12 12 13 15 15 17 17 11 10 10];
fitparam_orig = [50 20 20 50 50 50 50 50 50 50 50 50 50 50 50 20 20 20 20 20 50 50 50];
colors_orig = {[12 63 216],[25 127 229],[25 127 229],[38 191 242],[38 191 242],[0 0 204],[225 155 0],[237 186 0],[0 51 0],[0 51 0],[0 255 0],[0 187 0],[0 187 0],[0 119 0],[0 119 0],[214 124 0],[202 94 0],[202 94 0],[191 63 0],[191 63 0],[248 216 0],[51 255 255],[51 255 255]};


make_eachplot = false;
export_plots = false;

dmat = smartload(dmatname);

all_F2s = cell(0,1);
all_F3s = cell(0,1);
all_colors = [];


if strcmp(hem,'L')
    mapinds = [7 8 16 18 20 21 10 11 13 15 23 6 5 3 1];
else
    mapinds = [7 8 16 17 19 21 9 11 12 14 22 6 4 2 1];
end


vertices_inline_data = ft_read_cifti_mod(vertices_inline_file);

winners = ft_read_cifti_mod(winnerfile);
winners.data = winners.data(:,1);

if strcmp(hem,'L')
    winners.data(29697:end) = 0;
    vertices_inline_data.data(29697:end) = 0;
elseif strcmp(hem,'R')
    winners.data(1:29696) = 0;
    vertices_inline_data.data(1:29696) = 0;
end

vertex_nums = unique(vertices_inline_data.data);
vertex_nums(vertex_nums==0) = [];
vertices_inline = zeros(size(vertex_nums));
for v = 1:length(vertex_nums)
    vertices_inline(v) = find(vertices_inline_data.data==vertex_nums(v));
end



winnerinds = find(winners.data~=0);

[~,winnerlineinds] = min(dmat(winnerinds,vertices_inline),[],2);

networklinevals = networks.data(winnerinds);

% winnerinds(winnerlineinds==1) = [];
% winnerlineinds(winnerlineinds==1) = [];
% if strcmp(subname,'SIC02')
%     winnerinds(winnerlineinds==38) = [];
%     winnerlineinds(winnerlineinds==38) = [];
% end

winnerlineinds_spots = winnerlineinds(networklinevals==1.5);
winnerlineinds_spots((winnerlineinds_spots > min(winnerlineinds(networklinevals==10))) & (winnerlineinds_spots < max(winnerlineinds(networklinevals==10)))) = [];
winnerlineinds_spots((winnerlineinds_spots > min(winnerlineinds(networklinevals==11))) & (winnerlineinds_spots < max(winnerlineinds(networklinevals==11)))) = [];
winnerlineinds_spots((winnerlineinds_spots < max(winnerlineinds(networklinevals==17)))) = [];
winnerlineinds_spots_only_sorted = unique(sort(winnerlineinds_spots));


%winnerlineinds_spots_only_sorted = sort(setdiff(winnerlineinds_spots,winnerlineinds(networklinevals>=10)),'ascend');
spots_diff_inds = find(diff([0 ; winnerlineinds_spots_only_sorted])>10);
winnerlineinds_spots_edges = sort([winnerlineinds_spots_only_sorted(spots_diff_inds); winnerlineinds_spots_only_sorted(spots_diff_inds(2:end)-1)]);

activationmaps = ft_read_cifti_mod(activationfile);
activationmaps.data = activationmaps.data(:,mapinds);
activationmaps.mapname = activationmaps.mapname(mapinds);
mapnum_to_color_mapping = mapnum_to_color_mapping_orig(mapinds);
fitparam = fitparam_orig(mapinds);
colors = colors_orig(mapinds);


Rsq_onegauss = zeros(size(activationmaps.data,2),1);
Rsq_twogauss = Rsq_onegauss;
Rsq_threegauss = Rsq_onegauss;
Fstats_1V2 = Rsq_onegauss;
pvals_1V2 = Rsq_onegauss;
Fstats_2V3 = Rsq_onegauss;
pvals_2V3 = Rsq_onegauss;



F2s = cell(1,size(activationmaps.data,2));
F3s = cell(1,size(activationmaps.data,2));

Lowessresults = cell(2,size(activationmaps.data,2));



for m = 1:size(activationmaps.data,2)
    
    if ~isnan(mapnum_to_color_mapping(m))
        activationvals = activationmaps.data(winnerinds,m);
        
        winnerlineinds_temp = winnerlineinds + rand(size(winnerlineinds)).*.00001;
        [~,sorti] = sort(winnerlineinds_temp,'ascend');
        activationvals_temp = [activationvals; zeros(10,1)];
        winnerlineinds_temp = [winnerlineinds_temp; [max(winnerlineinds_temp)+1 : max(winnerlineinds_temp)+10]'];
        
        dataout = lowess([winnerlineinds_temp activationvals_temp],.11);
        
        %[sorted, sorti] = sort(winnerlineinds);
        Lowessresults{1,m} = dataout(:,1);
        Lowessresults{2,m} = dataout(:,3);
        
        if make_eachplot
            
            figure;
            set(gcf,'Position',[813 30 1102 805])
            set(gcf,'Color',[1 1 1]);%set(gcf,'Color',[.9 .9 .9])
            set(gca,'Color',[1 1 1])
            plot(winnerlineinds.*2,activationvals,'k.','MarkerSize',15);
            hold on;
            
            
            plot(dataout(:,1).*2,dataout(:,3),'b-','LineWidth',6)
            
            set(gca,'FontSize',30)
            xlim([0 212])
            axis1 = gca;
            axis1.XTick = axis1.XTick+12;
            axis1.XTickLabel = {'200' '150' '100' '50' '0'};
            if export_plots
                export_fig([subname '_' activationmaps.mapname{m} '_modelfits.pdf'],gcf)
            end
            
            title(activationmaps.mapname{m})
            
            
            
        end
        
        
    end
end

all_colors = [all_colors colors];

%%

min_linepos = 0;
max_linepos = 200;
boxwidth = 2;


figure;
set(gcf,'Position',[555 30 1799 805])
set(gcf,'Color',[1 1 1]);%set(gcf,'Color',[.9 .9 .9])
set(gca,'Color',[1 1 1])
axis1 = gca;
for m = 1:size(Lowessresults,2)
    %plot(axis1,all_F2s{m}(1:106),212-(sort(1:106).*2),'Color',all_colors{m} ./ 255,'LineWidth',5)
    %plot(axis1,sort(1:106).*2,all_F2s{m}(1:106),'Color',all_colors{m} ./ 255,'LineWidth',4)
    %plot(axis1,all_F3s{m}(1:106),212-(sort(1:106).*2),'Color',all_colors{m} ./ 255,'LineWidth',5)
    plot(axis1,(Lowessresults{1,m}.*2),Lowessresults{2,m},'Color',all_colors{m} ./ 255,'LineWidth',4)
    hold on
end



for e = 1:length(winnerlineinds_spots_edges)
    h = vline(winnerlineinds_spots_edges(e) .* 2);
    h.Color = [0.5, 0, 0.4];
    h.LineWidth = 4;
end



rectangle('Position',[0 axis1.YLim(2)-boxwidth winnerlineinds_spots_edges(1).*2 boxwidth],'FaceColor',[0 .4 0])
rectangle('Position',[winnerlineinds_spots_edges(1).*2 axis1.YLim(2)-boxwidth (winnerlineinds_spots_edges(2) - winnerlineinds_spots_edges(1)).*2 boxwidth],'FaceColor',[0.5, 0, 0.4])
rectangle('Position',[winnerlineinds_spots_edges(2).*2 axis1.YLim(2)-boxwidth (winnerlineinds_spots_edges(3) - winnerlineinds_spots_edges(2)).*2 boxwidth],'FaceColor',[.2 1 1])
rectangle('Position',[winnerlineinds_spots_edges(3).*2 axis1.YLim(2)-boxwidth (winnerlineinds_spots_edges(4) - winnerlineinds_spots_edges(3)).*2 boxwidth],'FaceColor',[0.5, 0, 0.4])
rectangle('Position',[winnerlineinds_spots_edges(4).*2 axis1.YLim(2)-boxwidth (winnerlineinds_spots_edges(5) - winnerlineinds_spots_edges(4)).*2 boxwidth],'FaceColor',[1 .5 0])
rectangle('Position',[winnerlineinds_spots_edges(5).*2 axis1.YLim(2)-boxwidth max_linepos boxwidth],'FaceColor',[0.5, 0, 0.4])


xlim([min_linepos max_linepos])
set(gca,'XTick',[0    50   100   150   200])
set(gca,'XTickLabel',{'200','150','100','50','0'})
%xlim([0.01 axis1.XLim(2)+2])

% for m = 1:length(all_F2s)
%     if ~(strcmp(subname,'SIC02') && all(all_colors{m}==[51 255 255])) && ~(strcmp(subname,'SIC01') && all(all_colors{m}==[0 51 0]))
%         %plot(axis1,[axis1.XLim(2)-2 axis1.XLim(2)],212-[all_F2s{m}.b1 all_F2s{m}.b1].*2,'Color',all_colors{m} ./ 255,'LineWidth',4)
%         arrow([axis1.XLim(2) 212-(all_F2s{m}.b1).*2],[axis1.XLim(2)-2 212-(all_F2s{m}.b1).*2],'Color',all_colors{m} ./ 255,'LineWidth',2)
%     end
%     if ~( all(all_colors{m}==[0 255 0])) && ~(strcmp(subname,'SIC01') && all(all_colors{m}==[51 255 255])) && ~(strcmp(subname,'SIC01') && all(all_colors{m}==[248 216 0]))
%         %plot(axis1,[axis1.XLim(2)-2 axis1.XLim(2)],212-[all_F2s{m}.b2 all_F2s{m}.b2].*2,'Color',all_colors{m} ./ 255,'LineWidth',4)
%         arrow([axis1.XLim(2) 212-(all_F2s{m}.b2).*2],[axis1.XLim(2)-2 212-(all_F2s{m}.b2).*2],'Color',all_colors{m} ./ 255,'LineWidth',2)
%     end
%     hold on
% end
set(gca,'FontSize',30)

set(gca,'Color',[.9 .9 .9])

box off




