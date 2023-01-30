subname = 'SIC01';

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

Face = true;
LE = true;
UE = false;

mapnum_to_color_mapping_orig = [1 2 2 3 3 4 5 6 7 7 8 9 9 12 12 13 15 15 17 17 11 10 10];
fitparam_orig = [50 20 20 50 50 50 50 50 50 50 50 50 50 50 50 20 20 20 20 20 50 50 50];
colors_orig = {[12 63 216],[25 127 229],[25 127 229],[38 191 242],[38 191 242],[0 0 204],[225 155 0],[237 186 0],[0 51 0],[0 51 0],[0 255 0],[0 187 0],[0 187 0],[0 119 0],[0 119 0],[214 124 0],[202 94 0],[202 94 0],[191 63 0],[191 63 0],[248 216 0],[51 255 255],[51 255 255]};


make_eachplot = true;
export_plots = true;

dmat = smartload(dmatname);

all_F2s = cell(0,1);
all_colors = [];




%% Face

if Face

vertices_inline_file = '/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/SIC02_precast_infomap_wacky2_subcortreg_ignoreverts/SIC02_precentral_v3_border_alt_continuousline_bothhems.dtseries.nii';

if strcmp(hem,'L')
    mapinds = [7 8 16 18 20 21 ];
else
    mapinds = [7 8 16 17 19 21];
end
max_position = 10000;
min_position = 30;


%%


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

% winnerinds(winnerlineinds==1) = [];
% winnerlineinds(winnerlineinds==1) = [];
% if strcmp(subname,'SIC02')
%     winnerinds(winnerlineinds==38) = [];
%     winnerlineinds(winnerlineinds==38) = [];
% end

winnerinds_orig = winnerinds;
winnerlineinds_orig = winnerlineinds;

winnerinds(winnerlineinds<min_position) = [];
winnerlineinds(winnerlineinds<min_position) = [];

activationmaps = ft_read_cifti_mod(activationfile);
activationmaps.data = activationmaps.data(:,mapinds);
activationmaps.mapname = activationmaps.mapname(mapinds);
mapnum_to_color_mapping = mapnum_to_color_mapping_orig(mapinds);
fitparam = fitparam_orig(mapinds);
colors = colors_orig(mapinds);

Rsq_onegauss = zeros(size(activationmaps.data,2),1);
Rsq_twogauss = Rsq_onegauss;
Fstats = Rsq_onegauss;
pvals = Rsq_onegauss;
freedman_lane_ps = Rsq_onegauss;




F2s = cell(1,size(activationmaps.data,2));

for m = 1:size(activationmaps.data,2)
    
    if ~isnan(mapnum_to_color_mapping(m))
        activationvals_orig = activationmaps.data(winnerinds_orig,m);
        activationvals = activationmaps.data(winnerinds,m);
        
        
        
        startpoints = [15 50 5 13 23 fitparam(m)];
        
        fitrange = [20 :20 :80];
        Rsqtest = zeros(length(fitrange),1);
        for f = 1:length(fitrange)
            test_startpoints = startpoints;
            test_startpoints(6) = fitrange(f);
            try
                
                [~, G] = fit(winnerlineinds,activationvals,'gauss1','StartPoint',test_startpoints(4:end),'Lower',[-Inf -Inf -Inf]);
                Rsqtest(f) = G.adjrsquare;
            catch
                try
                    [~, G] = fit(winnerlineinds,activationvals,'gauss1','StartPoint',test_startpoints(4:end),'Lower',[0 0 0]);
                Rsqtest(f) = G.adjrsquare;
                catch
                end
            end
        end
        [Rsqmax,fmaxi] = max(Rsqtest);
        if Rsqmax>0
            startpoints(6) = fitrange(fmaxi);
            
            if strcmp(subname,'SIC02') && strcmp(hem,'L') && (mapinds(m)==8 || mapinds(m)==21)
                [F, G] = fit(winnerlineinds,activationvals,'gauss1','Lower',[0 70 -Inf],'Upper',[Inf 100 Inf]);
                
            else
                
                try
                    [F, G] = fit(winnerlineinds,activationvals,'gauss1','StartPoint',startpoints(4:end),'Lower',[-Inf -Inf -Inf]);
                catch
                    [F, G] = fit(winnerlineinds,activationvals,'gauss1','StartPoint',startpoints(4:end),'Lower',[0 0 0]);
                end
            end
            Rsq_onegauss(m) = G.adjrsquare;
            worked = 1;
        else
            worked = 0;
        end
        
        
        if strcmp(subname,'SIC02') && strcmp(hem,'L') && (mapinds(m)==18 || mapinds(m)==20)
            [F2, G2] = fit(winnerlineinds,activationvals,'gauss2','Lower',[0 80 -Inf 0 0 -Inf],'Upper',[Inf 100 Inf Inf Inf Inf]);
        elseif strcmp(subname,'SIC01') && strcmp(hem,'L') && (mapinds(m)==20)
             [F2, G2] = fit(winnerlineinds,activationvals,'gauss2','Lower',[0 80 3 0 0 3],'Upper',[Inf 100 Inf Inf 80 Inf]);
        elseif strcmp(subname,'SIC01') && strcmp(hem,'L') && (mapinds(m)==21)
             [F2, G2] = fit(winnerlineinds,activationvals,'gauss2','Lower',[0 0 3 0 0 3],'Upper',[Inf 90 Inf Inf 90 Inf]);
        elseif strcmp(subname,'SIC02') && strcmp(hem,'L') && (mapinds(m)==7)
             [F2, G2] = fit(winnerlineinds,activationvals,'gauss2','Lower',[0 -Inf -Inf 0 -Inf -Inf],'Upper',[Inf Inf Inf Inf Inf 10]);
        else
            [F2, G2] = fit(winnerlineinds,activationvals,'gauss2','Lower',[0 -Inf -Inf 0 -Inf -Inf]);
        end
        
        
        if make_eachplot
            
            figure;
            set(gcf,'Position',[813 30 1102 805])
            set(gcf,'Color',[1 1 1]);%set(gcf,'Color',[.9 .9 .9])
            set(gca,'Color',[1 1 1])
            plot(winnerlineinds_orig.*2,activationvals_orig,'k.','MarkerSize',15);
            hold on;
            plot(sort(winnerlineinds).*2,F2(sort(winnerlineinds)),'b-','LineWidth',6)
            
            
            
            if worked
                Fstats(m) = ((G.sse-G2.sse) ./ (G.dfe - G2.dfe)) ./ (G2.sse ./ G2.dfe);
                pvals(m) = 1 - fcdf(Fstats(m),(G.dfe-G2.dfe),G2.dfe);
                %permuted_rsqs = freedman_lane_permutation_doublegaussian(winnerlineinds,activationvals, F2, 1000);
                %freedman_lane_ps(m) = nnz(permuted_rsqs > G2.rsquare) ./ numel(permuted_rsqs);
                plot(sort(winnerlineinds).*2,F(sort(winnerlineinds)),'r-','LineWidth',4)
                
            end
            set(gca,'FontSize',30)
            xlim([0 212])
            axis1 = gca;
            axis1.XTick = axis1.XTick+12;
            axis1.XTickLabel = {'200' '150' '100' '50' '0'};
            if export_plots
                axis1.XTickLabel = [];
                axis1.YTickLabel = [];
                export_fig([subname '_' activationmaps.mapname{m} '_modelfits_alldata_noaxlab.pdf'],gcf)
            end
            
            title(activationmaps.mapname{m})
            
            Rsq_twogauss(m) = G2.adjrsquare;
            
            
            disp([activationmaps.mapname{m} ': Fstat of difference=' num2str(Fstats(m)) '; p value of difference=' num2str(pvals(m)) '; Freedman-Lane p value=' num2str(freedman_lane_ps(m)) '; One gaussian Rsq=' num2str(Rsq_onegauss(m)) ', Two gaussian Rsq=' num2str(Rsq_twogauss(m))])
            
        end
        
        F2s{m} = F2;
        
    end
end

all_F2s = [all_F2s F2s];
all_colors = [all_colors colors];

end


%% LE

if LE

vertices_inline_file = '/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/SIC02_precast_infomap_wacky2_subcortreg_ignoreverts/SIC02_precentral_v3_border_alt_continuousline_bothhems.dtseries.nii';
if strcmp(hem,'L')
    mapinds = [10 11 13 15 ];%LHem   [22 6 4 2 1];%RHem        [1 3 5 6 23];%[1 3 5 6 7 8 10 11 13 15 16 18 20 21 23];
else
    mapinds = [9 11 12 14 ];
end
max_position = 70;


%%
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


winnerinds(winnerlineinds==1) = [];
winnerlineinds(winnerlineinds==1) = [];
if strcmp(subname,'SIC02')
    winnerinds(winnerlineinds==38) = [];
    winnerlineinds(winnerlineinds==38) = [];
end
winnerinds_orig = winnerinds;
winnerlineinds_orig = winnerlineinds;

winnerinds(winnerlineinds>max_position) = [];
winnerlineinds(winnerlineinds>max_position) = [];
% winnerinds(winnerlineinds<min_position) = [];
% winnerlineinds(winnerlineinds<min_position) = [];

activationmaps = ft_read_cifti_mod(activationfile);
activationmaps.data = activationmaps.data(:,mapinds);
activationmaps.mapname = activationmaps.mapname(mapinds);
mapnum_to_color_mapping = mapnum_to_color_mapping_orig(mapinds);
fitparam = fitparam_orig(mapinds);
colors = colors_orig(mapinds);

Rsq_onegauss = zeros(size(activationmaps.data,2),1);
Rsq_twogauss = Rsq_onegauss;
Fstats = Rsq_onegauss;
pvals = Rsq_onegauss;
freedman_lane_ps = Rsq_onegauss;




F2s = cell(1,size(activationmaps.data,2));

for m = 1:size(activationmaps.data,2)
    
    if ~isnan(mapnum_to_color_mapping(m))
        activationvals_orig = activationmaps.data(winnerinds_orig,m);
        activationvals = activationmaps.data(winnerinds,m);
        
        
        
        startpoints = [15 50 5 13 23 fitparam(m)];
        
        fitrange = [20 :20 :80];
        Rsqtest = zeros(length(fitrange),1);
        for f = 1:length(fitrange)
            test_startpoints = startpoints;
            test_startpoints(6) = fitrange(f);
            try
                
                [~, G] = fit(winnerlineinds,activationvals,'gauss1','StartPoint',test_startpoints(4:end),'Lower',[-Inf -Inf -Inf]);
                Rsqtest(f) = G.adjrsquare;
            catch
            end
        end
        [Rsqmax,fmaxi] = max(Rsqtest);
        if Rsqmax>0
            startpoints(6) = fitrange(fmaxi);
            
            [F, G] = fit(winnerlineinds,activationvals,'gauss1','StartPoint',startpoints(4:end),'Lower',[-Inf -Inf -Inf]);
            Rsq_onegauss(m) = G.adjrsquare;
            worked = 1;
        else
            worked = 0;
        end
        
        if strcmp(subname,'SIC02') && strcmp(hem,'R') && mapinds(m)==1
            [F2, G2] = fit(winnerlineinds,activationvals,'gauss2','Lower',[0 40 -Inf 0 -Inf -Inf]);
        elseif strcmp(subname,'SIC01') && strcmp(hem,'L') && mapinds(m)==11
            [F2, G2] = fit(winnerlineinds,activationvals,'gauss2','Lower',[0 100 -Inf 0 -Inf -Inf]);
            
        else
            [F2, G2] = fit(winnerlineinds,activationvals,'gauss2','Lower',[0 -Inf -Inf 0 -Inf -Inf]);
        end
        
        
        
        if make_eachplot
            
            figure;
            set(gcf,'Position',[813 30 1102 805])
            set(gcf,'Color',[1 1 1]);%set(gcf,'Color',[.9 .9 .9])
            set(gca,'Color',[1 1 1])
            plot(winnerlineinds_orig.*2,activationvals_orig,'k.','MarkerSize',15);
            hold on;
            plot(sort(winnerlineinds).*2,F2(sort(winnerlineinds)),'b-','LineWidth',6)
            
            
            %title(activationmaps.mapname{m})
            
            if worked
                Fstats(m) = ((G.sse-G2.sse) ./ (G.dfe - G2.dfe)) ./ (G2.sse ./ G2.dfe);
                pvals(m) = 1 - fcdf(Fstats(m),(G.dfe-G2.dfe),G2.dfe);
                %permuted_rsqs = freedman_lane_permutation_doublegaussian(winnerlineinds,activationvals, F2, 1000);
                %freedman_lane_ps(m) = nnz(permuted_rsqs > G2.rsquare) ./ numel(permuted_rsqs);
                plot(sort(winnerlineinds).*2,F(sort(winnerlineinds)),'r-','LineWidth',4)
                
            end
            set(gca,'FontSize',30)
            xlim([0 212])
            axis1 = gca;
            axis1.XTick = axis1.XTick+12;
            axis1.XTickLabel = {'200' '150' '100' '50' '0'};
            if export_plots
                axis1.XTickLabel = [];
                axis1.YTickLabel = [];
                export_fig([subname '_' activationmaps.mapname{m} '_modelfits_alldata_noaxlab.pdf'],gcf)
            end
            
            title(activationmaps.mapname{m})
            
            Rsq_twogauss(m) = G2.adjrsquare;
            
            disp([activationmaps.mapname{m} ': Fstat of difference=' num2str(Fstats(m)) '; p value of difference=' num2str(pvals(m)) '; Freedman-Lane p value=' num2str(freedman_lane_ps(m)) '; One gaussian Rsq=' num2str(Rsq_onegauss(m)) ', Two gaussian Rsq=' num2str(Rsq_twogauss(m))])
            
        end
        
        F2s{m} = F2;
        
    end
end

all_F2s = [all_F2s F2s];
all_colors = [all_colors colors];

end

%% UE

if UE

vertices_inline_file = '/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/SIC02_precast_infomap_wacky2_subcortreg_ignoreverts/SIC02_precentral_v3_border_continuousline_bothhems.dtseries.nii';

if strcmp(hem,'L')

    mapinds = [23 6 5 3 1];%LHem   
    
else
    
    mapinds = [22 6 4 2 1];%RHem        
end


max_position = 70;

%%

if strcmp(subname,'SIC02')
    
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
    
    
    
    winnerinds(winnerlineinds==1) = [];
    winnerlineinds(winnerlineinds==1) = [];
    if strcmp(subname,'SIC02')
        winnerinds(winnerlineinds==38) = [];
        winnerlineinds(winnerlineinds==38) = [];
    end
    
    winnerinds_orig = winnerinds;
    winnerlineinds_orig = winnerlineinds;
    
    winnerinds(winnerlineinds>max_position) = [];
    winnerlineinds(winnerlineinds>max_position) = [];
    % winnerinds(winnerlineinds<min_position) = [];
    % winnerlineinds(winnerlineinds<min_position) = [];
    
    winnerlineinds = winnerlineinds+6;
    winnerlineinds_orig = winnerlineinds_orig+6;
    
    
    
    activationmaps = ft_read_cifti_mod(activationfile);
    activationmaps.data = activationmaps.data(:,mapinds);
    activationmaps.mapname = activationmaps.mapname(mapinds);
    mapnum_to_color_mapping = mapnum_to_color_mapping_orig(mapinds);
    fitparam = fitparam_orig(mapinds);
    colors = colors_orig(mapinds);
    
    Rsq_onegauss = zeros(size(activationmaps.data,2),1);
    Rsq_twogauss = Rsq_onegauss;
    Fstats = Rsq_onegauss;
    pvals = Rsq_onegauss;
    freedman_lane_ps = Rsq_onegauss;
    
    
    
    
    F2s = cell(1,size(activationmaps.data,2));
    
    for m = 1:size(activationmaps.data,2)
        
        if ~isnan(mapnum_to_color_mapping(m))
            activationvals_orig = activationmaps.data(winnerinds_orig,m);
            activationvals = activationmaps.data(winnerinds,m);
            
            
            
            startpoints = [15 50 5 13 23 fitparam(m)];
            
            fitrange = [20 :20 :80];
            Rsqtest = zeros(length(fitrange),1);
            for f = 1:length(fitrange)
                test_startpoints = startpoints;
                test_startpoints(6) = fitrange(f);
                try
                    
                    [~, G] = fit(winnerlineinds,activationvals,'gauss1','StartPoint',test_startpoints(4:end),'Lower',[-Inf -Inf -Inf]);
                    Rsqtest(f) = G.adjrsquare;
                catch
                end
            end
            [Rsqmax,fmaxi] = max(Rsqtest);
            if Rsqmax>0
                startpoints(6) = fitrange(fmaxi);
                
                [F, G] = fit(winnerlineinds,activationvals,'gauss1','StartPoint',startpoints(4:end),'Lower',[-Inf -Inf -Inf]);
                Rsq_onegauss(m) = G.adjrsquare;
                worked = 1;
            else
                worked = 0;
            end
            
            if strcmp(subname,'SIC02') && strcmp(hem,'R') && mapinds(m)==1
                [F2, G2] = fit(winnerlineinds,activationvals,'gauss2','Lower',[0 40 -Inf 0 -Inf -Inf]);
            elseif strcmp(subname,'SIC02') && strcmp(hem,'L') && mapinds(m)==1
                [F2, G2] = fit(winnerlineinds,activationvals,'gauss2','Lower',[0 30 0 5 0 2]);
                
            else
                [F2, G2] = fit(winnerlineinds,activationvals,'gauss2');
            end
            
            
            
            if make_eachplot
                
                figure;
                set(gcf,'Position',[813 30 1102 805])
                set(gcf,'Color',[1 1 1]);%set(gcf,'Color',[.9 .9 .9])
                set(gca,'Color',[1 1 1])
                plot(winnerlineinds_orig.*2,activationvals_orig,'k.','MarkerSize',15);
                hold on;
                plot(sort(winnerlineinds).*2,F2(sort(winnerlineinds)),'b-','LineWidth',6)
                
                
                %
                
                if worked
                    Fstats(m) = ((G.sse-G2.sse) ./ (G.dfe - G2.dfe)) ./ (G2.sse ./ G2.dfe);
                    pvals(m) = 1 - fcdf(Fstats(m),(G.dfe-G2.dfe),G2.dfe);
                    %permuted_rsqs = freedman_lane_permutation_doublegaussian(winnerlineinds,activationvals, F2, 1000);
                    %freedman_lane_ps(m) = nnz(permuted_rsqs > G2.rsquare) ./ numel(permuted_rsqs);
                    plot(sort(winnerlineinds).*2,F(sort(winnerlineinds)),'r-','LineWidth',4)
                    
                end
                set(gca,'FontSize',30)
                xlim([0 212])
                axis1 = gca;
                axis1.XTick = axis1.XTick+12;
                axis1.XTickLabel = {'200' '150' '100' '50' '0'};
                if export_plots
                    axis1.XTickLabel = [];
                    axis1.YTickLabel = [];
                    export_fig([subname '_' activationmaps.mapname{m} '_modelfits_alldata_noaxlab.pdf'],gcf)
                end
                
                title(activationmaps.mapname{m})
                
                Rsq_twogauss(m) = G2.adjrsquare;
                
                disp([activationmaps.mapname{m} ': Fstat of difference=' num2str(Fstats(m)) '; p value of difference=' num2str(pvals(m)) '; Freedman-Lane p value=' num2str(freedman_lane_ps(m)) '; One gaussian Rsq=' num2str(Rsq_onegauss(m)) ', Two gaussian Rsq=' num2str(Rsq_twogauss(m))])
                
            end
            
            F2s{m} = F2;
            
        end
    end
end




if strcmp(subname,'SIC01')
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
    
    %dmat = smartload(dmatname);
    
    winnerinds = find(winners.data~=0);
    
    [~,winnerlineinds] = min(dmat(winnerinds,vertices_inline),[],2);
    
    winnerinds(winnerlineinds==1) = [];
    winnerlineinds(winnerlineinds==1) = [];
    
    winnerlineinds = winnerlineinds+6;
    
    activationmaps = ft_read_cifti_mod(activationfile);
    activationmaps.data = activationmaps.data(:,mapinds);
    activationmaps.mapname = activationmaps.mapname(mapinds);
    mapnum_to_color_mapping = mapnum_to_color_mapping_orig(mapinds);
    fitparam = fitparam_orig(mapinds);
    colors = colors_orig(mapinds);
    
    Rsq_onegauss = zeros(size(activationmaps.data,2),1);
    Rsq_twogauss = Rsq_onegauss;
    Fstats = Rsq_onegauss;
    pvals = Rsq_onegauss;
    freedman_lane_ps = Rsq_onegauss;
    
    
    F2s = cell(1,size(activationmaps.data,2));
    
    for m = 1:size(activationmaps.data,2)
        
        if ~isnan(mapnum_to_color_mapping(m))
            
            activationvals = activationmaps.data(winnerinds,m);
            
            
            startpoints = [15 50 5 13 23 fitparam(m)];
            
            fitrange = [20 :20 :80];
            Rsqtest = zeros(length(fitrange),1);
            for f = 1:length(fitrange)
                test_startpoints = startpoints;
                test_startpoints(6) = fitrange(f);
                try
                    
                    [~, G] = fit(winnerlineinds,activationvals,'gauss1','StartPoint',test_startpoints(4:end),'Lower',[-Inf -Inf -Inf]);
                    Rsqtest(f) = G.adjrsquare;
                catch
                end
            end
            [Rsqmax,fmaxi] = max(Rsqtest);
            if Rsqmax>0
                startpoints(6) = fitrange(fmaxi);
                
                [F, G] = fit(winnerlineinds,activationvals,'gauss1','StartPoint',startpoints(4:end),'Lower',[-Inf -Inf -Inf]);
                Rsq_onegauss(m) = G.adjrsquare;
                worked = 1;
            else
                worked = 0;
            end
            
            fitrange = [20 :20 :80];
            Rsqtest = zeros(length(fitrange));
            for f = 1:length(fitrange)
                for f2 = 1:length(fitrange)
                    test_startpoints = startpoints;
                    test_startpoints(3) = fitrange(f);
                    test_startpoints(6) = fitrange(f2);
                    [~, G2] = fit(winnerlineinds,activationvals,'gauss2','StartPoint',startpoints,'Lower',[0 -Inf -Inf 0 -Inf -Inf]);
                    Rsqtest(f,f2) = G2.adjrsquare;
                end
            end
            [~,fmaxi] = max(max(Rsqtest,[],2),[],1);
            [~,f2maxi] = max(max(Rsqtest,[],1),[],2);
            startpoints(3) = fitrange(fmaxi);
            startpoints(6) = fitrange(f2maxi);
            
            if strcmp(hem,'L') && mapinds(m)==3
                [F2, G2] = fit(winnerlineinds,activationvals,'gauss2','StartPoint',startpoints,'Lower',[3 55 .5 0 -Inf -Inf],'Upper',[Inf, 80, Inf, Inf, Inf, Inf]);
            elseif strcmp(hem,'R') && mapinds(m)==1
                [F2, G2] = fit(winnerlineinds,activationvals,'gauss2','StartPoint',startpoints,'Lower',[0 50 0 0 -Inf -Inf],'Upper',[Inf, 70, Inf, Inf, Inf, Inf]);
            else
                [F2, G2] = fit(winnerlineinds,activationvals,'gauss2','StartPoint',startpoints,'Lower',[0 -Inf -Inf 0 -Inf -Inf]);
            end
            
            
            
            if make_eachplot
                
                figure;
                set(gcf,'Position',[813 30 1102 805])
                set(gcf,'Color',[1 1 1]);%set(gcf,'Color',[.9 .9 .9])
                set(gca,'Color',[1 1 1])
                plot(winnerlineinds.*2,activationvals,'k.','MarkerSize',15);
                hold on;
                plot(sort(winnerlineinds).*2,F2(sort(winnerlineinds)),'b-','LineWidth',6)
                
                
                
                
                if worked
                    Fstats(m) = ((G.sse-G2.sse) ./ (G.dfe - G2.dfe)) ./ (G2.sse ./ G2.dfe);
                    pvals(m) = 1 - fcdf(Fstats(m),(G.dfe-G2.dfe),G2.dfe);
                    %permuted_rsqs = freedman_lane_permutation_doublegaussian(winnerlineinds,activationvals, F2, 1000);
                    %freedman_lane_ps(m) = nnz(permuted_rsqs > G2.rsquare) ./ numel(permuted_rsqs);
                    plot(sort(winnerlineinds).*2,F(sort(winnerlineinds)),'r-','LineWidth',4)
                    
                end
                set(gca,'FontSize',30)
                xlim([0 212])
                axis1 = gca;
                axis1.XTick = axis1.XTick+12;
                axis1.XTickLabel = {'200' '150' '100' '50' '0'};
                if export_plots
                    axis1.XTickLabel = [];
                    axis1.YTickLabel = [];
                    export_fig([subname '_' activationmaps.mapname{m} '_modelfits_alldata_noaxlab.pdf'],gcf)
                end
                
                title(activationmaps.mapname{m})
                
                Rsq_twogauss(m) = G2.adjrsquare;
                
                disp([activationmaps.mapname{m} ': Fstat of difference=' num2str(Fstats(m)) '; p value of difference=' num2str(pvals(m)) '; Freedman-Lane p value=' num2str(freedman_lane_ps(m)) '; One gaussian Rsq=' num2str(Rsq_onegauss(m)) ', Two gaussian Rsq=' num2str(Rsq_twogauss(m))])
                
            end
            
            F2s{m} = F2;
            
        end
    end
end


all_F2s = [all_F2s F2s];
all_colors = [all_colors colors];

end


%%

min_linepos = 212;
max_linepos = 0;

if Face
    min_linepos = min([min_linepos 0]);
    max_linepos = max([max_linepos 100]);
end
if LE
    min_linepos = min([min_linepos 120]);
    max_linepos = max([max_linepos 212]);
end
if UE
    min_linepos = min([min_linepos 80]);
    max_linepos = max([max_linepos 180]);
end


figure;
set(gcf,'Position',[813 30 566 805])
set(gcf,'Color',[1 1 1]);%set(gcf,'Color',[.9 .9 .9])
set(gca,'Color',[1 1 1])
axis1 = gca;
for m = 1:length(all_F2s)
    plot(axis1,all_F2s{m}(1:106),212-(sort(1:106).*2),'Color',all_colors{m} ./ 255,'LineWidth',5)
    %plot(axis1,sort(1:106).*2,all_F2s{m}(1:106),'Color',all_colors{m} ./ 255,'LineWidth',4)
    hold on
end
ylim([min_linepos max_linepos])
xlim([0.01 axis1.XLim(2)+2])
for m = 1:length(all_F2s)
    if ~(strcmp(subname,'SIC02') && all(all_colors{m}==[51 255 255])) && ~(strcmp(subname,'SIC01') && all(all_colors{m}==[0 51 0]))
        %plot(axis1,[axis1.XLim(2)-2 axis1.XLim(2)],212-[all_F2s{m}.b1 all_F2s{m}.b1].*2,'Color',all_colors{m} ./ 255,'LineWidth',4)
        arrow([axis1.XLim(2) 212-(all_F2s{m}.b1).*2],[axis1.XLim(2)-2 212-(all_F2s{m}.b1).*2],'Color',all_colors{m} ./ 255,'LineWidth',2)
    end
    if ~( all(all_colors{m}==[0 255 0])) && ~(strcmp(subname,'SIC01') && all(all_colors{m}==[51 255 255])) && ~(strcmp(subname,'SIC01') && all(all_colors{m}==[248 216 0]))
        %plot(axis1,[axis1.XLim(2)-2 axis1.XLim(2)],212-[all_F2s{m}.b2 all_F2s{m}.b2].*2,'Color',all_colors{m} ./ 255,'LineWidth',4)
        arrow([axis1.XLim(2) 212-(all_F2s{m}.b2).*2],[axis1.XLim(2)-2 212-(all_F2s{m}.b2).*2],'Color',all_colors{m} ./ 255,'LineWidth',2)
    end
    hold on
end
set(gca,'FontSize',30)
%ylim([min_linepos max_linepos])
%xlim([0.01 axis1.XLim(2)])
%xlim([0 212])

set(gca,'Color',[.9 .9 .9])
%set(gcf,'Color',[.9 .9 .9])

box off

%disp([Rsq_onegauss Rsq_twogauss])

%ft_write_cifti_mod([subname '_centersurround_model_curvefit.dtseries.nii'],Rsq_twogauss_map);
%ft_write_cifti_mod([subname '_homunculus_model_curvefit.dtseries.nii'],Rsq_onegauss_map);




