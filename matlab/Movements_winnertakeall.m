%%
subname = 'SIC01';
smoothname = '6.0';
%%

alldata_zstats = ft_read_cifti_mod([subname '_ZstatsMotorContrasts_MEonly_smooth' smoothname '.dscalar.nii']);




motorregions = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' subname '_precast_infomap_wacky2_subcortreg_ignoreverts/' subname '_rawassn_minsize10_regularized_motorrestricted.dtseries.nii']);
motor_spots_CON = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' subname '_precast_infomap_wacky2_subcortreg_ignoreverts/' subname '_rawassn_minsize10_regularized_CONandmotor_oneID_central_precentral.dtseries.nii']);
motor_spots_CON.data(59413:end) = 0;
motormask = true(size(motorregions.data));%motorregions.data>0 | motor_spots_CON.data>0; %motormask(motorregions.brainstructure(motorregions.brainstructure>0)>2) = false;

percentile_thresh = .9;


alldata_zstats_percentile = alldata_zstats;
alldata_zstats_percentile.data(:) = 0;
alldata_zstats_percentile.data(motormask,:) = rankorder(alldata_zstats.data(motormask,:),1,'ascend') ./ nnz(motormask);
ft_write_cifti_mod([subname '_ZstatsMotorContrastsPctile_inmask_smooth'],alldata_zstats_percentile)


[sortedvals, sorti] = sort(alldata_zstats_percentile.data,2,'descend');
sorti_bigenough = sorti .* (sortedvals>0);

sorti_bigenough = sorti_bigenough(:,any(sorti_bigenough>0,1));
out = alldata_copes;
out.data = sorti_bigenough;
recoloring = [1 1.5 2 2.5 3 4 5 6 6.6 7 7.5 8 9 11.5 12 13 15 16 17 18 11 10.4 10]; %1.33 11.25 %3.75 
out.data(out.data>0) = recoloring(out.data(out.data>0));

ft_write_cifti_mod([subname '_motor_allwinners'],out)
out.data = out.data .* (sortedvals>percentile_thresh);
ft_write_cifti_mod([subname '_motor_allwinners_threshp9'],out)


