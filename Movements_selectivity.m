subname = 'SIC02';
smoothname = '6.0';


%%

alldata_zstats = ft_read_cifti_mod([subname '_ZstatsMotorContrasts_MEonly_smooth' smoothname '.dscalar.nii']);
motions_tocompare = [6 8 11 19 20 21 22 23];
alldata_zstats.data = alldata_zstats.data(:,motions_tocompare);

motor_spots = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Evan/CIMT/Subnetworks/' subname '_precast_infomap_wacky2_subcortreg_ignoreverts/' subname '_rawassn_minsize10_regularized_motoronly_oneID_CS.dtseries.nii']);
motor_spots.data(59413:end) = 0;
motormask = true(size(motor_spots.data));


alldata_zstats_percentile = alldata_zstats;
alldata_zstats_percentile.data(:) = 0;
alldata_zstats_percentile.data(motormask,:) = rankorder(alldata_zstats.data(motormask,:),1,'ascend') ./ nnz(motormask);
%ft_write_cifti_mod([subname '_ZstatsMotorContrastsPctile_inmask_smooth'],alldata_zstats_percentile)


[sortedvals, sorti] = sort(alldata_zstats.data,2,'descend');



out = alldata_zstats;
out.data = (sortedvals(:,1) - sortedvals(:,2)) .* single(logical(motor_spots.data));% .* single(motor_spots.data~=9);
out.dimord = 'pos_time';

ft_write_cifti_mod([subname '_motor_winner_delta_CS'],out)
