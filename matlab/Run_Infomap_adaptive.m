function Run_Infomap_adaptive(rmat, dmatname, xdistance, thresholdarray, makebinary, outdir, ignoreverts, numreps,numpools,newinfomapversion)
%Run_Infomap_adaptive(rmat, dmatname, xdistance, thresholdarray, makebinary, outdir, [ignoreverts], [numreps], [numpools],[newinfomapversion])
%
% Run infomap on a matrix with a given distance exclusion, at various
% density thresholds, and write the results from all thresholds into a
% single text file named "rawassn.txt". This can take a long time for large
% matrices. It will run up to eight infomaps simultaneously if the Parallel
% Computing Toolbox is installed.
%
% Inputs:
%
% rmat - a correlation matrix to be infomapped. Can be a numeric matrix or
%  a cifti file that will be loaded
% dmatname - a .mat file containing a node-to-node distance matrix
% xdistance - the distance exclusion to apply, in mm (i.e., nodes closer
%  than xdistance are not allowed to be connected)
% thresholdarray - a vector of thresholds to apply to the matrix. Infomap
%  will be run separately for each threshold.
% makebinary - whether the matrix is binarized after thresholding. 1 = make
%  it binary; 0 = leave it weighted.
% outdir - the folder results will be written to. Will be created if it
%  doesn't exist.
% ignoreverts - an OPTIONAL cifti input containing vertices that should be
%  ignored (all connections set to zero)
% numreps - an OPTIONAL scalar input specifying the number of repetitions
%  infomap should run. Omit or leave empty ([]) to use the default of 100.
% numpools - an OPTIONAL scalar input specifying the number of parallel
%  pools to use. Omit or leave empty ([]) to use the default of 8.
%
% Requires the Resources scripts to be in your path (e.g.,
% /home/data/scripts/Resources/ and subfolders)
%
%EMG 06/25/20


if ~exist('numreps') || isempty(numreps)
    numreps = 100;
end

if ~exist('numpools') || isempty(numpools)
    numpools = 6;
end

if ~exist('newinfomapversion')
    newinfomapversion = true;
end

if ~exist(outdir)
    mkdir(outdir);
end

if ischar(ignoreverts)
    ignoreverts = ft_read_cifti_mod(ignoreverts);
    ignoreverts = find(ignoreverts.data);
end


dlmwrite([outdir '/thresholds.txt'],thresholdarray,'delimiter',' ')


[~,~] = system(['rm ' outdir '/pajek*']);


string = ['loading correlations...'];
fprintf(string);
%fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
%prevstring = string;

if ischar(rmat)
    rmat = ft_read_cifti_mod(rmat);
    rmat = rmat.data;
end
rmat = single(rmat);

warning off

string = ['applying distance exclusion...'];
fprintf(string);
%fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
%prevstring = string;
% tic

if isnumeric(dmatname)
    dmat = dmatname;
elseif strcmp(dmatname(end-9:end),'.dconn.nii')
    dmat = ft_read_cifti_mod(dmatname); dmat = dmat.data;
else
    dmat = smartload(dmatname);
end
clear dmatname
dmat = uint8(dmat);

% apply a distance exclusion?
if ~isempty(xdistance)
    if isnumeric(xdistance) && (xdistance>=0)
        rmat(dmat < xdistance) = 0;
    else
        error('xdistance is not >=0 or is not numeric.\n');
    end
end

clear dmat
%toc

%prevent subcort-to-subcort connections
rmat(59413:end,59413:end) = 0;

rmat(ignoreverts,:) = 0;
rmat(:,ignoreverts) = 0;

num_nodes = size(rmat,1);



string = ['calculating r thresholds...'];
fprintf(string);

for t = 1:length(thresholdarray)
    connectionmats{t} = false(num_nodes,num_nodes);
    num_top_bythresh(t) = ceil(num_nodes .* thresholdarray(t));
end


for i = 1:num_nodes
    if any(rmat(:,i))
        [~,sortinds] = sort(rmat(:,i),'descend');
        for t = 1:length(thresholdarray)
            connectionmats{t}(sortinds(1:num_top_bythresh(t)),i) = true;
            connectionmats{t}(i,sortinds(1:num_top_bythresh(t))) = true;
        end
    end
end
clear sortinds
if makebinary
    clear rmat
end



string = ['saving pajek files...'];
fprintf(string);
for t = 1:length(thresholdarray)
    inds = find(triu(connectionmats{t},1));
    pajekfile = [outdir '/pajek_col' num2str(t) '.net'];
    if makebinary
        mat2pajek_byindex(connectionmats{t},inds,pajekfile);
    else
        mat2pajek_byindex(rmat,inds,pajekfile);
    end
end

clear rmat connectionmats inds

string = ['running infomap...'];
fprintf(string);
if newinfomapversion
    if ~any(size(gcp('nocreate')))
        parpool(numpools);
    end
else
    eval(['if ~matlabpool(''size''); matlabpool open ' num2str(numpools) '; end'])
end
parfor i=1:length(thresholdarray)
    pajekfile = [ outdir '/pajek_col' num2str(i) '.net' ];
    rawclrs = run_infomap_on_pajekfile(pajekfile,numreps);
    dlmwrite([outdir '/rawassn_col' num2str(i) '.txt'],rawclrs,'\t')
end
if newinfomapversion
    %delete(gcp('nocreate'))
else
    eval('matlabpool close')
end

for i = 1:length(thresholdarray)
    pajekfile = [ outdir '/pajek_col' num2str(i) '.net' ];
    delete(pajekfile);
    delete([pajekfile(1:end-4) '.clu']);
end



for i = 1:length(thresholdarray)
    rawclrs_all(:,i) = load([outdir '/rawassn_col' num2str(i) '.txt']);
end
% write the raw assignments as .txt
dlmwrite([outdir '/rawassn.txt'],rawclrs_all,'\t');
delete([outdir '/rawassn_col*.txt'])




end


