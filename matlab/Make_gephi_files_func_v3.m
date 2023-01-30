function Make_gephi_files_func_v3(graph,assignments,assignmentnames,outname,otherfactors,otherfactornamnes,colormap)
%Make_gephi_files_func(graph,assignments,assignmentnames,outname,[otherfactors],[otherfactornames],[colormap])
warning off


templatefile = '/data/nil-bluearc/GMT/Evan/MSC/Analysis_V1/spring_embed/gephi_template.xls';

if ~exist('colormap','var')
    colormap = [1 0 0;0 0 .8;1 1 0;1 .8 .6;0 1 0;1 .6 1;0 .6 .6;0 0 0;.35 0 .65;.2 1 1;1 .5 0;.65 .25 1;0 .25 .6;.6 1 .6;.2 .3 1;1 1 1;0 .4 0; repmat([.25 .25 .25],50,1)];
end


    

if ~iscell(assignments)
    assignments = {assignments};
end
colors = cell(length(assignments),1);

for a = 1:length(assignments)

for i = 1:length(assignments{a})
    if assignments{a}(i)<1
        colors{a}(i,:) = [.67 .67 .67];
    elseif (assignments{a}(i) - floor(assignments{a}(i)))==0
        colors{a}(i,:) = colormap(assignments{a}(i),:);
    else
        
        base = floor(assignments{a}(i));
        factor = assignments{a}(i)-base;
        
        colors{a}(i,:) = (colormap(base,:) .* (1-factor)) + (colormap(base+1,:) .* factor);
        
    end
end

rgb{a} = cell(0,1);
for i = 1:length(colors{a})
    rgb{a}{i,1} = ['"' num2str(colors{a}(i,1)*255) ',' num2str(colors{a}(i,2)*255) ',' num2str(colors{a}(i,3)*255) '"'];
end
end
temp = dataset('XLSFile',templatefile);
IDs = (1:1:length(colors{a}))';
if exist('otherfactors') && ~isempty(otherfactors)
    sumfactors = sum(otherfactors,2);
    [~,sortedi] = sort(sumfactors,'ascend');
    for col = 1:size(otherfactors,2)
        factorcol = otherfactors(:,col);
        eval(['temp.' otherfactornamnes{col} ' = factorcol(sortedi);']);
    end
      
    for a = 1:length(assignments)
        rgb{a} = rgb{a}(sortedi);
    end
    assignments{1} = assignments{1}(sortedi);
    IDs = IDs(sortedi);
    
    
end

temp.ID = IDs;
for a = 1:length(assignments)
    if a==1
        evalc(['temp.Color = rgb{a};']);
    end
    evalc(['temp.' assignmentnames{a} '_color = rgb{a};']);
    evalc(['temp.' assignmentnames{a} '_value = assignments{a};']);
end


export(temp,'File',[outname '_gephi_nodes.csv'],'delimiter',',')




edgestemplatefile = '/data/nil-bluearc/GMT/Evan/MSC/Analysis_V1/spring_embed/gephi_edges_template.xls';
temp = dataset('XLSFile',edgestemplatefile);


graph = graph .* tril(true(size(graph)),-1);
[sources,targets] = find(graph);
weights = graph(find(graph));

temp.Source = sources;
temp.Target = targets;
temp.Type = repmat({'Undirected'},length(sources),1);
temp.Weight = weights;
temp.ID = [0:(length(sources)-1)]';
export(temp,'File',[outname '_gephi_edges.csv'],'delimiter',',')
         
            


