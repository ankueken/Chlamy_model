function Model_splitted = split_rxns(Model)

% Input: cobra model

if ~isfield(Model, 'rxns')
    disp('Error: Field "rxns" is missing.')
elseif ~isfield(Model, 'S')
    disp('Error: Field "S" is missing.')
elseif ~isfield(Model, 'rev')
    disp('Error: Field "rev" is missing.')
elseif ~isfield(Model, 'lb')
    disp('Error: Field "lb" is missing.')    
elseif ~isfield(Model, 'ub')
    disp('Error: Field "ub" is missing.')
elseif ~isfield(Model, 'c')
    disp('Error: Field "c" is missing.')
end

% splitted Model
Model_splitted.S = zeros(size(Model.S,1),0);
Model_splitted.rxns = cell(0);
Model_splitted.lb = zeros(0);
Model_splitted.ub = zeros(0);
Model_splitted.rev = zeros(0);
Model_splitted.c = zeros(0);

if isfield(Model, 'rxnNames')
    Model_splitted.rxnNames = cell(0);
end
if isfield(Model, 'rxnECNumbers')
    Model_splitted.rxnECNumbers = cell(0);
end
if isfield(Model, 'rules')
    Model_splitted.rules = cell(0);
end
if isfield(Model, 'rxnGeneMat')
    Model_splitted.rxnGeneMat = zeros(0,size(Model.rxnGeneMat,2));
end
if isfield(Model, 'grRules')
    Model_splitted.grRules = cell(0);
end
if isfield(Model, 'subSystems')
    Model_splitted.subSystems = cell(0);
end
if isfield(Model, 'rxnReferences')
    Model_splitted.rxnReferences = cell(0);
end
if isfield(Model, 'rxnNotes')
    Model_splitted.rxnNotes = cell(0);
end
if isfield(Model, 'confidenceScore')
    Model_splitted.confidenceScore = zeros(0);
end
if isfield(Model, 'proteins')
    Model_splitted.proteins = cell(0);
end

for i=1:length(Model.rxns)
    Model_splitted.S(:,end+1) = Model.S(:,i);
    Model_splitted.rxns{end+1,1} = Model.rxns{i};
	if Model.lb(i)>0
		Model_splitted.lb(end+1,1) = Model.lb(i);
	else
		Model_splitted.lb(end+1,1) = 0;
	end
    Model_splitted.ub(end+1,1) = Model.ub(i);
    Model_splitted.rev(end+1,1) = 0;
    Model_splitted.c(end+1,1) = Model.c(i);
    
    if isfield(Model, 'rxnNames')
        Model_splitted.rxnNames{end+1,1} = Model.rxnNames{i};
    end
    if isfield(Model, 'rxnECNumbers')
        Model_splitted.rxnECNumbers{end+1,1} = Model.rxnECNumbers{i};
    end
    if isfield(Model, 'rules')
        Model_splitted.rules{end+1,1} = Model.rules{i};
    end
    if isfield(Model, 'rxnGeneMat')
        Model_splitted.rxnGeneMat(end+1,:) = Model.rxnGeneMat(i,:);
    end
    if isfield(Model, 'grRules')
        Model_splitted.grRules{end+1,1} = Model.grRules{i};
    end
    if isfield(Model, 'subSystems')
        Model_splitted.subSystems{end+1,1} = Model.subSystems{i};
    end
    if isfield(Model, 'rxnReferences')
        Model_splitted.rxnReferences{end+1,1} = Model.rxnReferences{i};
    end
    if isfield(Model, 'rxnNotes')
        Model_splitted.rxnNotes{end+1,1} = Model.rxnNotes{i};
    end
    if isfield(Model, 'confidenceScore')
        Model_splitted.confidenceScore{end+1,1} = Model.confidenceScore{i};
    end
    if isfield(Model, 'proteins')
        Model_splitted.proteins{end+1,1} = Model.proteins{i};
    end
    
    if Model.rev(i)==1
        Model_splitted.S(:,end+1) = Model.S(:,i)*-1;
        Model_splitted.rxns{end+1,1} = strcat(Model.rxns{i}, '_rev');
        Model_splitted.lb(end+1,1) = 0;
        Model_splitted.ub(end+1,1) = Model.lb(i)*-1;
        Model_splitted.rev(end+1,1) = 0;
        Model_splitted.c(end+1,1) = Model.c(i);

        if isfield(Model, 'rxnNames')
            Model_splitted.rxnNames{end+1,1} = Model.rxnNames{i};
        end
        if isfield(Model, 'rxnECNumbers')
            Model_splitted.rxnECNumbers{end+1,1} = Model.rxnECNumbers{i};
        end
        if isfield(Model, 'rules')
            Model_splitted.rules{end+1,1} = Model.rules{i};
        end
        if isfield(Model, 'rxnGeneMat')
            Model_splitted.rxnGeneMat(end+1,:) = Model.rxnGeneMat(i,:);
        end
        if isfield(Model, 'grRules')
            Model_splitted.grRules{end+1,1} = Model.grRules{i};
        end
        if isfield(Model, 'subSystems')
            Model_splitted.subSystems{end+1,1} = Model.subSystems{i};
        end
        if isfield(Model, 'rxnReferences')
            Model_splitted.rxnReferences{end+1,1} = Model.rxnReferences{i};
        end
        if isfield(Model, 'rxnNotes')
            Model_splitted.rxnNotes{end+1,1} = Model.rxnNotes{i};
        end
        if isfield(Model, 'confidenceScore')
            Model_splitted.confidenceScore{end+1,1} = Model.confidenceScore{i};
        end
        if isfield(Model, 'proteins')
            Model_splitted.proteins{end+1,1} = Model.proteins{i};
        end
        
    end
end
Model_splitted.b = zeros(size(Model_splitted.S,1),1);
Model_splitted.mets = Model.mets;
if isfield(Model, 'metNames')
    Model_splitted.metNames = Model.metNames;
end
if isfield(Model, 'metFormulas')
    Model_splitted.metFormulas = Model.metFormulas;
end
if isfield(Model, 'genes')
    Model_splitted.genes = Model.genes;
end

end
