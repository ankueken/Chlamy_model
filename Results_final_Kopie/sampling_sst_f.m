function samplesMinFlux = sampling_sst_f(Chlamy_splitted, ExpDataProteinDistr,obj,cond)

% obj:  0 - no objective function
%      -1 - flux minimization

if strcmp(cond,'lowCO2')
    idx_c = 14;
elseif strcmp(cond,'highCO2')
    idx_c = 7;
end

Enzyme_ids_hs=[];
Enzyme_ids_y=[];
for i=2
    idx = find(strcmp(Chlamy_splitted.mets, strcat(ExpDataProteinDistr{1,i},' [s]')));
    Enzyme_ids_hs(i-1,1) = idx;
    idx = find(strcmp(Chlamy_splitted.mets, strcat(ExpDataProteinDistr{1,i},' [y]')));
    Enzyme_ids_y(i-1,1) = idx;
end

samplesFull=zeros(size(Chlamy_splitted.rxns,1),0); samplesMinFlux=zeros(size(Chlamy_splitted.rxns,1),0);
run=0;
while size(samplesMinFlux,2) < 10000
run=run+1;

% sample ratio from measured range
a = str2double(ExpDataProteinDistr{idx_c,2:end})-str2double(ExpDataProteinDistr{idx_c,2:end}).*(cell2mat(ExpDataProteinDistr{idx_c+1,2:end})/100);
a(a<0) = 0;
b = str2double(ExpDataProteinDistr{idx_c,2:end})+str2double(ExpDataProteinDistr{idx_c,2:end}).*(cell2mat(ExpDataProteinDistr{idx_c+1,2:end})/100);
b(b>1) = 1;
R = a + (b-a).*rand(1,length(b));

% find reaction around enzyme
id_y=[]; id_s=[]; ratio=[];
for i=1:size(Enzyme_ids_hs,1)
    enzyme_consuming_rxns_y = find(Chlamy_splitted.S(Enzyme_ids_y(i),:)~=0); % pyrenoid
    enzyme_consuming_rxns_s = find(Chlamy_splitted.S(Enzyme_ids_hs(i),:)~=0); % stroma
    
    id_y(end+1:end+length(enzyme_consuming_rxns_y)) = enzyme_consuming_rxns_y;
    id_s(end+1:end+length(enzyme_consuming_rxns_s)) = enzyme_consuming_rxns_s;
    ratio(end+1:end+length(enzyme_consuming_rxns_s)) = R(i); % ratio = y/s
end

A_ratio = zeros(length(ratio),size(Chlamy_splitted.S,2));
for i=1:length(ratio)
	A_ratio(i,id_y(i)) = 1-ratio(i);
    A_ratio(i,id_s(i)) = -ratio(i);
end

% 1 setup the model structure
Model.S = [Chlamy_splitted.S; A_ratio];
Model.rxns = Chlamy_splitted.rxns;
Model.mets = Chlamy_splitted.mets;

for i=1:length(ratio)
    Model.mets(end+1) = strcat('flux-' ,Model.rxns(id_y(i)),'/',Model.rxns(id_s(i)));
end

Model.rev = zeros(size(Model.rxns));

% 2 set objective
Model.c = ones(size(Chlamy_splitted.rxns))*obj;

% 3 set constraint bounds
Model.b = zeros(size(Model.mets));

% 4 set flux bounds
Model.lb=zeros(size(Model.rxns));
Model.ub=ones(size(Model.rxns))*1000;
Model.ub(94:176)=0;

% specify diffusion bounds (all non enzymatic reactions)
Model.ub(find(cellfun(@isempty,strfind(Chlamy_splitted.rxns,'Diffusion'))==0)) = 0;
Model.ub(find(cellfun(@isempty,strfind(Chlamy_splitted.rxns,'Diffusion RuBP'))==0)) = 1000;
Model.ub(find(cellfun(@isempty,strfind(Chlamy_splitted.rxns,'Diffusion 3PGA'))==0)) = 1000;
Model.ub(find(cellfun(@isempty,strfind(Chlamy_splitted.rxns,'Diffusion CO2'))==0)) = 1000;
Model.ub(219:226) = 1000;

Model.lb(strcmp(Model.rxns,'Diffusion CO2 c-hs')) = 398;
Model.lb(strcmp(Model.rxns,'Diffusion CO2 c-hs_rev')) = 0;
Model.ub(strcmp(Model.rxns,'Diffusion CO2 c-hs')) = 398;
Model.ub(strcmp(Model.rxns,'Diffusion CO2 c-hs_rev')) = 0;

Model.ub(209) = 796; % FNR
Model.ub(210) = 1194; % ATPase
Model.lb(209) = 796; % FNR
Model.lb(210) = 1194; % ATPase

Sol = optimizeCbModel(Model);

if ~isempty(Sol.x) 
   disp(size(samplesMinFlux))
%     samplesMinFlux(:,end+1) = Sol.x;
    
    options.nFiles = 5;
    options.nPointsPerFile = 100;
    options.nPointsReturned = 100;

    [Model_reduced,samples]=sampleCbModel(Model,'SampleModelsLow','ACHR',options);

    [MissingRxns,idxMR]=setdiff(Model.rxns,Model_reduced.rxns);

    samplesMinFlux(size(Model.rxns,1),end+1:end+size(samples,2)) = 0;
    samplesMinFlux(setdiff(1:size(Model.rxns,1),idxMR),end-(size(samples,2)-1):end) = samples;
    samplesMinFlux = unique(samplesMinFlux','rows')';
    
end
end
