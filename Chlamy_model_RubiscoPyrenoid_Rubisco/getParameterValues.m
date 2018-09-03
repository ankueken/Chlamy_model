function [X,K,V,CORR,PVAL,absolute_difference,chi_square_individual,chi_square,miC,maC] = ...
    getParameterValues(Chlamy_splitted,ExpDataProteinDistr,ExpDataConcentrations,samplesMinFlux,idx,cond)

if strcmp(cond,'lowCO2')
    idx_c = 6;
elseif strcmp(cond,'highCO2')
    idx_c = 3;
end

X=[];K=zeros(length(Chlamy_splitted.rxns),0);V=[];chi_square=[];absolute_difference=[];CORR=[];chi_square_individual=[];PVAL=[];
miC=[];maC=[];
k_id_s = [7 16 22 31 34 46 49 55 58 64 10];
k_id_y = [93 102 108 117 120 132 135 141 144 150 96];
kcat = [5.8 223 3372 64.5 0.42 25.2 81 0.0148 0.1382 3440 84.7];

M = cell(size(ExpDataConcentrations,2),1); % metabolites stroma + pyrenoid
Mc = cell(size(ExpDataConcentrations,2),1); % metabolite-enzyme-complexes stroma + pyrenoid
for i=1:size(ExpDataConcentrations,2)
    m = find(cellfun(@isempty,strfind(Chlamy_splitted.mets, strcat(char(ExpDataConcentrations{1,i}),' [s]')))==0);
    M{i}(end+1:end+length(m),1) = m;
    m = find(cellfun(@isempty,strfind(Chlamy_splitted.mets, strcat(char(ExpDataConcentrations{1,i}),' [y]')))==0);
    M{i}(end+1:end+length(m),1) = m;
    m = intersect(find(cellfun(@isempty,strfind(Chlamy_splitted.mets, char(ExpDataConcentrations{1,i})))==0),find(cellfun(@isempty,strfind(Chlamy_splitted.mets, 'complex [s]'))==0));
    Mc{i}(end+1:end+length(m),1) = m;
    m = intersect(find(cellfun(@isempty,strfind(Chlamy_splitted.mets, char(ExpDataConcentrations{1,i})))==0),find(cellfun(@isempty,strfind(Chlamy_splitted.mets, 'complex [y]'))==0));
    Mc{i}(end+1:end+length(m),1) = m;
    if i==15
        m1 = find(cellfun(@isempty,strfind(Chlamy_splitted.mets, 'X5P [s]'))==0);
        m2 = find(cellfun(@isempty,strfind(Chlamy_splitted.mets, 'Ru5P [s]'))==0);
        M{i}(end+1:end+length([m1;m2]),1) = [m1;m2];
        m1 = find(cellfun(@isempty,strfind(Chlamy_splitted.mets, 'X5P [y]'))==0);
        m2 = find(cellfun(@isempty,strfind(Chlamy_splitted.mets, 'Ru5P [y]'))==0);
        M{i}(end+1:end+length([m1;m2]),1) = [m1;m2];
        m1 = intersect(find(cellfun(@isempty,strfind(Chlamy_splitted.mets, 'X5P'))==0),find(cellfun(@isempty,strfind(Chlamy_splitted.mets, 'complex [s]'))==0));
        m2 = intersect(find(cellfun(@isempty,strfind(Chlamy_splitted.mets, 'Ru5P'))==0),find(cellfun(@isempty,strfind(Chlamy_splitted.mets, 'complex [s]'))==0));
        Mc{i}(end+1:end+length([m1;m2]),1) = [m1;m2];
        m1 = intersect(find(cellfun(@isempty,strfind(Chlamy_splitted.mets, 'X5P'))==0),find(cellfun(@isempty,strfind(Chlamy_splitted.mets, 'complex [y]'))==0));
        m2 = intersect(find(cellfun(@isempty,strfind(Chlamy_splitted.mets, 'Ru5P'))==0),find(cellfun(@isempty,strfind(Chlamy_splitted.mets, 'complex [y]'))==0));
        Mc{i}(end+1:end+length([m1;m2]),1) = [m1;m2];
    end
end

TE_s = cell(size(ExpDataProteinDistr,2)-1,1); % enzymes stroma
TE_y = cell(size(ExpDataProteinDistr,2)-1,1); % enzymes pyrenoid
TEc_s = cell(size(ExpDataProteinDistr,2)-1,1); % enzymes complex stroma
TEc_y = cell(size(ExpDataProteinDistr,2)-1,1); % enzymes complex pyrenoid
for i=1:size(ExpDataProteinDistr,2)-1
    m = find(cellfun(@isempty,strfind(Chlamy_splitted.mets, strcat(char(ExpDataProteinDistr{1,i+1}),' [s]')))==0);
    TE_s {i}(end+1:end+length(m),1) = m;
    m = find(cellfun(@isempty,strfind(Chlamy_splitted.mets, strcat(char(ExpDataProteinDistr{1,i+1}),' [y]')))==0);
    TE_y {i}(end+1:end+length(m),1) = m;
    m = intersect(find(cellfun(@isempty,strfind(Chlamy_splitted.mets, char(ExpDataProteinDistr{1,i+1})))==0),find(cellfun(@isempty,strfind(Chlamy_splitted.mets, 'complex [s]'))==0));
    TEc_s {i}(end+1:end+length(m),1) = m;
    m = intersect(find(cellfun(@isempty,strfind(Chlamy_splitted.mets, char(ExpDataProteinDistr{1,i+1})))==0),find(cellfun(@isempty,strfind(Chlamy_splitted.mets, 'complex [y]'))==0));
    TEc_y {i}(end+1:end+length(m),1) = m;
end

% find reaction around enzyme
id_y=[]; id_s=[];
for i=1:size(TE_s,1)
    enzyme_consuming_rxns_y = find(Chlamy_splitted.S(TE_y{i},:)~=0); % pyrenoid
    enzyme_consuming_rxns_s = find(Chlamy_splitted.S(TE_s{i},:)~=0); % stroma
    
    id_y(end+1:end+length(enzyme_consuming_rxns_y)) = enzyme_consuming_rxns_y;
    id_s(end+1:end+length(enzyme_consuming_rxns_s)) = enzyme_consuming_rxns_s;
end

for SAMPLE = idx
    
    % remove reactions with no flux
    redModel = removeRxns(Chlamy_splitted,Chlamy_splitted.rxns(samplesMinFlux(:,SAMPLE)<1e-22),[],false);
    
    %% stoichiometry
    
    A = redModel.S';
    A(A>0) = 0;
    A = abs(A);
    alpha = A;
    
    kDiff = zeros(0,size(alpha,1));
    for i=1:length(id_y)
        is = find(strcmp(redModel.rxns,Chlamy_splitted.rxns(id_s(i))));
        iy = find(strcmp(redModel.rxns,Chlamy_splitted.rxns(id_y(i))));
        
        if ~isempty(is) && ~isempty(iy)
            kDiff(end+1,[is iy]) = [1 -1];
        end
    end
    
    Csum = zeros(size(M,1),size(alpha,2));
    for i=1:size(M,1)
        Csum(i,[M{i}' Mc{i}']) = 1;
    end
    
    Model2.S = [alpha	eye(size(alpha,1)); % alpha*x + k = v
                Csum zeros(size(Csum,1),size(alpha,1));
                zeros(size(kDiff,1),size(alpha,2)) kDiff;
                zeros(size(kDiff,1),size(alpha,2)) kDiff;];
    
    %% mets, rxns
    
    Model2.rxns=cell(size(Model2.S,2),1);
    for i=1:size(Model2.S,2)
        Model2.rxns{i} = strcat('R',num2str(i));
    end
    
    Model2.mets=cell(size(Model2.S,1),1);
    for i=1:size(Model2.S,1)
        Model2.mets{i} = strcat('M',num2str(i));
    end
    
    %% b vector
    
    posFlux = samplesMinFlux(samplesMinFlux(:,SAMPLE)>=1e-22,SAMPLE);
    if strcmp(cond,'lowCO2')
        Model2.b = [log(posFlux); log(str2double(ExpDataConcentrations{idx_c,:})'); 3.8*ones(size(kDiff,1),1); -3.8*ones(size(kDiff,1),1)];
    else
        Model2.b = [log(posFlux); log(str2double(ExpDataConcentrations{idx_c,:})'); 4.4*ones(size(kDiff,1),1); -4.4*ones(size(kDiff,1),1)];
    end
    %% csense: Constraints sense ('E' equality, 'G' greater than, 'L' less than)
    
    Model2.csense = [repmat('E',length(log(posFlux)),1); repmat('L',length(log(str2double(ExpDataConcentrations{idx_c,:})')),1); repmat('L',size(kDiff,1),1); repmat('G',size(kDiff,1),1)];
    
    %% objective
    
    % known kcat values
    for i=1:length(k_id_s)
        if ~isempty(find(strcmp(redModel.rxns,Chlamy_splitted.rxns(k_id_s(i)))))
            k_red_id_s(i) = find(strcmp(redModel.rxns,Chlamy_splitted.rxns(k_id_s(i))));
        else
            k_red_id_s(i) = NaN;
        end
        
        if ~isempty(find(strcmp(redModel.rxns,Chlamy_splitted.rxns(k_id_y(i)))))
            k_red_id_y(i) = find(strcmp(redModel.rxns,Chlamy_splitted.rxns(k_id_y(i))));
        else
            k_red_id_y(i) = NaN;
        end
    end
    
    % since lower bound is fixed - minimize distance to lb
    kc = zeros(size(alpha,1),1);
    kc([k_red_id_s(~isnan(k_red_id_s)) k_red_id_y(~isnan(k_red_id_y))]) = -1;
    
    % concentration related objective coefficients
    cc = -ones(size(alpha,2),1);
    
    for i=2:size(TE_s,1)
        cc(TE_s{i}) = -1;
        cc(TE_y{i}) = -1;
        cc(TEc_s{i}) = -1;
        cc(TEc_y{i}) = -1;
    end

    Model2.c = [cc; kc];
       
    %% bounds
    
    % bounds on concentrations
    x_lb = log(ones(size(alpha,2),1)*0.01);
    x_ub = log(ones(size(alpha,2),1)*1e5);
    
    for i=1:size(M,1)
        x_ub(M{i}) = log(str2double(ExpDataConcentrations{idx_c,i}));
        if i==2
            x_lb(M{i}) = log(str2double(ExpDataConcentrations{idx_c,i})*0.5);
        end
        if strcmp(cond,'lowCO2')
            if i==4 || i==5 || i==6 || i==7 || i==8 || i==9 || i==10 || i==11 || i==12 || i==13 || i==14
                x_lb(M{i}) = log(str2double(ExpDataConcentrations{idx_c,i})*0.99);
            end
        else
            if i==4 || i==5 || i==7 || i==8 || i==10 || i==11 || i==12 || i==13 || i==14
                x_lb(M{i}) = log(str2double(ExpDataConcentrations{idx_c,i})*0.99);
            elseif i==6
                x_lb(M{i}) = log(str2double(ExpDataConcentrations{idx_c,i})*0.95);
            elseif i==9 
                x_lb(M{i}) = log(str2double(ExpDataConcentrations{idx_c,i})*0.9);
            end
        end
    end
    
    x_lb(70:128) = log(1e-50*ones(length(70:128),1));
    x_ub(70:128) = log(1e-50*ones(length(70:128),1)); % pyrenoid metabolites
    
    noNANk = k_red_id_s(~isnan(k_red_id_s));
    % bounds on kcat values
    k_lb = log(1e-50*ones(size(alpha,1),1));
    k_lb(noNANk([1 3 4 5 6 7 8 9 10 11])) = log(kcat([1 3 4 5 6 7 8 9 10 11]));
    
    k_ub = log(ones(size(alpha,1),1)*1e4);
    k_ub(noNANk([1 3 4 6 7 10 11])) = log(kcat([1 3 4 6 7 10 11]));
       
    Model2.lb = [x_lb; k_lb];
    Model2.ub = [x_ub; k_ub];
    
    %% solve
    
    Sol=optimizeCbModel(Model2);
    
%     table([1:15]',ExpDataConcentrations{1,:}',Csum*exp(Sol.x(1:size(Csum,2))),str2double(ExpDataConcentrations{idx_c,:})')
 
    if Sol.stat==1
        
        x0 = exp(Sol.x(1:length(redModel.mets)));
        kinetic_param = exp(Sol.x(length(redModel.mets)+1:length(redModel.mets)+length(redModel.rxns)));
        
        [~,Diff,~]=solve_ode_test(redModel,0:1,x0,kinetic_param);
        if Diff<1e-6
              
            for i=1:size(M,1) % 
                Simulated(i)=sum(x0([M{i};Mc{i}]));
            end
                        
            [CORR(end+1,1),PVAL(end+1,1)] = corr(Simulated',str2double(ExpDataConcentrations{idx_c,:})');
            absolute_difference(:,end+1) = Simulated - str2double(ExpDataConcentrations{idx_c,:});
            y = (Simulated - str2double(ExpDataConcentrations{idx_c,:})).^2;
            s = str2double(ExpDataConcentrations{idx_c+1,:}).^2;
            chi_square_individual(:,end+1) = (y(s~=0)./s(s~=0))';
            chi_square(end+1,1)=sum(sum(y(s~=0)./s(s~=0)));

            X(:,end+1)=x0;
            K(samplesMinFlux(:,SAMPLE)>=1e-22,end+1)=kinetic_param;
            samplesMinFlux(samplesMinFlux(:,SAMPLE)<1e-22,SAMPLE) = 0;
            V(:,end+1)=samplesMinFlux(:,SAMPLE);
            [mi,ma] = MyCobraFVA(Model2);
            miC(:,end+1) = exp(mi);
            maC(:,end+1) = exp(ma);
            
        end
    end
end
end
