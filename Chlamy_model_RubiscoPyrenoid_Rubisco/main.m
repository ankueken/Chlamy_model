addpath(genpath('..'))
load('../Chlamy_model.mat') % load model
load('ExpDataProteinDistr_0.mat') % load experimental data on enzyme distributions

% sample steady-state flux distributions - include threshold on relative enzyme activity
samples_l_00 = sampling_sst_f(Chlamy_splitted, ExpDataProteinDistr,0,'lowCO2');
samples_h_00 = sampling_sst_f(Chlamy_splitted, ExpDataProteinDistr,0,'highCO2');

save('samples_0.mat','samples_l_00','samples_h_00')

% load('samples_0.mat')
% samples with realistic rubisco distribution
l_idx=find((samples_l_00(93,:)./samples_l_00(7,:))>1.5);
samples_l_00=samples_l_00(:,l_idx);
h_idx=find((samples_h_00(93,:)./samples_h_00(7,:))<0.3);
samples_h_00=samples_h_00(:,h_idx);

% keep unique solutions
samples_l_00 = unique(samples_l_00','rows')';
samples_h_00 = unique(samples_h_00','rows')';

%%
load('ExpDataConcentrations_0.mat') % load experimental data on metabolites

[X_l_00,K_l_00,V_l_00,CORR_l_00,PVAL_l_00,absolute_difference_l_00,chi_square_individual_l_00,chi_square_l_00,mi_l_00,ma_l_00] = ...
    getParameterValues(Chlamy_splitted,ExpDataProteinDistr,ExpDataConcentrations,samples_l_00,1:size(samples_l_00,2),'lowCO2');

[X_h_00,K_h_00,V_h_00,CORR_h_00,PVAL_h_00,absolute_difference_h_00,chi_square_individual_h_00,chi_square_h_00,mi_h_00,ma_h_00] = ...
    getParameterValues(Chlamy_splitted,ExpDataProteinDistr,ExpDataConcentrations,samples_h_00,1:size(samples_h_00,2),'highCO2');

save('Result_get_param_original.mat')

%%
% load('Result_get_param_original.mat')

[X_l_00s,K_l_00s,V_l_00s,CORR_l_00s,PVAL_l_00s,absolute_difference_l_00s,chi_square_individual_l_00s,chi_square_l_00s] = ...
    sortResults(X_l_00,K_l_00,V_l_00,CORR_l_00,PVAL_l_00,absolute_difference_l_00,chi_square_individual_l_00,chi_square_l_00);

[X_h_00s,K_h_00s,V_h_00s,CORR_h_00s,PVAL_h_00s,absolute_difference_h_00s,chi_square_individual_h_00s,chi_square_h_00s] = ...
    sortResults(X_h_00,K_h_00,V_h_00,CORR_h_00,PVAL_h_00,absolute_difference_h_00,chi_square_individual_h_00,chi_square_h_00);

X_l_00s(X_l_00s<=1e-20)=0;X_h_00s(X_h_00s<=1e-20)=0;

%% chi square
chi_square_h_00s=sum(chi_square_individual_l_00s([1 2 4:7 9:11 15],1:5000));
chi_square_l_00s=sum(chi_square_individual_h_00s([1 2 4:7 9:11 15],1:5000));

mean(chi_square_l_00s(1:1000))
mean(chi_square_h_00s(1:1000))

sum(chi_square_l_00s(1:1000)<=7.962)/1000
sum(chi_square_h_00s(1:1000)<=7.962)/1000

%% individual chi square

subplot(2,1,1)
boxplot(chi_square_individual_l_00s(:,1:1000)','labels',ExpDataConcentrations{1,:}')
set(gca,'XTickLabelRotation',45)
ylabel('Chi-square')
title('LC*')
subplot(2,1,2)
boxplot(chi_square_individual_h_00s(:,1:1000)','labels',ExpDataConcentrations{1,:}')
set(gca,'XTickLabelRotation',45)
ylabel('Chi-square')
title('HC')

%% correlation

subplot(1,2,1)
[H,b]=hist(CORR_l_00s(1:1000));
bar(b,H/sum(H),1)
xlabel('Pearson correlation coefficient')
ylabel('Fraction')
title('LC*')
subplot(1,2,2)
[H,b]=hist(CORR_h_00s(1:1000));
bar(b,H/sum(H),1)
xlabel('Pearson correlation coefficient')
title('HC')

%% net flux
Net_samples_l_00 = Vnet_calc(Chlamy_splitted, V_l_00s(:,1:1000));
Net_samples_h_00 = Vnet_calc(Chlamy_splitted, V_h_00s(:,1:1000));

%% concentration
Ms = cell(size(ExpDataConcentrations,2),1); % metabolites stroma + pyrenoid
Mcs = cell(size(ExpDataConcentrations,2),1); % metabolite-enzyme-complexes stroma + pyrenoid
My = cell(size(ExpDataConcentrations,2),1); % metabolites stroma + pyrenoid
Mcy = cell(size(ExpDataConcentrations,2),1); % metabolite-enzyme-complexes stroma + pyrenoid
for i=1:size(ExpDataConcentrations,2)
    m = find(cellfun(@isempty,strfind(Chlamy_splitted.mets, strcat(char(ExpDataConcentrations{1,i}),' [s]')))==0);
    Ms{i}(end+1:end+length(m),1) = m;
    m = find(cellfun(@isempty,strfind(Chlamy_splitted.mets, strcat(char(ExpDataConcentrations{1,i}),' [y]')))==0);
    My{i}(end+1:end+length(m),1) = m;
    m = intersect(find(cellfun(@isempty,strfind(Chlamy_splitted.mets, char(ExpDataConcentrations{1,i})))==0),find(cellfun(@isempty,strfind(Chlamy_splitted.mets, 'complex [s]'))==0));
    Mcs{i}(end+1:end+length(m),1) = m;
    m = intersect(find(cellfun(@isempty,strfind(Chlamy_splitted.mets, char(ExpDataConcentrations{1,i})))==0),find(cellfun(@isempty,strfind(Chlamy_splitted.mets, 'complex [y]'))==0));
    Mcy{i}(end+1:end+length(m),1) = m;
    if i==15
        m1 = find(cellfun(@isempty,strfind(Chlamy_splitted.mets, 'X5P [s]'))==0);
        m2 = find(cellfun(@isempty,strfind(Chlamy_splitted.mets, 'Ru5P [s]'))==0);
        Ms{i}(end+1:end+length([m1;m2]),1) = [m1;m2];
        m1 = find(cellfun(@isempty,strfind(Chlamy_splitted.mets, 'X5P [y]'))==0);
        m2 = find(cellfun(@isempty,strfind(Chlamy_splitted.mets, 'Ru5P [y]'))==0);
        My{i}(end+1:end+length([m1;m2]),1) = [m1;m2];
        m1 = intersect(find(cellfun(@isempty,strfind(Chlamy_splitted.mets, 'X5P'))==0),find(cellfun(@isempty,strfind(Chlamy_splitted.mets, 'complex [s]'))==0));
        m2 = intersect(find(cellfun(@isempty,strfind(Chlamy_splitted.mets, 'Ru5P'))==0),find(cellfun(@isempty,strfind(Chlamy_splitted.mets, 'complex [s]'))==0));
        Mcs{i}(end+1:end+length([m1;m2]),1) = [m1;m2];
        m1 = intersect(find(cellfun(@isempty,strfind(Chlamy_splitted.mets, 'X5P'))==0),find(cellfun(@isempty,strfind(Chlamy_splitted.mets, 'complex [y]'))==0));
        m2 = intersect(find(cellfun(@isempty,strfind(Chlamy_splitted.mets, 'Ru5P'))==0),find(cellfun(@isempty,strfind(Chlamy_splitted.mets, 'complex [y]'))==0));
        Mcy{i}(end+1:end+length([m1;m2]),1) = [m1;m2];
    end
end

for i=1:15
    PData_l(i,:) = [mean(sum(X_l_00s(Ms{i},1:1000),1)); mean(sum(X_l_00s(Mcs{i},1:1000),1));mean(sum(X_l_00s(My{i},1:1000),1)); mean(sum(X_l_00s(Mcy{i},1:1000),1))];
    PDatas_l(i,:) = [std(sum(X_l_00s(Ms{i},1:1000),1)); std(sum(X_l_00s(Mcs{i},1:1000),1));std(sum(X_l_00s(My{i},1:1000),1)); std(sum(X_l_00s(Mcy{i},1:1000),1))];
    Tot_C_l(i,:) = [mean(sum(X_l_00s([Ms{i};Mcs{i};My{i};Mcy{i}],1:1000),1)) std(sum(X_l_00s([Ms{i};Mcs{i};My{i};Mcy{i}],1:1000),1))];
end
for i=1:15
    PData_h(i,:) = [mean(sum(X_h_00s(Ms{i},1:1000),1)); mean(sum(X_h_00s(Mcs{i},1:1000),1));mean(sum(X_h_00s(My{i},1:1000),1)); mean(sum(X_h_00s(Mcy{i},1:1000),1))];
    PDatas_h(i,:) = [std(sum(X_h_00s(Ms{i},1:1000),1)); std(sum(X_h_00s(Mcs{i},1:1000),1));std(sum(X_h_00s(My{i},1:1000),1)); std(sum(X_h_00s(Mcy{i},1:1000),1))];
    Tot_C_h(i,:) = [mean(sum(X_h_00s([Ms{i};Mcs{i};My{i};Mcy{i}],1:1000),1)) std(sum(X_h_00s([Ms{i};Mcs{i};My{i};Mcy{i}],1:1000),1))];
end

figure
MAP = [0.95 0.87 0.73;0.45 0.26 0.26;0.93 0.84 0.84;0.42 0.25 0.39];
subplot(2,1,1)
hold on
m=[1 2 4 5 6 7 9 10 11 15];
bar(1:10,[Tot_C_h(m,1)'; str2double(ExpDataConcentrations{3,m})]','grouped')
errorbar([1:10]-0.15,Tot_C_h(m,1),Tot_C_h(m,2),'k.')
errorbar([1:10]+0.15,str2double(ExpDataConcentrations{3,m})',str2double(ExpDataConcentrations{4,m}),'k.')
colormap(MAP)
set(gca,'XTick',1:10,'XTickLabel',ExpDataConcentrations{1,m},'XTickLabelRotation',45,'YScale', 'log')
ylabel('Concentration [\muM]')
legend('Model', 'Experiment','Orientation','horizontal','Location','Best')
title('HC')
subplot(2,1,2)
hold on
m=[1 2 4 5 6 7 9 10 11 15];
bar(1:10,[Tot_C_l(m,1)'; str2double(ExpDataConcentrations{6,m})]','grouped')
errorbar([1:10]-0.15,Tot_C_l(m,1),Tot_C_l(m,2),'k.')
errorbar([1:10]+0.15,str2double(ExpDataConcentrations{6,m})',str2double(ExpDataConcentrations{7,m}),'k.')
colormap(MAP)
set(gca,'XTick',1:10,'XTickLabel',ExpDataConcentrations{1,m},'XTickLabelRotation',45,'YScale', 'log')
ylabel('Concentration [\muM]')
title('LC*')

%% enzymes
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

for i=1:size(ExpDataProteinDistr,2)-1
    EData_l(i,:) = [mean(sum(X_l_00s(TE_s{i},1:1000),1)); mean(sum(X_l_00s(TEc_s{i},1:1000),1));mean(sum(X_l_00s(TE_y{i},1:1000),1)); mean(sum(X_l_00s(TEc_y{i},1:1000),1))];
    EDatas_l(i,:) = [std(sum(X_l_00s(TE_s{i},1:1000),1)); std(sum(X_l_00s(TEc_s{i},1:1000),1));std(sum(X_l_00s(TE_y{i},1:1000),1)); std(sum(X_l_00s(TEc_y{i},1:1000),1))];
    Tot_E_l(i,:) = [mean(sum(X_l_00s([TE_s{i};TEc_s{i};TE_y{i};TEc_y{i}],1:1000),1)) std(sum(X_l_00s([TE_s{i};TEc_s{i};TE_y{i};TEc_y{i}],1:1000),1))];
    Ratio_E_l(i,1) = mean(sum([sum(X_l_00s(TE_y{i},1:1000),1)],1) ./ sum([sum(X_l_00s(TE_s{i},1:1000),1)],1));
    Ratio_E_l(i,2) = std(sum([sum(X_l_00s(TE_y{i},1:1000),1); sum(X_l_00s(TEc_y{i},1:1000),1)],1) ./ sum([sum(X_l_00s(TE_s{i},1:1000),1); sum(X_l_00s(TEc_s{i},1:1000),1)],1));
end
for i=1:size(ExpDataProteinDistr,2)-1
    EData_h(i,:) = [mean(sum(X_h_00s(TE_s{i},1:1000),1)); mean(sum(X_h_00s(TEc_s{i},1:1000),1));mean(sum(X_h_00s(TE_y{i},1:1000),1)); mean(sum(X_h_00s(TEc_y{i},1:1000),1))];
    EDatas_h(i,:) = [std(sum(X_h_00s(TE_s{i},1:1000),1)); std(sum(X_h_00s(TEc_s{i},1:1000),1));std(sum(X_h_00s(TE_y{i},1:1000),1)); std(sum(X_h_00s(TEc_y{i},1:1000),1))];
    Tot_E_h(i,:) = [mean(sum(X_h_00s([TE_s{i};TEc_s{i};TE_y{i};TEc_y{i}],1:1000),1)) std(sum(X_h_00s([TE_s{i};TEc_s{i};TE_y{i};TEc_y{i}],1:1000),1))];
    Ratio_E_h(i,1) = mean(sum([sum(X_h_00s(TE_y{i},1:1000),1); sum(X_h_00s(TEc_y{i},1:1000),1)],1) ./ sum([sum(X_h_00s(TE_s{i},1:1000),1); sum(X_h_00s(TEc_s{i},1:1000),1)],1));
    Ratio_E_h(i,2) = std(sum([sum(X_h_00s(TE_y{i},1:1000),1); sum(X_h_00s(TEc_y{i},1:1000),1)],1) ./ sum([sum(X_h_00s(TE_s{i},1:1000),1); sum(X_h_00s(TEc_s{i},1:1000),1)],1));
end
i=1
mean(sum([sum(X_l_00s(TE_y{i},1:1000),1); sum(X_l_00s(TEc_y{i},1:1000),1)],1) ./ sum([sum(X_l_00s(TE_s{i},1:1000),1); sum(X_l_00s(TEc_s{i},1:1000),1); sum(X_l_00s(TE_y{i},1:1000),1); sum(X_l_00s(TEc_y{i},1:1000),1)],1))
mean(sum([sum(X_h_00s(TE_y{i},1:1000),1); sum(X_h_00s(TEc_y{i},1:1000),1)],1) ./ sum([sum(X_h_00s(TE_s{i},1:1000),1); sum(X_h_00s(TEc_s{i},1:1000),1); sum(X_h_00s(TE_y{i},1:1000),1); sum(X_h_00s(TEc_y{i},1:1000),1)],1))

mean(Net_samples_l_00.flux(20,1:1000)./ sum([Net_samples_l_00.flux(3,1:1000); Net_samples_l_00.flux(20,1:1000)],1))
mean(Net_samples_h_00.flux(20,:)./ sum([Net_samples_h_00.flux(3,:); Net_samples_h_00.flux(20,:)],1))

table(ExpDataProteinDistr{1,2:end}', Ratio_E_l, Ratio_E_h)

%% kcat
k_id_s = [7 16 22 31 34 46 49 55 58 64 10];
k_id_y = [93 102 108 117 120 132 135 141 144 150 96];
kcat = [5.8 223 3372 64.5 0.42 25.2 81 0.0148 0.1382 3440 84.7];

%% 
X_l_00s(:,1:1000); % concentration
V_l_00s(:,1:1000); % flux 
K_l_00s(:,1:1000); % kcat

table([mean(K_h_00s(k_id_s([10 9 8 7 5 4 6 3 2 11 1]),1:1000)')';...
   mean(K_h_00s(k_id_y(1),1:1000))],[std(K_h_00s(k_id_s([10 9 8 7 5 4 6 3 2 11 1]),1:1000)')';...
   std(K_h_00s(k_id_y(1),1:1000))],[mean(K_l_00s(k_id_s([10 9 8 7 5 4 6 3 2 11 1]),1:1000)')';...
   mean(K_l_00s(k_id_y(1),1:1000))],[std(K_l_00s(k_id_s([10 9 8 7 5 4 6 3 2 11 1]),1:1000)')';...
   std(K_l_00s(k_id_y(1),1:1000))])

load('../free_bound_met_indices.mat')
for i=1:length(Msf)
    ACsfl(i,:) = X_l_00s([Msf{i}],1:1000);
    ACyfl(i,:) = X_l_00s([Myf{i}],1:1000);
    ACsbl(i,:) = sum(X_l_00s([Msb{i}],1:1000),1);
    ACybl(i,:) = sum(X_l_00s([Myb{i}],1:1000),1);
end
load('../free_bound_met_indices.mat')
for i=1:length(Msf)
    ACsfh(i,:) = X_h_00s([Msf{i}],1:1000);
    ACyfh(i,:) = X_h_00s([Myf{i}],1:1000);
    ACsbh(i,:) = sum(X_h_00s([Msb{i}],1:1000),1);
    ACybh(i,:) = sum(X_h_00s([Myb{i}],1:1000),1);
end

N = {'RuBP';'3PGA';'BPGA';'GAP';'DHAP';'FBP';'F6P';'E4P';'Xu5P';'SBP';'S7P';'R5P';'Ru5P';'G1P';'G6P';'ADPG';'ATP';'ADP';'NADP';'NADPH';'CO2'};
MAP = [1 0.84 0; 0.87 0.49 0;0.68 0.92 1;0.2 0.3 0.49];
i=1
subplot(2,4,1)
bar([mean(ACsfh(i,:)) mean(ACyfh(i,:))])
hold on
errorbar(1,mean(ACsfh(i,:)),std(ACsfh(i,:)), std(ACsfh(i,:)),'k.')
errorbar(2,mean(ACyfh(i,:)),std(ACyfh(i,:)), std(ACyfh(i,:)),'k.')
ylabel('free RuBP [\muM]')
set(gca, 'XTickLabel',{'stroma','pyrenoid'})
subplot(2,4,2)
bar([mean(ACsfl(i,:)) mean(ACyfl(i,:))])
hold on
errorbar(1,mean(ACsfl(i,:)),std(ACsfh(i,:)), std(ACsfl(i,:)),'k.')
errorbar(2,mean(ACyfl(i,:)),std(ACyfh(i,:)), std(ACyfl(i,:)),'k.')
set(gca, 'XTickLabel',{'stroma','pyrenoid'})

subplot(2,4,3)
bar([mean(ACsbh(i,:)) mean(ACybh(i,:))])
hold on
errorbar(1,mean(ACsbh(i,:)),std(ACsbh(i,:)), std(ACsbh(i,:)),'k.')
errorbar(2,mean(ACybh(i,:)),std(ACybh(i,:)), std(ACybh(i,:)),'k.')
ylabel('bound RuBP [\muM]')
set(gca, 'XTickLabel',{'stroma','pyrenoid'})
subplot(2,4,4)
bar([mean(ACsbl(i,:)) mean(ACybl(i,:))])
hold on
errorbar(1,mean(ACsbl(i,:)),std(ACsbl(i,:)), std(ACsbl(i,:)),'k.')
errorbar(2,mean(ACybl(i,:)),std(ACybl(i,:)), std(ACybl(i,:)),'k.')
set(gca, 'XTickLabel',{'stroma','pyrenoid'})

i=2
subplot(2,4,5)
bar([mean(ACsfh(i,:)) mean(ACyfh(i,:))])
hold on
errorbar(1,mean(ACsfh(i,:)),std(ACsfh(i,:)), std(ACsfh(i,:)),'k.')
errorbar(2,mean(ACyfh(i,:)),std(ACyfh(i,:)), std(ACyfh(i,:)),'k.')
ylabel('free 3PGA [\muM]')
set(gca, 'XTickLabel',{'stroma','pyrenoid'})
subplot(2,4,6)
bar([mean(ACsfl(i,:)) mean(ACyfl(i,:))])
hold on
errorbar(1,mean(ACsfl(i,:)),std(ACsfh(i,:)), std(ACsfl(i,:)),'k.')
errorbar(2,mean(ACyfl(i,:)),std(ACyfh(i,:)), std(ACyfl(i,:)),'k.')
set(gca, 'XTickLabel',{'stroma','pyrenoid'})

subplot(2,4,7)
bar([mean(ACsbh(i,:)) mean(ACybh(i,:))])
hold on
errorbar(1,mean(ACsbh(i,:)),std(ACsbh(i,:)), std(ACsbh(i,:)),'k.')
errorbar(2,mean(ACybh(i,:)),std(ACybh(i,:)), std(ACybh(i,:)),'k.')
ylabel('bound 3PGA [\muM]')
set(gca, 'XTickLabel',{'stroma','pyrenoid'})
subplot(2,4,8)
bar([mean(ACsbl(i,:)) mean(ACybl(i,:))])
hold on
errorbar(1,mean(ACsbl(i,:)),std(ACsbl(i,:)), std(ACsbl(i,:)),'k.')
errorbar(2,mean(ACybl(i,:)),std(ACybl(i,:)), std(ACybl(i,:)),'k.')
set(gca, 'XTickLabel',{'stroma','pyrenoid'})

bar([mean(Net_samples_h_00.flux([2 37 38],:)')' mean(Net_samples_l_00.flux([2 37 38],:)')'])
hold on
errorbar([1:3]-0.15,mean(Net_samples_h_00.flux([2 37 38],:)')',std(Net_samples_h_00.flux([2 37 38],:)')','k.')
errorbar([1:3]+0.15,mean(Net_samples_l_00.flux([2 37 38],:)')',std(Net_samples_l_00.flux([2 37 38],:)')','k.')
MAP = [0.87 0.49 0;0.2 0.3 0.49];
colormap(MAP)
legend('HC', 'LC*','Orientation','horizontal')
set(gca,'XTickLabel',{'CO_2';'RuBP';'3PGA'})
ylabel({'Net flux of metabolite transport'; 'from stroma into pyrenoid'})


%% delta G - free and bound concentration
% a A + b B --> c C + d D
% Keq from equilibrator
% delta G = -RT * (ln Keq - ln [C][D]/[A][B])

clear X_l X_h
Names = readtable('../Escher_met_list.csv');
for i=1:size(Names,1)-1
    XN = Names{i+1,2:end};
    X_l(i,:) = sum(X_l_00s(XN(~isnan(XN)),1:1000),1);
    X_h(i,:) = sum(X_h_00s(XN(~isnan(XN)),1:1000),1);
end

RT = 8.3144598*298.15;

%% stroma - low
% Rubisco
deltaG2_ls(1,:) = -RT*(log(1.3*10e4) - log(X_l(2,:).^2./(X_l(21,:).*X_l(1,:)))); 
% PGK
deltaG2_ls(2,:) = -RT*(log(5.8*10e-4) - log(X_l(3,:).*X_l(18,:)./(X_l(17,:).*X_l(2,:)))); 
% GAPDH
deltaG2_ls(3,:) = -RT*(log(34.4) - log(X_l(4,:).*X_l(19,:)./(X_l(3,:).*X_l(20,:))));
% TPI
deltaG2_ls(4,:) = -RT*(log(9.12) - log(X_l(5,:)./X_l(4,:)));
% FBA
deltaG2_ls(5,:) = -RT*(log(1.3*10e3) - log(X_l(6,:)./(X_l(4,:).*X_l(5,:))));
% SBA
deltaG2_ls(6,:) = -RT*(log(145) - log(X_l(10,:)./(X_l(5,:).*X_l(8,:))));
% FBPase
deltaG2_ls(7,:) = -RT*(log(101) - log(X_l(7,:)./X_l(6,:)));
% SBPase
deltaG2_ls(8,:) = -RT*(log(1.5*10e3) - log(X_l(11,:)./X_l(10,:)));
% TRK 1
deltaG2_ls(9,:) = -RT*(log(0.0177) - log((X_l(9,:).*X_l(8,:))./(X_l(4,:).*X_l(7,:))));
% TRK 2
deltaG2_ls(10,:) = -RT*(log(90.216) - log((X_l(9,:).*X_l(12,:))./(X_l(4,:).*X_l(11,:))));
% RPE
deltaG2_ls(11,:) = -RT*(log(0.257) - log(X_l(13,:)./X_l(9,:)));
% RPI
deltaG2_ls(12,:) = -RT*(log(0.458) - log(X_l(13,:)./X_l(12,:)));
% PRK
deltaG2_ls(13,:) = -RT*(log(2*10e3) - log((X_l(18,:).*X_l(1,:))./(X_l(17,:).*X_l(13,:))));

%% pyrenoid - low
% Rubisco
deltaG2_lp(1,:) = -RT*(log(1.3*10e4) - log(X_l(23,:).^2./(X_l(22,:).*X_l(42,:))));

%% stroma - high
% Rubisco
deltaG2_hs(1,:) = -RT*(log(1.3*10e4) - log(X_h(2,:).^2./(X_h(21,:).*X_h(1,:)))); 
% PGK
deltaG2_hs(2,:) = -RT*(log(5.8*10e-4) - log(X_h(3,:).*X_h(18,:)./(X_h(17,:).*X_h(2,:)))); 
% GAPDH
deltaG2_hs(3,:) = -RT*(log(34.4) - log(X_h(4,:).*X_h(19,:)./(X_h(3,:).*X_h(20,:))));
% TPI
deltaG2_hs(4,:) = -RT*(log(9.12) - log(X_h(5,:)./X_h(4,:)));
% FBA
deltaG2_hs(5,:) = -RT*(log(1.3*10e3) - log(X_h(6,:)./(X_h(4,:).*X_h(5,:))));
% SBA
deltaG2_hs(6,:) = -RT*(log(145) - log(X_h(10,:)./(X_h(5,:).*X_h(8,:))));
% FBPase
deltaG2_hs(7,:) = -RT*(log(101) - log(X_h(7,:)./X_h(6,:)));
% SBPase
deltaG2_hs(8,:) = -RT*(log(1.5*10e3) - log(X_h(11,:)./X_h(10,:)));
% TRK 1
deltaG2_hs(9,:) = -RT*(log(0.0177) - log((X_h(9,:).*X_h(8,:))./(X_h(4,:).*X_h(7,:))));
% TRK 2
deltaG2_hs(10,:) = -RT*(log(90.216) - log((X_h(9,:).*X_h(12,:))./(X_h(4,:).*X_h(11,:))));
% RPE
deltaG2_hs(11,:) = -RT*(log(0.257) - log(X_h(13,:)./X_h(9,:)));
% RPI
deltaG2_hs(12,:) = -RT*(log(0.458) - log(X_h(13,:)./X_h(12,:)));
% PRK
deltaG2_hs(13,:) = -RT*(log(2*10e3) - log((X_h(18,:).*X_h(1,:))./(X_h(17,:).*X_h(13,:))));

%% pyrenoid -high
% Rubisco
deltaG2_hp(1,:) = -RT*(log(1.3*10e4) - log(X_h(23,:).^2./(X_h(22,:).*X_h(42,:))));

enzymes={'Rubisco_p','Rubisco_s','PGK','GAPDH','TPI','FBA','SBA','FBPase','SBPase','TRK1','TRK2','RPE','RPI','PRK'};

figure
subplot(4,4,1)
[H,b]=hist(deltaG2_lp);
bar(b,H/11)
colormap([0.5 0.5 0.5])
title(enzymes(1))

for i=2:length(enzymes)
    subplot(4,4,i)
    [H,b]=hist(deltaG2_ls(i-1,:));
    bar(b,H/11)
    title(enzymes(i))
end
