addpath(genpath('..'))
load('../Chlamy_model.mat') % load model
load('ExpDataProteinDistr_t.mat') % load experimental data on enzyme distributions

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
load('ExpDataConcentrations_0.mat')

[X_l_00,K_l_00,V_l_00,CORR_l_00,PVAL_l_00,absolute_difference_l_00,chi_square_individual_l_00,chi_square_l_00,mi_l_00,ma_l_00] = ...
    getParameterValues(Chlamy_splitted,ExpDataProteinDistr,ExpDataConcentrations,samples_l_00,1:size(samples_l_00,2),'lowCO2');

[X_h_00,K_h_00,V_h_00,CORR_h_00,PVAL_h_00,absolute_difference_h_00,chi_square_individual_h_00,chi_square_h_00,mi_h_00,ma_h_00] = ...
    getParameterValues(Chlamy_splitted,ExpDataProteinDistr,ExpDataConcentrations,samples_h_00,1:size(samples_h_00,2),'highCO2');

save('Result_get_param_original.mat')
%%
%load('Result_get_param_original.mat')
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

