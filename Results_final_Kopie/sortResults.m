function [Xn,Kn,Vn,CORRn,PVALn,absolute_difference_n,chi_square_individual_n,chi_square_n] = sortResults(X,K,V,CORR,PVAL,absolute_difference,chi_square_individual,chi_square)

[~,ind] = sort(chi_square);

Xn = X(:,ind);
Kn = K(:,ind);
Vn = V(:,ind);
CORRn = CORR(ind);
PVALn = PVAL(ind);
absolute_difference_n = absolute_difference(:,ind);
chi_square_individual_n = chi_square_individual(:,ind);
chi_square_n = chi_square(ind);

end
