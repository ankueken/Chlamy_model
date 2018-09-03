function [xt,Diff,v,T,X] = solve_ode_test(model,time,x0,kinetic_param)

% This function solves the ODEs 
%
%
% INPUTS:
% -------
%                model: model structure that containts the information of the model
%                tSpan: simulation time
%                   x0: The initial concentration of the metabolites and enzymes
%        kinetic_param: kinetic parameters of rxns
% 
%
% OUTPUTS:
% ------
%           conc: metabolite concentrations and enzyme fractions:
%                 Columns correspond to the time points 
%                 Rows correspond to metabolites 
%              v: Rate of rxns in the model at each time point:
%                 Columns correspond to time points
%                 Row correspond to reactions
%

S = model.S;

options=odeset('NonNegative',1:size(S,1));   
  % Solve system of ODEs
  [T,X]=ode15s(@(t,x)chlamy_ode(t,x,kinetic_param,S),time,x0,options);
%   X(1,:)=x0;
%   X(2,:)=chlamy_ode(1,x0,kinetic_param,S);
  
  % check for steady-state
%  Diff = sum(abs((X(2,:)+eps)./(X(1,:)+eps)));
   Diff = sum((abs(X(1,:)-X(end,:))));
   xt=X(1,:)';
%  %disp(max((abs(X(2,:)))))
% 
%   %---------  Concentrations ---------------
%   % Rows in conc correspond to metabolites and columns correspond to time points
%   conc = (X(1,:)+X(2,:))'; 
  
  %---------- Compute reaction rates ----------------
  % Total number of time points
%    time_point_num=length(T);
%   
%    v = zeros(length(kinetic_param),time_point_num);
  
  % First compute the rate of reactions 
v = kinetic_param .* prod(repmat(xt,1,length(kinetic_param)) .^ (abs( S .* (S<0)) ))';
v=v';
end
