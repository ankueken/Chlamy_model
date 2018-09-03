function xprime = chlamy_ode(t,x,kinetic_param,S)
%
% This function computes the conservation of mass for all metabolites and
% enzymes in the form of free and complex (PMID:24928774).
%
%
% INPUTS:
% -------
%    t_interval: The time interval for simulations
%             x: Vector of initial concentrations of metabolites and enzymes
% kinetic_param: Vector of elementary kinetic parameters
%             S: Stoichiometric matrix of decomposed reactions
%
%
% OUTPUTS:
% ---------
%        xprime: Value of dx/dt
%

xprime = zeros(length(x),1);

% Reaction rates
v = kinetic_param .* prod(repmat(x,1,length(kinetic_param)) .^ (abs( S .* (S<0)) ))';

% Loop over each metabolite
xprime = S*v;
