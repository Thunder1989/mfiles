function priors = sensorPriors
% returns a structure of prior parameters for the MMPP code

priors.aL = 1; priors.bL = 1;       % lambda0, baseline rate
priors.aD = zeros(1,7)+10;           % day effect dirichlet params
priors.aH = zeros(4*24,7)+1;          % time of day effect dirichlet params

priors.z01 = .01*10000; priors.z00 = .99*10000;     % z(t) event process
priors.z10 = .5*10000; priors.z11 = .5*10000;     
priors.aE = 2; priors.bE = 1/2;       % gamma(t), or lomax, for event # process
%larger aE is better, smaller bE is better
priors.MODE = 0;