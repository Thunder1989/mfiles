function priors = sensorPriors
% returns a structure of prior parameters for the MMPP code

priors.aL = 1; priors.bL = 1;       % lambda0, baseline rate
priors.aD = zeros(1,7)+5;           % day effect dirichlet params
priors.aH = zeros(4*24,7)+1;          % time of day effect dirichlet params

priors.z01 = .01*10000; priors.z00 = .99*10000;     % z(t) event process
priors.z10 = .1*10000; priors.z11 = .9*10000;     
priors.aE = 3; priors.bE = 1/4;       % gamma(t), or NBin, for event # process

priors.MODE = 0;