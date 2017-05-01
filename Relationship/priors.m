function priors = sensorPriors
% returns a structure of prior parameters for the MMPP code

priors.Bsigma = 10;	%X^B sigma
priors.Dmu = zeros(1,7);	% day effect mu
priors.Dsigma = zeros(1,7)+10;	% day effect sigma
priors.Hmu = zeros(4*24,7);	% ToD effect mu
priors.Hsigma = zeros(4*24,7)+10;	% ToD effect sigma

priors.z01 = .01*10000; priors.z00 = .99*10000;     % z(t) event process
priors.z10 = .5*10000; priors.z11 = .5*10000;     

priors.Esigma = 10;       % gamma(t), or lomax, for event # process
priors.sigma0 = 10;       % gamma(t), or lomax, for event # process
%larger aE is better, smaller bE is better

priors.MODE = 0;