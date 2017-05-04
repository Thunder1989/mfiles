function priors = sensorPriors
% returns a structure of prior parameters for the MMPP code

priors.Bsigma = 10;	%X^B sigma

%day prior hyperparameter
priors.Dmu = zeros(1,7);	
priors.Dsigma = zeros(1,7)+10;

%ToD prior hyperparameter
priors.Hmu = zeros(4*24,7);
priors.Hsigma = zeros(4*24,7)+10;

%transition prior hyperparameter
priors.z01 = .01*10000; priors.z00 = .99*10000;     % z(t) event process
priors.z10 = .25*10000; priors.z11 = .75*10000;     

priors.Esigma = 10;
%NE prior hyperparameter
priors.mu0 = 0;
priors.sigma0 = 10;

priors.MODE = 0;