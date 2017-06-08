% created 18/10/05
echo on
clc

%   ************************************************
%   * Simple Switching Kalman Filter demonstration *
%   ************************************************
%
%   First, we generate some random data with different dynamics.
%
%   We choose three AR(2) processes with different characteristics,
%   and for each sample from the model:
%
%      y_t = alpha_1.y_{t-1} + alpha_2.y{t-2} + e_t,
%
%   where var(e) = epsilon.
%
pause  % Press a key to continue.

alpha_a = [.1 -.9];
epsilon_a = 1;
alpha_b = [.8 .1];
epsilon_b = 1;
alpha_c = [1.8 -1];
epsilon_c = .1;

T = 200;
p = 2;

echo off

y_a = randn(T,1)*sqrt(epsilon_a);
for t=p+1:T
  y_a(t) = alpha_a*y_a(t-1:-1:t-2) + randn*sqrt(epsilon_a);
end

y_b = randn(T,1)*sqrt(epsilon_b);
for t=p+1:T
  y_b(t) = alpha_b*y_b(t-1:-1:t-2) + randn*sqrt(epsilon_b);
end

y_c = randn(T,1)*sqrt(epsilon_c);
for t=p+1:T
  y_c(t) = alpha_c*y_c(t-1:-1:t-2) + randn*sqrt(epsilon_c);
end

h = figure;
subplot(1,3,1)
plot(y_a)
subplot(1,3,2)
plot(y_b)
subplot(1,3,3)
plot(y_c)
echo on

%
%
%  Next, concatenate them and add some noise...
%
pause  % Press a key to continue.

y = [y_a; y_b; y_c];
noisevar = 2;
y = y + randn(size(y))*sqrt(noisevar);

clf
plot(y)
echo on
clc

%  Now we represent these dynamic processes in vector state space form,
%
%    x_t = A.x_{t-1} + beta_t
%    y_t = H.x_t + gamma_t
%
%    where cov(beta) = Q, and cov(gamma) = R.
%
%  We also set the other SKF parameters: the probability of transitions
%  between dynamics, the initial estimate of the hidden state, and the
%  initial probability of being in each dynamic (or 'switch state').
%
pause  % Press a key to continue.

A = zeros(p,p,3);
Q = zeros(p,p,3);
R = zeros(1,1,3);
H = zeros(1,p,3);

A(:,:,1) = [alpha_a;1 0];
H(:,:,1) = [1 0];
Q(:,:,1) = [epsilon_a 0;0 0];
R(:,:,1) = noisevar;

A(:,:,2) = [alpha_b;1 0];
H(:,:,2) = [1 0];
Q(:,:,2) = [epsilon_b 0;0 0];
R(:,:,2) = noisevar;

A(:,:,3) = [alpha_c;1 0];
H(:,:,3) = [1 0];
Q(:,:,3) = [epsilon_c 0;0 0];
R(:,:,3) = noisevar;

Z = [.99 .01 .01;.01 .99 .01;.01 .01 .99];
x_0 = [0 0]';
pi = [1 1 1]/3;


%
%
%  With this state space representation, the Switching Kalman Filter
%  can now be used to infer the probability that each data point was
%  produced by any of these three dynamics, given the data points
%  preceeding it.
%
pause  % Press a key to continue.

S = skf(y,A,H,Q,R,x_0,Z,pi);

subplot(2,1,1)
plot(y)
title('Test data, with multiple dynamic states')
subplot(2,1,2)
[a i] = max(S,[],2);
plot(i,'r+');
set(gca,'YTick',[1 2 3],'YLim',[.9 3.1]);
title('Most probable state at each time step, inferred by SKF')
clc

%  The SKF returns the probabilities of each switch state over time. This
%  plot in red shows the most likely switch state at each time point. The
%  model can be seen to differentiate between the dynamics.
%  Well done, everyone.
%
pause   % Press any key to end the demonstration.

echo off
close(h)