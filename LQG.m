clc

load('dataset.mat');    % dataset loading
J = get_linearization();% matrices of the linearized system
s = tf('s');
A = J.A;
B = J.B;
C = J.C;
D = J.D;
G = J.C*(s*eye(5)-J.A)^(-1)*J.B;    % transfer function

% Model dimensions:
p = size(J.C,1); % no. of outputs (y)
[n,m] = size(J.B); % no. of states and inputs (u)


% 1) Design state feedback regulator
Q = eye(n); % weight on integrated error
R = eye(m); % input weight
Kr = lqr(J.A,J.B,Q,R); % optimal state-feedback regulator

% 2) Design Kalman filter % don’t estimate integrator states
Bnoise = [0 0;0 0;1 0;0 0;0 1]; % process noise model (Gd)
W = log_vars.W; % process noise weight
V = log_vars.V; % measurement noise weight
Estss = ss(J.A,[J.B Bnoise],J.C,0); 
[Kess, Ke] = kalman(Estss,W,V); % Kalman filter

K_lqg =[A-B*Kr-Ke*C, Ke;
        -Kr, zeros(2,4)];
K_LQG =ss(A-B*Kr-Ke*C,Ke,-Kr,zeros(2,4));
% Test simulation on simulink: model.slx