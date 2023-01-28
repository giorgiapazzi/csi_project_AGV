clc
%Inizializzazione
init

load('dataset.mat');    % dataset loading
J = get_linearization();% matrices of the linearized system
s = tf('s');
A1 = J.A;
B1 = J.B;
C1 = J.C;
D1 = J.D;
G_act = blkdiag(G_tau_phi.NominalValue,G_tau_psi.NominalValue);
G = C1*(s*eye(5)-A1)^(-1)*B1*G_act;    % transfer function
G = minreal(G);
[A, B, C, D] = ssdata(G)
% Model dimensions:
p = size(C,1); % no. of outputs (y)
[n,m] = size(B); % no. of states and inputs (u)
[A_act,B_act,C_act,D_act] = ssdata((ss(G_act)));
log_vars.A_act = A_act; 
log_vars.B_act = B_act;
log_vars.C_act = C_act; 
log_vars.D_act = D_act; 
% 1) Design state feedback regulator
Q = eye(n); % weight on integrated error
R = eye(m); % input weight
Kr = lqr(A,B,Q,R); % optimal state-feedback regulator

% 2) Design Kalman filter % don’t estimate integrator states
Bnoise = [0 0;0 0;1 0;0 0;0 1;0 0; 0 0;0 0;0 0]; % process noise model (Gd)
W = log_vars.W; % process noise weight
V = log_vars.V; % measurement noise weight
Estss = ss(A,[B Bnoise],C,0); 
[Kess, Ke] = kalman(Estss,W,V); % Kalman filter
Kr_sim = Kr(:,1:5);
save('dataset.mat','log_vars')
%K_lqg =[A-B*Kr-Ke*C, Ke;
%         -Kr, zeros(2,4)];
%K_LQG =ss(A-B*Kr-Ke*C,Ke,-Kr,zeros(2,4));
% Test simulation on simulink: 
% For linear system: LQG_no_integrator.slx
% For non linear system: LQG_non_linear.slx