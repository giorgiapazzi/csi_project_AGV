clc

load('dataset.mat');    % dataset loading
J = get_linearization_lqg();    % matrices of the linearized system
s = tf('s');
G = J.C*(s*eye(6)-J.A)^(-1)*J.B;    % transfer function

% Model dimensions:
p = size(J.C,1); % no. of outputs (y)
[n,m] = size(J.B); % no. of states and inputs (u)
Znp=zeros(n,p); Zpn=zeros(p,n);
Znn=zeros(n,n); Zpp=zeros(p,p);

% 1) Design state feedback regulator
% augment plant with integrator:
A = [J.A Znp;-J.C Zpp]; 
B = [J.B;-J.D]; 
C = [J.C Zpp];
Q = [Znn Znp;Zpn eye(p,p)]; % weight on integrated error
R = eye(m); % input weight

% % Controllo delle proprietà strutturali del sistema aumentato
% Gi = C*(s*eye(10)-A)^(-1)*B; % sistema di trasferimento del sistema aumentato
% rank(ctrb(A,B));    % rango della matrice di raggiungibilità
% rank(obsv(A,C));    % rango della matrice di osservabilità
% eig(A); % autovalori della matrice A del sistema aumentato
% [Ao, Bo, Co] = obsvf(A,B,C);    % sistema in forma standard di osservabilità
% [Abar,Bbar,Cbar,T] = ctrbf(A,B,C) % sistema in forma standard di raggiungibilità

Kr=lqr(A,B,Q,R); % LQR
Krp=Kr(1:m,1:n); % state feedback
Kri=Kr(1:m,n+1:n+p); % integrator feedback

% 2) Design Kalman filter % don’t estimate integrator states
Bnoise = [0 0;0 0;0 0;1 0;0 0;0 1]; % process noise model (Gd)
W = log_vars.W; % process noise weight
V = log_vars.V; % measurement noise weight
Estss = ss(J.A,[J.B Bnoise],J.C,0); 
[Kess, Ke] = kalman(Estss,W,V); % Kalman filter

% Test simulation on simulink: model_lqg_int.slx