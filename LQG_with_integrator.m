clc
load('dataset.mat');
J = get_linearization();
s = tf('s');
G = J.C*(s*eye(6)-J.A)^(-1)*J.B
% Model dimensions:
p = size(J.C,1) % no. of outputs (y)
[n,m] = size(J.B) % no. of states and inputs (u)
Znp=zeros(n,p); Zpn=zeros(p,n);
Znn=zeros(n,n); Zpp=zeros(p,p);

% 1) Design state feedback regulator
A = [J.A Znp;-J.C Zpp]; 
B = [J.B;-J.D]; % augment plant with integrator
C = [J.C Zpp];

Q=[Znn Znp;Zpn eye(p,p)]; % weight on integrated error
R=eye(m); % input weight
rank(ctrb(A,B))
rank(obsv(A,C))
eig(A)
[Ao, Bo, Co] = obsvf(A,B,C)
[Abar,Bbar,Cbar] = ctrbf(A,B,C)
Kr=lqr(A,B,Q,R); % LQR
Krp=Kr(1:m,1:n); % state feedback
Kri=Kr(1:m,n+1:n+p); % integrator feedback

% 2) Design Kalman filter % don’t estimate integrator states
Bnoise = [0 0;0 0;0 0;1 0;0 0;0 1]; % process noise model (Gd)
W = log_vars.W; % process noise weight
V = log_vars.V; % measurement noise weight
Estss = ss(J.A,[J.B Bnoise],J.C,0); 
[Kess, Ke] = kalman(Estss,W,V); % Kalman filter
%Ke = lqe(J.A,Bnoise,J.C,W,V); % Kalman filter gain
% 
% % 3) Form overall controller
% Ac=[Zmm Zmn;-b*Kri a-b*Krp-Ke*c]; % integrator included
% Bcr = [eye(m); Znm]; Bcy = [-eye(m); Ke];
% Cc = [-Kri -Krp]; Dcr = Zmm; Dcy = Zmm;
% Klqg2 = ss(Ac,[Bcr Bcy],Cc,[Dcr Dcy]); % Final 2 dof controller from [r y]' to u
% Klqg = ss(Ac,-Bcy,Cc,-Dcy); % feedback part of controller from -y to u
