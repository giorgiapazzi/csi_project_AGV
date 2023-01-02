clc
load('dataset.mat');
J = get_linearization();
s = tf('s');
G = J.C*(s*eye(6)-J.A)^(-1)*J.B
% Model dimensions:
p = size(J.C,1) % no. of outputs (y)
[n,m] = size(J.B) % no. of states and inputs (u)


% 1) Design state feedback regulator
Q=[eye(n)]; % weight on integrated error
R=eye(m); % input weight
Kr=lqr(J.A,J.B,Q,R); % optimal state-feedback regulator
%Krp=Kr(1:m,1:n); % state feedback

% 2) Design Kalman filter % don’t estimate integrator states
Bnoise = eye(n); % process noise model (Gd)
W = log_vars.W; % process noise weight
V = log_vars.V; % measurement noise weight
Estss = ss(J.A,[J.B Bnoise],J.C,[J.D zeros(4,6)]); 
[Kess, Ke] = kalman(Estss,W,V); % Kalman filter
%Ke = lqe(J.A,Bnoise,J.C,W,V); % Kalman filter gain


% 3) Form overall controller
Klqg=[J.A-J.B*Kr-Ke*J.C, Ke;
      -Kr, zeros(2,4)]; % integrator included

% % Simulation
% sys1 = feedback(G,Klqg); step(sys1,50); % 1 dof simulation
% sys = feedback(G*Klqg2,1,2,1,+1); % 2 dof simulation
% sys2 = sys*[1; 0]; hold; step(sys2,50);

% % 3) Form controller from [r y]’ to u
% Ac=[0mm 0mn;-b*Kri a-b*Krp-Ke*c]; % integrators included
% Bc = [eye(m) -eye(m);0nm +Ke];
% Cc = [-Kri -Krp]; Dc = [0mm 0mm];
% Klqg = ss(Ac,Bc,Cc,Dc); % Final 2-DOF controller. Includes
% % integrators, Kr and Kalman filte