close all

Jo = get_linearization_Hinf();    % matrices of the linearized system
s = tf('s');
G1 = Jo.C*(s*eye(6)-Jo.A)^(-1)*Jo.B;    % transfer function
G = minreal(G1);

% Definizione dei parametri delle matrici peso.
M1 = 1.5;
wB1 = 0.399;
A1 = 1e-4;

% Matrici di peso
Wu = eye(2);
wP1 = makeweight(30,0.09,0.1);
wP2 = makeweight(30,0.75,0.5);
wP3 = makeweight(30,0.6,0.5);
wT1 = makeweight(0.1,0.5,30);
wT3 = makeweight(0.1,30,30);
wT = blkdiag(wT1,wT3);
WP =  [wP1 0;0 wP3];
% bodemag(WP1); hold on; bodemag(wP3);
% wP2 = (s/M1+wB2)/(s+wB2*A1);
% wP3 = (s/M1+wB3)/(s+wB3*A1);
% wP4 = (s/M1+wB4)/(s+wB4*A1);
% WP =  blkdiag(wP1, wP2, wP3, wP4);

% Controllore H inf
[khinf,ghinf,gopt] = mixsyn(G,WP,Wu,wT);


% Valori singolari
K = khinf;
S = inv(eye(2)+G*K);
% sigma(S)
log_vars.K = K;


figure(2);
bodemag(S); hold on; bodemag(1/WP);
% Risposta al gradino con r = [1 0]'o r = [0 1]'
% figure(2); step(eye(2)-S)
% 
% % Risposta del sistema con riferimenti r = [1 -1]'
% r = [1 -1]';
% 
% Tr = (eye(2)-S)*r;
% figure(3); step(Tr)
save('dataset','log_vars');
















