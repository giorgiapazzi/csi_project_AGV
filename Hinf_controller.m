close all

Jo = get_linearization();    % matrices of the linearized system
s = tf('s');
G1 = Jo.C*(s*eye(6)-Jo.A)^(-1)*Jo.B;    % transfer function
G = minreal(G1);

% Definizione dei parametri delle matrici peso.
M1 = 1.5;

% Matrici di peso
Wu = 0.5.*eye(2);
WP1 = makeweight(70,0.399,0.5);
WP2 = makeweight(100,0.537,0.5);
WP3 = makeweight(45,0.781,0.5);
WP4 = makeweight(70,2.1,0.5);
WP =  blkdiag(wP1, wP2, wP3, wP4);
% wP1 = (s/M1+wB1)/(s+wB1*A1);
% wP2 = (s/M1+wB2)/(s+wB2*A1);
% wP3 = (s/M1+wB3)/(s+wB3*A1);
% wP4 = (s/M1+wB4)/(s+wB4*A1);
% WP =  blkdiag(wP1, wP2, wP3, wP4);

% Controllore H inf
[khinf,ghinf,gopt] = mixsyn(G,WP,Wu,[]);


% Valori singolari
K = khinf;
S = inv(eye(2)+G*K);
sigma(S)

% Risposta al gradino con r = [1 0]'o r = [0 1]'
figure(2); step(eye(2)-S)

% Risposta del sistema con riferimenti r = [1 -1]'
r = [1 -1]';

Tr = (eye(2)-S)*r;
figure(3); step(Tr)
















