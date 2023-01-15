%% Nominal plant
s = tf('s');
Jo = get_linearization_Hinf();    % matrices of the linearized system
s = tf('s');
Gnom = Jo.C*(s*eye(6)-Jo.A)^(-1)*Jo.B;    % transfer function
G = minreal(Gnom);


%% Controller
load('dataset');
K = log_vars.K;
[A B C D] = ssdata(K);
% Ks = C*(s*eye(6)-A)^(-1)*B;
% Ks = minreal(Ks);
[num1,den1] = ss2tf(A,B,C,D,1);
[num2,den2] = ss2tf(A,B,C,D,2);
K1 = tf(num1(1,:),den1);
K2 = tf(num1(2,:),den1);
K3 = tf(num2(1,:),den2);
K4 = tf(num2(2,:),den2);
Kt = [K1 K3;K2 K4];
Kt = minreal(Kt);

%% Weighting filter for uncertainty modelling
w1 = makeweight(0.20,35,10);
w2 = makeweight(0.25,40,10);
Wi = blkdiag(w1,w2);

wP = 0.04*(s+10)/(s+0.005);
WP = blkdiag(wP,wP);

wu = 4e-2*(0.01*s+1)/(0.005*s+1);
Wu = blkdiag(wu,wu);

%% Generalized plant P
systemnames = 'G WP Wu Wi';
inputvar = '[ydel{2}; ref{2}; w{2}; u{2}]';
outputvar = '[Wi ; WP ; Wu ; ref-G-w]';
input_to_G = '[u+ydel]';
input_to_WP = '[ref-G-w]';
input_to_Wu = '[u]';
input_to_Wi = '[u]';
sysoutname = 'P';
cleanupsysic = 'yes';
sysic;
