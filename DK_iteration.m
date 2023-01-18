% Uses the Robust Control toolbox

%
% Weights.
%
% Wp = 0.5*tf([10 1],[10 1.e-5])*eye(2);
wi = 10^4*1/(s)^5/(1+10^5*s)*(1+s*10^2)^6;
Wi = blkdiag(wi,wi);
WP = log_vars.WP;
%
systemnames = 'sys WP WU Wi';
inputvar = '[udel{2}; ref{4}; u{2}]';
outputvar = '[Wi ; WP ; WU ; ref-sys]';
input_to_sys = '[u+udel]';
input_to_WP = '[sys]';
input_to_WU = '[u]';
input_to_Wi = '[u]';
sysoutname = 'P';
cleanupsysic = 'yes';
sysic;
P = minreal(ss(P));
%
% Initialize.
%
omega = logspace(-3,3,61);
blk = [1 1; 1 1; 2 2];
nmeas = 2; nu = 2; d0 = 1;
%delta in questo caso è diag{delta_i, delta_p}
%delta_i è un blocco diagonale 2x2 ed è per questo che ho [1 1; 1 1];
%delta_P invece è una matrice piena (non diagonale)
D = append(d0,d0,tf(eye(2)),tf(eye(2)));
%
% START ITERATION.
%
% STEP 1: Find H-infinity optimal controller
% with given scalings:
%
[K,Nsc,gamma,info] = hinfsyn(D*P*inv(D),nmeas,nu,....
                   'method','lmi','Tolgam',1e-3);
Nf = frd(lft(P,K),omega);
%
% STEP 2: Compute mu using upper bound:
%
[mubnds,Info] = mussv(Nf,blk,'c'); 
bodemag(mubnds(1,1),omega);
murp = norm(mubnds(1,1),inf,1e-6);
%
% STEP 3: Fit resulting D-scales:
%
[dsysl,dsysr] = mussvunwrap(Info);
dsysl = dsysl/dsysl(3,3);
d1 = fitfrd(genphase(dsysl(1,1)),4);