% Uses the Robust Control toolbox
close all
global rp w_rp

load('dataset');
s = tf('s');
sys = log_vars.sys;
S = log_vars.S;
T = log_vars.T;
K = log_vars.K;
J = get_linearization();
A = J.A;
B = J.B;
A_i = J.A_i;
B_i = J.B_i;
C = J.C;
D = J.D;
% autovalori = sigma(T);
% sigma(T)
omega = logspace(-1,6,302);
%
% Weights.
%
% Wp = 0.5*tf([10 1],[10 1.e-5])*eye(2);
rp_tau = w_rp/(rp);
wi = rp_tau*rp*s/(1+rp*s);
Wi = blkdiag(wi,wi);
WP = log_vars.WP;
WU = log_vars.WU;


%% Generalized plant P with Wi, Wu and Wp
systemnames = 'sys WP WU Wi';
inputvar = '[udel{2}; w{4}; u{2}]';
outputvar = '[Wi ;WU; WP ; -w-sys]';
input_to_sys = '[u+udel]';
input_to_WP = '[sys]';
input_to_WU = '[u]';
input_to_Wi = '[u]';
sysoutname = 'P';
cleanupsysic = 'yes';
sysic;
P = minreal(ss(P));

%
%% Initialize.
%
omega = logspace(-3,3,61);
blk = [1 1; 1 1; 1 1; 1 1; 1 1; 1 1];
nmeas = 4; nu = 2; d0 = 1; 
%delta in questo caso è diag{delta_i, delta_p}
%delta_i è un blocco diagonale 2x2 ed è per questo che ho [1 1; 1 1];
%delta_P invece è una matrice piena (non diagonale)
D_left = append(d0,d0,d0*tf(eye(6)),tf(eye(4)));
D_right = append(d0,d0,d0*tf(eye(6)));
%
% START ITERATION.
%
% STEP 1: Find H-infinity optimal controller
% with given scalings:
%

    [K,Nsc,gamma,info] = hinfsyn(D_left*P*inv(D_right),nmeas,nu,....
                   'method','lmi','Tolgam',1e-3);

    Nf = frd(lft(P,K),omega);
%
gamma_prec = gamma+1; 
gamma_corr = gamma;
while (gamma_prec > gamma_corr)
% STEP 2: Compute mu using upper bound:
    %Verifica della robusta stabilità
    [mubnds,Info] = mussv(Nf(1:2,1:2),[1 1;1 1],'c'); 
    bodemag(mubnds(1,1),omega);
    murs = norm(mubnds(1,1),inf,1e-6);
    %Verifica della performance nominale
    [mubnds_pn,Info_np] = mussv(Nf(3:8,3:6),[4 6],'c');
    bodemag(mubnds_pn(1,1),omega);
    munp = norm(mubnds_pn(1,1),inf,1e-6);
    %Verifica della robusta performance
    [mubnds_rp,Info_rp] = mussv(Nf,[1 1;1 1;4 6],'c');
    bodemag(mubnds_rp(1,1),omega);
    murp = norm(mubnds_rp(1,1),inf,1e-6);
%   
% STEP 3: Fit resulting D-scales:
%
    [dsysl,dsysr] = mussvunwrap(Info_rp);
    dsysl = dsysl/dsysl(3,3);
    func_order_4 = fitfrd(genphase(dsysl(1,1)),2); 
    %viene generata la fase interpolando con una funzione del 4° ordine
%     func_order_4=func_order_4.C*(inv(s*eye(4)-func_order_4.A))*func_order_4.B+func_order_4.D; 
    % poiché viene restituita in forma di stato viene trasfromata in 
    % funzione di trasferimento prima di metterla in Dk
    d0 = tf(minreal(func_order_4));
    
%     func_order_4_p = fitfrd(genphase(dsysl_p(1,1)),4);
%     func_order_4_p=func_order_4_p.C*(inv(s*eye(4)-func_order_4_p.A))*func_order_4_p.B+func_order_4_p.D; 
%     D_right=func_order_4_p;
    D_left = append(d0,d0,d0*tf(eye(6)),tf(eye(4)));
    D_right = append(d0,d0,d0*tf(eye(6)));
    
     [K,Nsc,gamma,info] = hinfsyn(D_left*P*inv(D_right),nmeas,nu,....
                   'method','lmi','Tolgam',1e-3);

    Nf = frd(lft(P,K),omega);[K,Nsc,gamma,info] = hinfsyn(D_left*P*inv(D_right),nmeas,nu,....
                   'method','lmi','Tolgam',1e-3);

    Nf = frd(lft(P,K),omega);

    gamma_prec = gamma_corr;
    gamma_corr = gamma;
end