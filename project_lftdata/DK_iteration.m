close all
global rp w_rp

s = tf('s');
J = get_linearization();    % richiama la funzione per la linearizzazione del sistema
A = J.A;
B = J.B;
A_i = J.A_i;    % Matrice della dinamica A incerta
B_i = J.B_i;    % Matrice B incerta
C = J.C;
D = J.D;

%Costruzione della funzione di trasferimento incerta
%Gp = C*(s*eye(5)-A_i)^(-1)*B_i;
Gp = ss(A_i,B_i,C,D);
G_i= minreal(tf(Gp));
[Ap Bp Cp Dp] = ssdata(G_i);
Gi = minreal(ss(Ap,Bp,Cp,Dp));

%Definizione dei pesi di performance
% Costruzione della WP, peso sulla performance
M = 1;  % picco massimo di S che da prassi garantisce buoni margini di guadagno sul sistema
AP = 10^-1; % errore massimo a regime
wBp = 1;    % frequenza minima di banda per la performance
wP = 5*10^-2*(s/M+wBp)/(s+wBp*AP);  % peso sulla performance
%wP = (s/(M)^1/2+wBp)^2/(s+wBp*(AP)^1/2)^2; % wp per maggiore pendenza 
WP = blkdiag(wP,wP,wP,wP);

%Costruzione della WU, peso sullo sforzo di controllo
wBu = 1;
WU = 1/10*tf(eye(2));   % peso sullo sforzo di controllo

%Costruzione della WT, peso sul rumore di misura
Wt = makeweight(10^-2,20,500);
WT = blkdiag(Wt,Wt,Wt,Wt);


%Funzione di trasferimento del sistema nominale
SYS = ss(A,B,C,D);   
Gnom = minreal(tf(SYS));
[Anom Bnom Cnom Dnom] = ssdata(Gnom);
sys = minreal(ss(Anom,Bnom,Cnom,Dnom));


%% Generalized plant P with WP, WU and WT
systemnames = 'Gp WP WU WT';
inputvar = '[w{4}; u{2}]';
outputvar = '[WP ; WU; WT; -w-Gp]';
input_to_Gp = '[u]';
input_to_WP = '[-w-Gp]';
input_to_WU = '[u]';
input_to_WT = '[Gp]';
sysoutname = 'P_i';
cleanupsysic = 'yes';
sysic;

[P, Delta, blk] = lftdata(P_i);
                          

%% DK-iteration tramite musyn
     
% Il comando musyn prende la mixed-mu M in ingresso, sapendo che M = lft(delta,N)
% dove qui al posto della N si ha la P
nmeas = 4; nu = 2;  
omega = logspace(-3,3,61);
opts = musynOptions('Display','full','MaxIter',20,'TolPerf',0.001,'FrequencyGrid',omega)
[K_DK,CLPperf,info_mu] = musyn(P_i,nmeas,nu,opts);
[A_DK B_DK C_DK D_DK] = ssdata(K_DK);

%Verifica della nominale stabilità
N = lft(P,K_DK);
eig(N);


%% DK ITERATION MANUALE 
% funzione di interpolazione scelta di ordine 2
%con la funzione wi di grado 1 funziona bene fino alla quarta iterazione
%poi si perde ma arriva a muRP<1

% omega = logspace(-3,3,61);
% nmeas = 4; nu = 2; d0 = 1; 
% %delta in questo caso è diag{delta_i, delta_p}
% %delta_i è un blocco diagonale 2x2 ed è per questo che ho [1 1; 1 1];
% %delta_P invece è una matrice piena (non diagonale)
% D_left = append(tf(eye(14)));
% D_right = append(d0,d0,d0,d0,d0,d0,d0,d0,d0,tf(eye(6)));
% %
% % START ITERATION.
% %
% % STEP 1: Find H-infinity optimal controller
% % with given scalings:
% %
% 
%     [K,Nsc,gamma,info] = hinfsyn(P_i,nmeas,nu,....
%                    'method','lmi','Tolgam',1e-3);
%     
% 
%     Nf = frd(lft(P,K),omega);
% %
% gamma_prec = gamma+1; 
% gamma_corr = gamma;
% N_it = 0;
% while (N_it<10)
% % STEP 2: Compute mu using upper bound:
%     %Verifica della robusta stabilità
%     [mubnds,Info] = mussv(Nf(1:9,1:9),blk,'c'); 
%     bodemag(mubnds(1,1),omega);
%     murs = norm(mubnds(1,1),inf,1e-6);
%     %Verifica della performance nominale
%     [mubnds_pn,Info_np] = mussv(Nf(10:end,10:end),[4 10],'c');
%     bodemag(mubnds_pn(1,1),omega);
%     munp = norm(mubnds_pn(1,1),inf,1e-6);
%     %Verifica della robusta performance
%     [mubnds_rp,Info_rp] = mussv(Nf,[9 0;4 10],'c');
%     bodemag(mubnds_rp(1,1),omega);
%     murp = norm(mubnds_rp(1,1),inf,1e-6)
% %   
% % STEP 3: Fit resulting D-scales:
% %
%     [dsysl,dsysr] = mussvunwrap(Info_rp);
%     dsysl = dsysl/dsysl(3,3);
%     func_order_4 = fitfrd(genphase(dsysl(1,1)),1); 
%     %viene generata la fase interpolando con una funzione del 4° ordine
%     %func_order_4=func_order_4.C*(inv(s*eye(4)-func_order_4.A))*func_order_4.B+func_order_4.D; 
%     % poiché viene restituita in forma di stato viene trasfromata in 
%     % funzione di trasferimento prima di metterla in Dk
%     d0 = tf(minreal(func_order_4));
%     
% %     func_order_4_p = fitfrd(genphase(dsysl_p(1,1)),4);
% %     func_order_4_p=func_order_4_p.C*(inv(s*eye(4)-func_order_4_p.A))*func_order_4_p.B+func_order_4_p.D; 
% %     D_right=func_order_4_p;
%     D_left = append(d0,d0,d0,d0,d0,d0,d0,d0,d0,tf(eye(14)));
%     D_right = append(d0,d0,d0,d0,d0,d0,d0,d0,d0,tf(eye(6)));
%     
%      [K,Nsc,gamma,info] = hinfsyn(D_left*P*inv(D_right),nmeas,nu,....
%                    'method','lmi','Tolgam',1e-3);
% 
%     Nf = frd(lft(P,K),omega);
% 
% %     gamma_prec = gamma_corr;
% %     gamma_corr = gamma;
%     N_it = N_it+1;
%   
% end

%% RS con robuststab
%altrimenti is può far applicare la robusta stabilità direttamente da
%matlab, 
% looptranfer = loopsens(Gp, K_DK);
% Ti = looptranfer.Ti;
% Tif = ufrd(Ti, omega);
% opt = robopt('Display','on');
% con il comando successivo mi dice sulla robusta stabilità
%in particolare stabmarg mi indica gli upper e lower bound, l'inverso del
%lower bound deve essere uguale al muRSinf ottenuto con il comando mussv
%destabunc mi indica esattamente l'incertezza massima tollerabile dal
%sistema oltre la quale non è robustamente stabile
% [stabmarg, destabunc, report] = robuststab(Tif,opt) %incertezze massime che mi porterebbero all'instabilità
