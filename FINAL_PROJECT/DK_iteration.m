%Inizializzazione
init
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
WU = 1/20*tf(eye(2));   % peso sullo sforzo di controllo

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
     
nmeas = 4; nu = 2;  % numero delle uscite e degli ingressi dell'impianto
omega = logspace(-3,3,61);
opts = musynOptions('Display','full','MaxIter',20,'TolPerf',0.001,'FrequencyGrid',omega)
[K_DK,CLPperf,info_mu] = musyn(P_i,nmeas,nu,opts);  % myxsin si applica direttamente all'impianto incerto

% Per validazione in simulink: dk.slx
[A_DK B_DK C_DK D_DK] = ssdata(K_DK);


%% Analisi delle prestazioni
%Verifica della nominale stabilità
N = lft(P,K_DK);
eig(N) % Gli autovalori di N devono essere tutti a parte reale negativa
% In alternativa si possono controllare gli autovalori a ciclo chiuso
CL = feedback(sys*K_DK,eye(4));
eig(CL)    % Gli autovalori a ciclo chiuso devono essere tutti a parte reale negativa

% RS con robstab e RP con robperf
opt = robopt('Display','on');
Delta_P = ultidyn('Delta_P', [4 10]); % Incertezze fittizie
delta = blkdiag(Delta, Delta_P);    % Matrice delle incertezze finale

F = lft(delta, N);
Ff = ufrd(F, omega);
M = lft(Delta,N);   % Si ottiene la fdt tra udel e ydel
Mf = ufrd(M,omega); % risposta in frequenza

% Robusta stabilità
[stabmarg, destabunc, info] = robstab(Mf, opt)

% Robusta prestazione
[stabmarg, destabunc, info] = robustperf(Ff, opt)

