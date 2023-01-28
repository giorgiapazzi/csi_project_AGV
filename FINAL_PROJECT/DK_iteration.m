% Inizializzazione: 
init

global rp w_rp


s = tf('s');
J = get_linearization();
% Sistema nominale
A = J.A;
B = J.B;
C = J.C;
D = J.D;
SYS = ss(A,B,C,D);   
Gnom = minreal(tf(SYS));
[Anom Bnom Cnom Dnom] = ssdata(Gnom);
sys = minreal(ss(Anom,Bnom,Cnom,Dnom));
% Calcolo della pseudoinversa sinistra
G_inv = inv(sys'*sys)*sys';

omega = logspace(-1,6,302);

% Sistema incerto
A_i = J.A_i;
B_i = J.B_i;
Gp = C*(s*eye(5)-A_i)^(-1)*B_i; 

%%%%%%%%%%%%%%%%%
%Scelta dei pesi
%%%%%%%%%%%%%%%%%
Garray = usample(Gp,50);
orderWt = 2;
Garrayg = frd(Garray,logspace(-3,3,60));
[Usys,Info] = ucover(Garrayg,Gp.NominalValue,orderWt,'in');
wi = Info.W1;
sigma(G_inv*(Gp.NominalValue-Garray),{10^-3,10^3},'r--'); hold on; sigma(wi,{10^-3,10^3},'b-'); 
grid on; hold off;

%Funzione peso che sta sopra i valori singolari da 10^-5 in poi 
%muRSinf = 0.1322; muNPinf = 0.0729; muRPinf = 0.2841
%Non funziona con la dk
% wi = 1/(1+s*10^10)^2*1/(1+s*10^6)^2*(1+s*10^3)^4*(1+s*10^17)^10*1/(1+s*10^15)^10;

Wi = [wi wi;wi wi]; % funzione peso per l'incertezza piena Delta (matrice 2x2)


%Peso sulla performance
M = 2;      % picco massimo di S che da prassi garantisce buoni margini di guadagno sul sistema
AP = 3;     % errore massimo a regime
wBp = 1;    % frequenza minima di banda per la performance
% wP = (s/(M)^1/2+wBp)^2/(s+wBp*(AP)^1/2)^2;    % wp per maggiore pendenza 
wP = (s/M+wBp)/(s+wBp*AP);  % peso sulla performance
WP = blkdiag(wP,wP,wP,wP);

%Peso sulla Wu
Wu = tf(1);     % peso sullo sforzo di controllo
WU = blkdiag(Wu,Wu);

%Peso sulla Wt
wBt = 1;
wT = s/(s+wBt);     % peso sul rumore di misura
WT = blkdiag(wT,wT,wT,wT);


%% Generalized plant P with Wi, Wu and Wp
systemnames = 'sys Wi WP WU WT ';
inputvar = '[udel{2}; w{4}; u{2}]';
outputvar = '[Wi ; WP ; WU; WT; -w-sys]';
input_to_sys = '[u+udel]';
input_to_Wi = '[udel]';
input_to_WP = '[-sys-w]';
input_to_WU = '[u]';
input_to_WT = '[sys]';
sysoutname = 'P';
cleanupsysic = 'yes';
sysic;
P = minreal(ss(P));

Deltai = ultidyn('Deltai',[2 2]);
DeltaP = ultidyn('DeltaP',[4 10]);
Delta = blkdiag(Deltai,DeltaP);
                         

%% DK-iteration tramite musyn
% Il comando musyn prende la mixed-mu M in ingresso, sapendo che M = lft(delta,N)
% dove qui al posto della N si ha la P
nmeas = 4; nu = 2;  
omega = logspace(-1,6,302);
Fu = lft(Deltai,P);

opts = musynOptions('Display','full','MaxIter',100,'TolPerf',0.001,'FrequencyGrid',omega)
[K_DK,CLPperf,info_mu] = musyn(Fu,nmeas,nu,opts);

%Funzione di trasferimento a ciclo chiuso
looptrans = loopsens(sys,K_DK);
%Plot
figure(1)
sigma(looptrans.So, 'b'); hold on; sigma(1/WP,'b-.'); 
legend('S','gopt/W1','location','south'); hold off;
figure(2)
sigma(looptrans.To,'r'); hold on; sigma(1/WT,'r-.'); 
legend('T','gopt/WT'); hold off;
figure(3)
sigma(K_DK*looptrans.So,'g'); hold on; sigma(1/WU,'g-.'); 
legend('KS','gopt/WU'); hold off; 
% figure(4)
% sigma(S,'b',K*S,'r',T,'g',gopt/WP,'b-.',gopt/WU,'m-.',gopt/WT,'g-.',{1e-3,1e3})
% legend('S','KS','T','GAM/W1','GAM/W2','GAM/W3','Location','SouthWest');



%% Per la simulazione su SIMULINK: dk.slx
[A_DK B_DK C_DK D_DK] = ssdata(K_DK);


%% RS con robuststab
looptranfer = loopsens(Gp, K_DK);
Ti = looptranfer.Ti;
Tif = ufrd(Ti, omega);
opt = robopt('Display','on');
N = lft(P,K_DK);
F = lft(Delta,N);   % per performance robusta
Ff = ufrd(F,omega);
M = lft(Deltai,N);  % per stabilità robusta
Mf = ufrd(M,omega);

% Robusta stabilità
% stabmarg indica gli upper e lower bound
% destabunc indica esattamente l'incertezza massima tollerabile dal
% sistema oltre la quale non è robustamente stabile
%[stabmarg, destabunc, report] = robuststab(Tif,opt) 
[stabmarg, destabunc, report] = robstab(Mf,opt) 

% Robusta performance
[stabmarg, destabunc, info] = robustperf(Ff, opt)