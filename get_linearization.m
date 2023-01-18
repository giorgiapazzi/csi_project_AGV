function j_matrix = get_linearization_lqg()
    clc
    global  l L IPy IPz Mv mp ma IAy d IGz IAz a b a_i b_i rp_i ra_i rp ra
   
    %% Vettore di stato X = [y theta phi_dot psi psi_dot]
    %Tau=[tau_phi tau_psi]  ingressi
    Tau=sym('tau',[2 1],'real');
    
    %X vettore di stato
    X = sym('x',[5 1],'real');

    % definizione della h (modello di osservazione)
    % vettore di variabili misurate [X(3),X(4),X(5),X(6)]
    h = [X(1),X(2),X(3),X(4)]';
    
    % definizione della f (funzione di transizione di stato)
    % Matrice della dinamica
    M=[a*rp^2*cos(X(4))^2+b*(rp/L)^2*sin(X(4))^2+IPy -IPz*(rp/L)*sin(X(4));-IPz*(rp/L)*sin(X(4)) IPz];
    N=[(rp^2*b/L^2-a*rp^2)*cos(X(4))*sin(X(4))*X(5)*X(3);-IPz*rp/L*cos(X(4))*X(5)*X(3)];
    q_ddot = inv(M)*([Tau(1); Tau(2)]-N);
    
%     M_i=[a_i*rp_i^2*cos(X(4))^2+b_i*(rp_i/L)^2*sin(X(4))^2+IPy -IPz*(rp_i/L)*sin(X(4));-IPz*(rp_i/L)*sin(X(4)) IPz];
%     N_i=[(rp_i^2*b_i/L^2-a_i*rp_i^2)*cos(X(4))*sin(X(4))*X(5)*X(3);-IPz*rp_i/L*cos(X(4))*X(5)*X(3)];
%     q_ddot_i = inv(M_i)*([Tau(1); Tau(2)]-N_i);    
    %X_dot
    %x1_dot = rp*X(4)*cos(X(3))*cos(X(5));
    x1_dot = rp*X(3)*sin(X(2))*cos(X(4));
    x2_dot = -rp*X(3)*sin(X(4))/L;
    x3_dot = q_ddot(1);
    x4_dot = X(5);
    x5_dot = q_ddot(2);
    f_cont = [x1_dot x2_dot x3_dot x4_dot x5_dot]';  % funzione di transizione di stato
    
%     x1_dot_i = rp_i*X(3)*sin(X(2))*cos(X(4));
%     x2_dot_i = -rp_i*X(3)*sin(X(4))/L;
%     x3_dot_i = q_ddot_i(1);
%     x4_dot_i = X(5);
%     x5_dot_i = q_ddot_i(2);
%     f_cont_i = [x1_dot_i x2_dot_i x3_dot_i x4_dot_i x5_dot_i]';  % funzione di transizione di stato
%     
    A = simplify(jacobian(f_cont, X));
    B = simplify(jacobian(f_cont, Tau));
    C = simplify(jacobian(h, X));
    D = simplify(jacobian(h, Tau));
    
%     A_i = simplify(jacobian(f_cont_i, X));
%     B_i = simplify(jacobian(f_cont_i, Tau));

    x_eq = [0 0 1 0 0];
    tau_eq = [0 0];
    A = double(subs(A, [X(1);X(2);X(3);X(4);X(5);Tau(1);Tau(2)], [x_eq'; tau_eq']));                                                                                                                                                                                            
    A_i = [0, (rp_i), 0, 0, 0;
           0, 0, 0,-(rp_i), 0;
           0, 0, 0, 0, 0;
           0, 0, 0, 0, 1;
           0, 0, 0, 0, (20*(rp_i)*(1150*(rp_i)^2 + 9))/(23000*(rp_i)^2 + 180)];
 
    B = double(subs(B, [X(1);X(2);X(3);X(4);X(5);Tau(1);Tau(2)], [x_eq'; tau_eq']));
    B_i = [0 0;
           0 0;
           1200/(23000*(rp_i)^2 + 180) 0;
           0, 0;
           0, 100/23]
    C = double(subs(C, [X(1);X(2);X(3);X(4);X(5);Tau(1);Tau(2)], [x_eq'; tau_eq']))
    D = double(subs(D, [X(1);X(2);X(3);X(4);X(5);Tau(1);Tau(2)], [x_eq'; tau_eq']))
%     A_i = double(subs(A_i, [X(1);X(2);X(3);X(4);X(5);Tau(1);Tau(2)], [x_eq'; tau_eq']));
%     B_i = double(subs(B_i, [X(1);X(2);X(3);X(4);X(5);Tau(1);Tau(2)], [x_eq'; tau_eq']));
%     
    % Controllo proprietà strutturali delle matrici del sistema
    % eig(A);   % autovalori del sistema
    R = rank(ctrb(A,B)) % rango matrice di raggiungibilità
    O = rank(obsv(A,C))  % rango matrice di osservabilità
    [Ao, Bo, Co] = obsvf(A,B,C);    % sistema in forma standard di osservabilità

    j_matrix = struct();
    j_matrix.A = A;
    j_matrix.B = B;
    j_matrix.C = C;
    j_matrix.D = D;
    j_matrix.A_i = A_i;
    j_matrix.B_i = B_i;
end

