function j_matrix = get_linearization_lqg()
    clc
    global rp l L IPy IPz Mv mp ma IAy ra d IGz IAz a b 
   
    %% Vettore di stato X = [x y theta phi_dot psi psi_dot]
    %Tau=[tau_phi tau_psi]  ingressi
    Tau=sym('tau',[2 1],'real');
    
    %X vettore di stato
    X = sym('x',[6 1],'real');
    
    % definizione della h (modello di osservazione)
    % vettore di variabili misurate [X(3),X(4),X(5),X(6)]
    h = [X(1),X(2),X(4),X(5)]';
    
    % definizione della f (funzione di transizione di stato)
    % Matrice della dinamica
    M=[a*rp^2*cos(X(5))^2+b*(rp/L)^2*sin(X(5))^2+IPy -IPz*(rp/L)*sin(X(5));-IPz*(rp/L)*sin(X(5)) IPz];
    N=[(rp^2*b/L^2-a*rp^2)*cos(X(5))*sin(X(5))*X(6)*X(4);-IPz*rp/L*cos(X(5))*X(6)*X(4)];
    q_ddot = inv(M)*([Tau(1); Tau(2)]-N);
        
    %X_dot
    x1_dot = rp*X(4)*cos(X(3))*cos(X(5));
    x2_dot = rp*X(4)*sin(X(3))*cos(X(5));
    x3_dot = -rp*X(4)*sin(X(5))/L;
    x4_dot = q_ddot(1);
    x5_dot = X(6);
    x6_dot = q_ddot(2);
    f_cont = [x1_dot x2_dot x3_dot x4_dot x5_dot x6_dot]';  % funzione di transizione di stato
    A = simplify(jacobian(f_cont, X));
    B = simplify(jacobian(f_cont, Tau));
    C = simplify(jacobian(h, X));
    D = simplify(jacobian(h, Tau));
    
    x_eq = [0 pi 0 0];
    tau_eq = [0 0];
    A = double(subs(A, [X(3);X(4);X(5);X(6);Tau(1);Tau(2)], [x_eq'; tau_eq']));
    B = double(subs(B, [X(3);X(4);X(5);X(6);Tau(1);Tau(2)], [x_eq'; tau_eq']));
    C = double(subs(C, [X(3);X(4);X(5);X(6);Tau(1);Tau(2)], [x_eq'; tau_eq']));
    D = double(subs(D, [X(3);X(4);X(5);X(6);Tau(1);Tau(2)], [x_eq'; tau_eq']));

    % Controllo proprietà strutturali delle matrici del sistema
    % eig(A);   % autovalori del sistema
    % R = rank(ctrb(A,B));  % rango matrice di raggiungibilità
    % O = rank(obsv(A,C));  % rango matrice di osservabilità
    [Ao, Bo, Co] = obsvf(A,B,C);    % sistema in forma standard di osservabilità

    j_matrix = struct();
    j_matrix.A = A;
    j_matrix.B = B;
    j_matrix.C = C;
    j_matrix.D = D;
end

