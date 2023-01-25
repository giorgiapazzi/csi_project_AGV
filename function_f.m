function x = function_f(x,tau)
    l = 0.4;
    L = 1;
    IPy = 0.15;
    IPz = 0.23;
    Mv = 15;
    mp = 0.5;
    ma = 0.5;
    IAy = 0.12;
    d = 0.6;
    IGz = 0.655;
    IAz = 0.8;
    rp = 0.3;
    ra = 0.3;
    a = Mv+mp+2*ma+2*IAy/(ra^2);
    b = 1/2*(Mv*l^2+mp*L^2+2*ma*d^2+IGz+IPz+2*IAz+2*IAy*(d/ra)^2);

    dt = 0.01;  % tempo di campionamento per la funzione f discretizzata

    % definizione della f (funzione di transizione di stato)
    % Matrice della dinamica
    M=[a*rp^2*cos(x(4))^2+b*(rp/L)^2*sin(x(4))^2+IPy -IPz*(rp/L)*sin(x(4));-IPz*(rp/L)*sin(x(4)) IPz];
    N=[(rp^2*b/L^2-a*rp^2)*cos(x(4))*sin(x(4))*x(5)*x(3);-IPz*rp/L*cos(x(4))*x(5)*x(3)];
    q_ddot = inv(M)*([tau(1); tau(2)]-N);

    % funzione f discretizzata:
    x(1) = x(1) + [rp*x(3)*sin(x(2))*cos(x(4))]*dt;
    x(2) = x(2) + [-rp*x(3)*sin(x(4))/L]*dt;
    x(3) = x(3) + [q_ddot(1)]*dt;
    x(4) = x(4) + [x(5)]*dt;
    x(5) = x(5) + [q_ddot(2)]*dt;
end