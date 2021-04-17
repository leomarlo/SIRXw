function f=PA_SIRXi_w(x,beta,gamma,w,kap,delta)
    % x(1) = rhos
    % x(2) = rhoi
    % x(3) = rhor
    % x(4) = rhosi
    % x(5) = rhoss
    
    f(1)= - beta * x(4);
    f(2)= - (gamma+kap) * x(2) + beta * x(4) ;
    f(3)= gamma * x(2) + delta * (1-x(1)-x(2) -x(3));
    f(4)= -(beta+gamma+w+kap)*x(4) +beta*x(4)*(2*x(5)-x(4))/x(1);
    f(5)= -2*beta*x(4)*x(5)/x(1) + w*(x(1)/(x(1)+x(3)))*x(4);

    f=f';
end