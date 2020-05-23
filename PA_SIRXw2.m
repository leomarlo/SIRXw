function f=PA_SIRXw2(x,beta,gamma,w,kap0,kap,delta)
    % x(1) = rhos
    % x(2) = rhor
    % x(3) = rhoxs
    % x(4) = rhoxi
    
    % x(5) = rhosi
    % x(6) = rhoss
    % x(7) = rhosxs
    % x(8) = rhoixs
    
    f(1)= - beta *x(5) - kap0*x(1) + delta*x(3);    
    f(2)= beta * (1-x(1)-x(2) -x(3)-x(4)) + delta* x(4) ; % There is a gamma missing!!
    f(3)= kap0* x(1) - delta *x(3);
    f(4)= (kap0+kap)* (1-x(1)-x(2) -x(3)-x(4))  - delta *x(4);
  
    f(5)= -(beta+gamma+w+2*kap0+kap)*x(5) +beta*x(5)*(2*x(6)-x(5))/x(1) + delat * x(8);
    f(6)= -2*beta*x(5)*x(6)/x(1) + w*(x(1)/(x(1)+x(2)))*x(5) - 2*kap0*x(6) + delta *x(7);
    f(7)= -delta*x(7) + 2*kap0*x(6) - beta * x(7)*(x(5)/x(1));
    f(8)= - delta*x(8) - gamma*x(8) - (kap0+kap)*x(8) + kap*x(5);

    f=f';
end