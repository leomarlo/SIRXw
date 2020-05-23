function f=PA_SIRXw(x,beta,gamma,w,kap0,kap)
    % x(1) = rhos
    % x(2) = rhoi
    % x(3) = rhor
    % x(4) = rhosi
    % x(5) = rhoss
    
    f(1)= - beta *x(4) - kap0*x(1);
    f(2)= beta * x(4) - (gamma+kap0+kap)*x(2)  ;
    f(3)= gamma*x(2);
    f(4)= -(beta+gamma+w+2*kap0+kap)*x(4) +beta*x(4)*(2*x(5)-x(4))/x(1);
    f(5)= -2*beta*x(4)*x(5)/x(1) + w*(x(1)/(x(1)+x(3)))*x(4) - 2*kap0*x(5);


    f=f';
end