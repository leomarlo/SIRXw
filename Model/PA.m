function f=PA(x,Model,beta,gamma,w)

    if Model==1
        fss= - 2*beta*x(3)*x(4)/(1-x(1)-x(2)) + w*x(3);
    elseif Model==2
        fss= - 2*beta*x(3)*x(4)/(1-x(1)-x(2)) + w*(1-x(1)-x(2))*x(3)/(1-x(1));
    else
    end

    fi= beta *x(3) - gamma*x(1);
    fr= gamma*x(1);
    fsi= -(beta+gamma+w)*x(3) +beta*x(3)*(2*x(4)-x(3))/(1-x(1)-x(2));


    f=[fi; fr;fsi;fss];
end