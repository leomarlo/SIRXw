%% check implicit formula

function f = implicitrinf(beta,rates,rho0,rhosi0,rhoss0,rinf)
    f = rho0 + (beta/(beta+rates)) * ( rhosi0 + (2 * rhoss0 -rhosi0 ) * log( (1-rho0)/(1-rinf)  ) );
end
