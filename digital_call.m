% define function for digital call
function call_price=digital_call(S0, K, r, sig, T)
    dminus=1/(sig*sqrt(T))*(log(S0/K)+(r-1/2*sig^2)*T);
call_price = exp(-r*T)*normcdf(dminus);

