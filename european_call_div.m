% define function for BS formula with dividend rate q
function call_price=european_call_div(S0, K, r, sig, T, q)
    dplus=1/(sig*sqrt(T))*(log(S0/K)+(r-q+1/2*sig^2)*T);
    dminus=1/(sig*sqrt(T))*(log(S0/K)+(r-q-1/2*sig^2)*T);

call_price = exp(-q*T)*S0*normcdf(dplus)-normcdf(dminus)*K*exp(-r*T);

