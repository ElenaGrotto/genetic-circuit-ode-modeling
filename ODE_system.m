function dydt = ODE_system(t,y,params)

A=y(1);
R=y(2);

dA=params.alfa_a + (params.beta_a*A^params.n)/(params.Ka^params.n+A^params.n)-params.delta*A*R-params.lambda_a*A;
dR=params.alfa_r + (params.beta_r*A^params.p)/(params.Kr^params.p+A^params.p)-params.lambda_r*R;

dydt=[dA; dR];

end