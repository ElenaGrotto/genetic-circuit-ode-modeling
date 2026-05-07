clc, clear, close all

%Parameters
params.alfa_a=0.00875;
params.alfa_r=0.025;
params.beta_a=7.5;
params.beta_r=2.5;
params.Kr=0.14*2.5*10^4;
params.delta=4.10^(-8);
params.lambda_a=10^(-4);
params.lambda_r=10^(-4);
params.Ka=0.2*2.5*10^4;
params.n=2;
params.p=5;

%Initial conditions
A0=0;
R0=0;
y0=[A0 R0];
tspan=[0 10^5];

%Exact solution
options_exact=odeset('RelTol',1e-10,'AbsTol',1e-12);
[t_exact,y_exact]=ode45(@(t,y) ODE_system(t,y,params),tspan,y0, options_exact);

%Solvers to test
solvers={@ode45, @ode15s, @ode23, @ode23s, @ode23tb};
solver_names={'ode45','ode15s','ode23','ode23s','ode23tb'};

times=zeros(size(solvers));
steps=zeros(size(solvers));
rel_err_A=zeros(size(solvers));
rel_err_R=zeros(size(solvers));
results=struct;

%Loop over solvers
for i=1:length(solvers)
    tic
    [t,y]=solvers{i}(@(t,y) ODE_system(t,y,params), tspan, y0);
    times(i)=toc;
    steps(i)=length(t);
    
    %Relative error on A and R
    A_ref=interp1(t_exact,y_exact(:,1),t);
    R_ref=interp1(t_exact,y_exact(:,2),t);
    rel_err_A(i)=norm(y(:,1)-A_ref)/norm(A_ref);
    rel_err_R(i)=norm(y(:,2)-R_ref)/norm(R_ref);

    results(i).t=t;
    results(i).y=y;
end

%Results table
Table=table(solver_names',times',steps',rel_err_A',rel_err_R', ...
    'VariableNames',{'Method','CPU time','Steps','RelError_A','RelError_R'});
disp('                === Solver comparison ===')
disp(Table);

%Plots
figure(1)
subplot(211)
hold on

for i=1:length(solvers)
    plot(results(i).t,results(i).y(:,1),'DisplayName',[solver_names{i} 'A(t)'])
end

plot(t_exact,y_exact(:,1),'k--','LineWidth',1.5,'DisplayName','Reference');
xlabel('Time'); ylabel('A(t)');
legend show;
title('Comparison of A(t) solutions');

subplot(212)
hold on
for i=1:length(solvers)
    plot(results(i).t,results(i).y(:,2),'DisplayName',[solver_names{i} 'R(t)'])
end
plot(t_exact,y_exact(:,2),'k--','LineWidth',1.5,'DisplayName','Reference');
xlabel('Time'); ylabel('R(t)');
legend show;
title('Comparison of R(t) solutions');

%Convergence test
tolerances = [1e-3 1e-4 1e-6 1e-8];
methods = {@ode45, @ode15s, @ode23s, @ode23tb};
names = {'ode45','ode15s','ode23s','ode23tb'};

rows = [];  

for m=1:length(methods)
    solver=methods{m};
    for tol=tolerances
        opts=odeset('RelTol',tol,'AbsTol',tol*1e-2);
        tic;
        [t,y]=solver(@(t,y) ODE_system(t,y,params), tspan, y0, opts);
        time=toc;
       
        A_ref=interp1(t_exact,y_exact(:,1),t);
        R_ref=interp1(t_exact,y_exact(:,2),t);
        errA=norm(y(:,1)-A_ref)/norm(A_ref);
        errR=norm(y(:,2)-R_ref)/norm(R_ref);
       
        rows=[rows; {names{m}, tol, length(t), time, errA, errR}];
    end
end

ConvTable=cell2table(rows, ...
    'VariableNames', {'Method','RelTol','Steps','CPU_time','RelError_A','RelEerror_R'});

disp('                        === Convergence study ===')
disp(ConvTable)
