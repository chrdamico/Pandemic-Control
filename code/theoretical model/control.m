    global x0 ptf tf step gamma ni beta mi alpha rho k eta a b c;

%{
This script controls a nonlinear epidemiological SIR model with varying
population described in the paper "Optimal Control for a SIR Epidemiological 
Model with Time-varying Populations". 2 control variables are available.
The numerical solution to the problem is obtained using the
forwards-backwards method. All the relevant parameters are the same as in
the paper.
%}



% State initial condition
x0=[1000;10;0];
% Costate final condition
ptf=[0;0;0];
% Time of simulation
tf=100;

% Numerical step in control computation
step=0.25;

% Model parameters
gamma=0.00683;
ni=0.00188;
beta=0.2426;
mi=0.005;
alpha=0.00002;
rho=0.007;
k=0.3;
eta=0.1;
a=5;
b=50;
c=300;

% Time initial condition
Tu=linspace(0,tf);
t=linspace(0,tf);

% Control initial guess
u=randn(size(Tu));
v=randn(size(Tu));

%% Computation loop
for i=1:200
    
    % Forward solution
   options = odeset('AbsTol', 1e-4, 'RelTol', 1e-4);
   [Tx,X] = ode23(@(t, x) stateEq(t,x,u,v,Tu), [0,tf], x0, options);
   
   x1 = X(:,1);
   x2 = X(:,2);
   x3 = X(:,3);
   
   % Backwards solution
   options = odeset('AbsTol', 1e-2, 'RelTol', 1e-2);
   [Tp,P] = ode23(@(t, p) costateEq(t, p, u, v, Tu, x1, x2, x3, Tx), [tf,0], ptf, options);
   
   p1=P(:,1);
   p1=interp1(Tp,p1,Tx);
   
   p2=P(:,2);
   p2=interp1(Tp,p2,Tx);
   
   p3=P(:,3);
   p3=interp1(Tp,p3,Tx);
   
   % Derivative of the Hamiltonina with respect to v used in control update and as
   % an exit condition
   dHu=pHu(x1,p1,p3,Tx,u,Tu);
   Hu_norm = dHu'*dHu;
   
   dHv=pHv(x2,p1,p2,Tx,v,Tu);
   Hv_norm = dHv'*dHv;
   
   if (Hu_norm+Hv_norm) < 1.0e-6
       break
   else
       % Control update
       u_old = u;
       v_old = v;
       u = AdjControlu(dHu, Tx, u_old, Tu, step);
       v = AdjControlv(dHv, Tx, v_old, Tu, step);
   end
    
end

% State plot
figure; 
plot(Tx, X(:,1),Tx, X(:,2),Tx, X(:,3));
xlabel('t'); 
legend(['S'; 'I'; 'R'])

% Control plot
figure; 
plot(Tu, u, Tu, v);
xlabel('t'); 
legend('u','v');


%% State equations
function dx = stateEq(t,x,u,v,Tu)
global x0 ptf tf step gamma ni beta mi alpha rho k eta a b c;
dx=zeros(3,1);
u=interp1(Tu,u,t);
v=interp1(Tu,v,t);
N = x(1)+x(2)+x(3);
dx(1)=gamma*N-ni*(x(1))-beta*((x(2)*x(1))/N)+rho*x(3)-k*x(1)*u;
dx(2)=beta*((x(2)*x(1))/N)-(ni+mi+alpha)*x(2)-eta*x(2)*v;
dx(3)=-ni*x(3)-rho*x(3)+k*x(1)*u+alpha*x(2)+eta*x(2)*v;
end

%% Costate equations
function dp = costateEq(t, p, u, v, Tu, x1, x2, x3, xt)
global x0 ptf tf step gamma ni beta mi alpha rho k eta a b c;
dp = zeros(3,1);
x1=interp1(xt, x1, t);
x2=interp1(xt, x2, t);
x3=interp1(xt, x3, t);

u = interp1(Tu, u, t);
v = interp1(Tu, v, t);

N = x1+x2+x3;

dp(3) = p(3)*(ni + k*u + (x2*beta)/N) - k*p(1)*u - (x2*beta*p(2))/N;
dp(2) = p(2)*(alpha + mi + ni + eta*v - (x1*beta)/N) - a - p(1)*(alpha + eta*v) + (x1*beta*p(3))/N;
dp(1) = p(1)*(ni + rho) - p(3)*rho;
end

%% dHu
function dHu = pHu(x1,p1,p3,Tx,u,Tu)
global b k;
u = interp1(Tu,u,Tx);
dHu = (2*b).*u + k*x1.*(p1-p3);
end

%% dHv
function dHv = pHv(x2,p1,p2,Tx,v,Tu)
global c eta;
v = interp1(Tu,v,Tx);
dHv = (2*c).*v+eta*x2.*(p1-p2);
end

%% Adjust control u
function u_new = AdjControlu(apHu, tx, u, Tu, step)
global b;
apHu = interp1(tx, apHu, Tu);
u_new = min(1, max(0,u-step*(apHu/(2*b))));
end

%% Adjust control v
function v_new = AdjControlv(bpHv, tx, v, Tv, step)
global c;
bpHv = interp1(tx, bpHv, Tv);
v_new = min(1, max(0,v-step*(bpHv/(2*c))));
end