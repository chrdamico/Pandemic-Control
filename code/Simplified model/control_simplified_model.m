clear all
close all 

%{
In this script, we control the model estimated from real world data in the
file "identification_simplified.m". The control is computed numerically
using the forwards-backwards method. Costs are the same as in the paper 
"Optimal Control for a SIR Epidemiological Modelwith Time-varying Populations".
Initial conditions and length of the simulation are the ones from
the solution of ODE23 in the identification file. Effectiveness of the cure
has been substantially decreased to make the control more reasonable
(otherwise, the epidemic was completely stopped in a few days using the
parameters in the paper).
%}



% Global parameters used inside the functions
global x0 ptf tf step a b a_cost b_cost c_cost eta 

% State initial condition
x0=[5392.9;81.7;25.4];

% Costate final condition
ptf=[0;0;0];

% Time of simulation
tf=521;

% Numerical step in control computation
step=0.05;

% Model parameters
a = 2.91172870512048e-06;
b = 0.000837104202224668;
% Costs
a_cost = 5;
c_cost = 5000000;
% Control effect
eta = 0.005;

% Time initial condition
Tu=linspace(0,tf);
t=linspace(0,tf);

% Control initial guess
v=randn(size(Tu));

%% Computation loop
for i=1:100
    i
    % Forward solution
   options = odeset('AbsTol', 1e-4, 'RelTol', 1e-4);
   [Tx,X] = ode23(@(t, x) stateEq(t,x,v,Tu), [0,tf], x0, options);
   
   x1 = X(:,1);
   x2 = X(:,2);
   x3 = X(:,3);
   
   % Backwards solution
   options = odeset('AbsTol', 1e-2, 'RelTol', 1e-2);
   [Tp,P] = ode23(@(t, p) costateEq(t, p, v, Tu, x1, x2, x3, Tx), [tf,0], ptf, options);
   
   p1=P(:,1);
   p1=interp1(Tp,p1,Tx);
   
   p2=P(:,2);
   p2=interp1(Tp,p2,Tx);
   
   p3=P(:,3);
   p3=interp1(Tp,p3,Tx);
   
   % Derivative of the Hamiltonina with respect to v used in control update and as
   % an exit condition
   dHv=pHv(x2,p1,p2,Tx,v,Tu);
   Hv_norm = dHv'*dHv;
   
   if (Hv_norm) < 1.0e-6
       break
   else
       % Control update
       v_old = v;
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
plot(Tu, v);
xlabel('t'); 
legend('v');


%% State equations
function dx = stateEq(t,x,v,Tu)
global x0 ptf tf step k eta a b;

dx=zeros(3,1);

v=interp1(Tu,v,t);

dx(1)= a*-x(1)*x(2);
dx(2)= a*x(1)*x(2) - b*x(2) - eta*x(2)*v;
dx(3)= b*x(2) + eta*x(2)*v;
end

%% Costate equations
function dp = costateEq(t, p, v, Tu, x1, x2, x3, xt)
global x0 ptf tf step k eta a b a_cost;
dp = zeros(3,1);
x1=interp1(xt, x1, t);
x2=interp1(xt, x2, t);
x3=interp1(xt, x3, t);

v = interp1(Tu, v, t);


dp(3) = 0;
dp(2) = p(2)*(b - a*x1 + eta*v) - p(1)*(b+eta*v) - 2*a_cost*x2 + a*x1*p(3);
dp(1) = - a*x2*p(2);

end



%% dHv
function dHv = pHv(x2,p1,p2,Tx,v,Tu)
global c_cost eta;
v = interp1(Tu,v,Tx);
dHv = (2*c_cost).*v+eta*x2.*(p1-p2);
end


%% Adjust control v
function v_new = AdjControlv(bpHv, tx, v, Tv, step)
global c_cost;
bpHv = interp1(tx, bpHv, Tv);
v_new = min(1, max(0,v-step*(bpHv/(2*c_cost))));
end