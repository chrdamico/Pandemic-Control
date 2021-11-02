% This script was used to derive all the required symbolic expressions used
% in the "control.m" file

%% Variables
syms l0 l1 l2 l3 dl1 dl2 dl3 N S I R a b c alpha beta gamma ni mi rho k eta u v real

%% Lambda0 and lamda vector
l0=1;
l=[l1, l2, l3]';

%% Lagrangian
L=a*I+b*u+c*v;

%% State vector
x=[S,I,R]';

%% f(x), g1(x) and g2(x)
f=[gamma*N-ni*S-beta*((I*S)/N)+rho*R;
    beta*((I*S)/N)-(ni+mi+alpha)*I;
    alpha*I-(ni+rho)*R];

g1=[-k*S,0,k*S]';

g2=[0,-eta*I,eta*I]';

%% Scalar product
C=dot(l,(f+g1*u+g2*v));

%% Hamiltonian
H = L+C;

%% Equazioni aggiunte
dl1=-diff(H,x(1));
dl2=-diff(H,x(2));
dl3=-diff(H,x(3));

%% lambda derivatives vector
dl=[dl1, dl2, dl3]';

%% H derivative with respect to u
dHu=diff(H,u);

%% H derivative with respect to v
dHv=diff(H,v);

%% [f,g1]
fg1=k*[-gamma*N-rho*(R+S);
    beta*((I*S)/N);
    gamma*N-beta*((I*S)/N)+rho*(R+S)];

%% Phi1 first order derivative
dphi1=dot(l,fg1);

%% [f+g1u+g2v,[f,g1]]
tot=(jacobian(fg1,x))*(f+g1*u+g2*v)-(jacobian(f+g1*u+g2*v,x))*fg1;

%% Phi1 second order derivative
ddphi1=-dot([0,a,0]',fg1)+dot(l,tot);

%% Singular solution for u
%% v = 0
ddphi1_0=subs(ddphi1,v,0);
using_0=solve(ddphi1_0==0,u);
%% v = 1
ddphi1_1=subs(ddphi1,v,1); 
using_1=solve(ddphi1_1==0,u);