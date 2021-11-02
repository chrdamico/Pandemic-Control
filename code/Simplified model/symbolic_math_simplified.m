%{
This script was used to generate symbolic quantities used in the
"control_simplified_model.m" file
%}

%% Variables
syms l0 l1 l2 l3 dl1 dl2 dl3 N S I R a_cost b_cost c_cost alpha beta gamma ni mi rho k eta u v real
syms a b d_cost real

%% Lambda0 and lamda vector
l0=1;
l=[l1, l2, l3]';

%% Lagrangian
L=a_cost*I^2+b_cost*u^2+c_cost*v^2;

%% State vector
x=[S,I,R]';

%% f(x), g1(x) and g2(x)
f=[-a*S*I;
    a*S*I-b*I;
    b*I];


g2=[0,-eta*I,eta*I]';

%% Scalar product
C=dot(l,(f+g2*v));

%% Hamiltonian
H = L+C;

%% Equazioni aggiunte
dl1=-diff(H,x(1));
dl2=-diff(H,x(2));
dl3=-diff(H,x(3));

%% lambda derivatives vector
dl=[dl1, dl2, dl3]';


%% H derivative with respect to v
dHv=diff(H,v);
