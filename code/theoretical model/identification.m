clf(figure)
clear
%{
This script tries to estimate the model of the system in the paper "Optimal Control for a 
SIR Epidemiological Modelwith Time-varying Populations" using the
previously solved optimal control data from the file "control.m". This
identification does not work properly at all, and even little added noise makes
the results terrible. This can be seen from the really high conditioning
number of the regressor. Becouse of this, a simplified model was used to
model real world ebola data as seen in the "_simplified" files.
%}

%% Simulating the system

% Initial conditions
tf = 50;
x0 = [1000, 30, 0];

% Time initial condition
Tu=linspace(0,tf);
t=linspace(0,tf);

% Hardcoded optimal input
Ud = [1,1,1,1,1,1,1,1,0.968975258664535,0.835972250677274,0.766964837833510,0.717829513121885,0.688082664419035,0.661921969580057,0.638695971959779,0.622984631993726,0.611370044533465,0.602617298313306,0.595554942386413,0.589400032755461,0.584507810452772,0.579915373612029,0.575802534788047,0.571680931330957,0.567761040096223,0.563800693242479,0.559852766702756,0.556066243477167,0.551740462726636,0.547742874924934,0.544159248011150,0.541052413638480,0.538249991170953,0.535290379136444,0.531848252686945,0.528467226414745,0.525107352792290,0.521788119300718,0.518008072047296,0.513998567107418,0.510391203251002,0.507370712135997,0.504306577694856,0.501209521409481,0.497825668681742,0.494251258560809,0.490686662094824,0.486914685684356,0.483355325085999,0.479803794143267,0.476255615169393,0.472523188240093,0.468847549691350,0.465209020881004,0.461508248694335,0.457814074782951,0.453887125904862,0.449946601860613,0.445903891548508,0.441747631591038,0.437304040602459,0.432949625464375,0.428620926569270,0.424005381125835,0.419236791988639,0.414509460765611,0.409501063193467,0.404447973462088,0.399278693734699,0.394493150684554,0.388104776125962,0.380020449222713,0.372599507872502,0.364868590219668,0.356969627687001,0.349020974346037,0.340737128507770,0.332504021351399,0.323423302639476,0.314955397219170,0.303499829360809,0.289776939766344,0.276388352482314,0.261936373376595,0.247604087776336,0.231914228615667,0.216151302555999,0.198846815751408,0.182220035983445,0.166062323348898,0.147826417323592,0.128746755353282,0.107713295625979,0.0850115244641498,0.0593391382567978,0.0459667537374715,0.0335465024411029,0.0185669742568344,0.00279537155668542,0];
Vd = [0.574915527109439,0.579912235572687,0.565646455896653,0.557697951387214,0.541008684217111,0.535348109386659,0.524077059274393,0.507722236495374,0.494231837310565,0.484731129123938,0.474088097068005,0.463688800765647,0.453288614119364,0.443483992735201,0.434118459959794,0.425145377477690,0.416582116303506,0.408460467291915,0.400727867346233,0.393333443987490,0.386426244012125,0.378126516991714,0.370291946271536,0.362911529223711,0.355973751902006,0.349562036366638,0.343253564725893,0.337224509249406,0.331859537875024,0.326668108653262,0.321242141687208,0.314833380966317,0.308556330715345,0.302557259818921,0.297245515268557,0.291949122767707,0.286729334977581,0.281724684409151,0.277366464437644,0.273151169205749,0.268817221440492,0.263818694411572,0.258386908152788,0.252909607053646,0.247996180269474,0.243521659364794,0.239234611532719,0.235085673105781,0.231043766434424,0.226691561228944,0.222386833952903,0.218020006845304,0.213628707100921,0.209368700738265,0.205107344359431,0.200782848297718,0.196614564703356,0.192780350412559,0.188962504990054,0.185309661259223,0.181412420160449,0.177100007411409,0.172975904969085,0.168511240050634,0.164331233236855,0.160327434171514,0.156352238062284,0.152533528164517,0.148754111934603,0.145245304262123,0.141275076850862,0.136773321231741,0.132572767445885,0.128347004702500,0.124159593837745,0.120060792283663,0.115975865589152,0.111921202674268,0.107840405915601,0.103849829301427,0.0995994561783935,0.0951413529032418,0.0907215774399359,0.0862363685582093,0.0817297069787457,0.0771630602616126,0.0725293005162590,0.0677961020057889,0.0630570584356491,0.0582503606757993,0.0532998239953349,0.0481709990741101,0.0428821438793109,0.0372927702512229,0.0314244616583221,0.0258339364563065,0.0199605146883699,0.0136413937261810,0.00682986996910930,0];

% Solution of state equations
options = odeset('AbsTol', 1e-4, 'RelTol', 1e-4);
[Tx,X] = ode23(@(t, x) stateEq(t,x,Ud,Vd,Tu), [0,tf], x0, options);

% Renaming variables
Sd = X(:,1);
Id = X(:,2);
Rd = X(:,3);
td = Tx';

%% Parameter estimates
% Interpolation to make the size of Ud and Vd the same as the size of the
% optimal inputs
Ud = interp1(Tu,Ud,Tx);
Vd = interp1(Tu,Vd,Tx);

% Number of data-points available to get an estimate
numdata = 35;

% Building the regressor and data vectors 
% d = A est where
% d is the data vector
% A is the regressor
% est is the vector containing the estimated parameters
d = zeros(numdata*3, 1);
A = zeros(numdata*3, 8);

for i = 2:1:numdata-1
    % Every three elements the same terms are present
    ii = (i-1)*3;
    % Data vector
    d(ii) = (Sd(i+1) - Sd(i-1))/(td(i+1) - td(i-1));
    d(ii+1) = (Id(i+1) - Id(i-1))/(td(i+1) - td(i-1));
    d(ii+2) = (Rd(i+1) - Rd(i-1))/(td(i+1) - td(i-1));


    % Regressor (parametrized in order wrt: gamma, ni, beta, rho, k, eta,
    % alpha, mi)
    A(ii,1) = -Sd(i)+Rd(i)+Id(i);
    A(ii,2) = -Sd(i);
    A(ii,3) = -(Sd(i)*Id(i))/(Sd(i)+Id(i)+Rd(i));
    A(ii,4) = Rd(i);
    A(ii,5) = -Sd(i)*Ud(i);
    A(ii,6) = 0;
    A(ii,7) = 0;
    A(ii,8) = 0;

    A(ii+1,1)= 0;
    A(ii+1,2)= -Id(i);
    A(ii+1,3)= (Sd(i)*Id(i))/(Sd(i)+Id(i)+Rd(i));
    A(ii+1,4)= 0;
    A(ii+1,5)= 0;
    A(ii+1,6)= -Id(i)*Vd(i);
    A(ii+1,7)= -Id(i);
    A(ii+1,8)= -Id(i);

    A(ii+2,1)= 0.0;
    A(ii+2,2)= -Rd(i);
    A(ii+2,3)= 0;
    A(ii+2,4)= -Rd(i);
    A(ii+2,5)= Sd(i)*Ud(i);
    A(ii+2,6)= Id(i)*Vd(i);
    A(ii+2,7)= Id(i);
    A(ii+2,7)= 0;
end

% Condition number of the regressor
cond(A)

% Computation of the estimate using the pseudoinverse of the regressor
est = pinv(A)*d;
disp(est)

%% Evolution of the system with the new parameters
% Interpolation to make the size of Ud and Vd the same as the size of the
% optimal inputs
Ud = interp1(Tx,Ud,Tu);
Vd = interp1(Tx,Vd,Tu);

% Solution of the system with new parameters
options = odeset('AbsTol', 1e-4, 'RelTol', 1e-4);
[Tx_est,X_est] = ode23(@(t, x) stateEq_est(t,x,Ud,Vd,Tu,est), [0,tf], x0, options);


%% plot
X_est = interp1(Tx_est, X_est, Tx);

plot(Tx, X_est(:,1),Tx, X_est(:,2),Tx, X_est(:,3), Tx, X(:,1),Tx, X(:,2),Tx, X(:,3));
legend(['S est'; 'I est'; 'R est'; 'S rea'; 'I rea'; 'R rea'])


%% State equations computed with the estimated parameters
function dx = stateEq_est(t,x,u,v,Tu, est)
% estimated parameters
% parametrized in order wrt: gamma, ni, beta, rho, k, eta, alpha, mi)
gamma = est(1);
ni = est(2);
beta = est(3);
rho = est(4);
k = est(5);
eta = est(6);
alpha = est(7);
mi = est(8);

dx=zeros(3,1);
u=interp1(Tu,u,t);
v=interp1(Tu,v,t);
N = x(1)+x(2)+x(3);
dx(1)=gamma*N-ni*(x(1))-beta*((x(2)*x(1))/N)+rho*x(3)-k*x(1)*u;
dx(2)=beta*((x(2)*x(1))/N)-(ni+mi+alpha)*x(2)-eta*x(2)*v;
dx(3)=-ni*x(3)-rho*x(3)+k*x(1)*u+alpha*x(2)+eta*x(2)*v;
end



%% State equations with previous parameters
function dx = stateEq(t,x,u,v,Tu)
% Model parameters
gamma=0.00683;
ni=0.00188;
beta=0.2426;
mi=0.005;
alpha=0.00002;
rho=0.007;
k=0.3;
eta = 0.1;

dx=zeros(3,1);
u=interp1(Tu,u,t);
v=interp1(Tu,v,t);
N = x(1)+x(2)+x(3);
dx(1)=gamma*N-ni*(x(1))-beta*((x(2)*x(1))/N)+rho*x(3)-k*x(1)*u;
dx(2)=beta*((x(2)*x(1))/N)-(ni+mi+alpha)*x(2)-eta*x(2)*v;
dx(3)=-ni*x(3)-rho*x(3)+k*x(1)*u+alpha*x(2)+eta*x(2)*v;
end



