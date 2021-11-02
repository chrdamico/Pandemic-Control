clear all; close all

%% DATA



% Data import from the excel file provided by the "Center for Disease and
% Cost prevention"
% The data was reported by the World Health Organization
[Id_raw, Deaths_raw, Dates] = import_epidemics_data('ebola_data.xlsx');

% The data is noisy and presents some discontinuity points
% thus, a low pass filter will be used to smoothen it
% The epidemic is close to finishing at around 180th time step 
% thus, only the first 180 datapoints will be used. Using more, makes
% the lowpass filter perform poorly since it's applied to a quasi-constant
% variable


% Datapoints to use
datapoints = 180;

% Low pass filter
Id = lowpass(Id_raw(1:datapoints), 0.50);
Deaths = lowpass(Deaths_raw(1:datapoints), 0.50);

% Total subjected population estimate
N = 5500;

% Changing notation in a more "Control" and "identification" setting

% Time 
td = Dates(1:datapoints);
tf = td(length(td));

% Epidemiological data

% Hypothesis: length of the disease = 0
% Thus it is possible to compute S, I and R
Rd =  Id - Deaths;
Sd = N - Id - Rd;
Id = Id; % Useless, still here ease of reading
numdata = length(Id);

% initial condition of the simulation 
x0 = [Sd(1), Id(1), Rd(1)];

% Building the regressor
d = ones(length(Rd), 1);
A = ones(size(Sd));
for i = 2:1:numdata-1    
    ii = (i-1)*3;
    
    % Numerical derivative of the data
    d(ii) = (Sd(i+1) - Sd(i-1))/(td(i+1) - td(i-1));
    d(ii+1) = (Id(i+1) - Id(i-1))/(td(i+1) - td(i-1));
    d(ii+2) = (Rd(i+1) - Rd(i-1))/(td(i+1) - td(i-1));
    
    % Regressor matrix
    A(ii,1) = -Sd(i)*Id(i); 
    A(ii,2) = 0;
    
    A(ii+1,1)= Sd(i)*Id(i); 
    A(ii+1,2) = -Id(i);
    
    A(ii+2,1)= 0.0;
    A(ii+2,2) = Id(i);
end

% Conditioning number
conditioning_number= cond(A)

% Estimation
est = pinv(A)*d

% To check the validity of the obtained solution, a solution of the model
% using it is compared to the actual data

% Solution of the model using estimated parameters
options = odeset('AbsTol', 1e-4, 'RelTol', 1e-4);
[Tx,X] = ode23(@(t, x) stateEq(t,x,est), [0,tf], x0, options);


% Plot of the model evolution
plot(Tx, X(:,1),Tx, X(:,2),Tx, X(:,3));
legend(['S est'; 'I est'; 'R est'])

% Single comaparisons of model and real data

figure()
plot(Tx, X(:,1), Dates(1:datapoints), Sd) 
legend(['S est'; 'S rea'])

figure()
plot(Tx, X(:,2), Dates(1:datapoints), Id) 
legend(['I est'; 'I rea'])

figure()
plot(Tx, X(:,3), Dates(1:datapoints), Rd) 
legend(['R est'; 'R rea'])

%% State equations 
function dx = stateEq(t,x, est)
% Model parameters
a = est(1);
b= est(2);

dx=zeros(3,1);

dx(1)=-a*x(1)*x(2);
dx(2)=a*x(1)*x(2)-b*x(2);
dx(3)=b*x(2);
end


