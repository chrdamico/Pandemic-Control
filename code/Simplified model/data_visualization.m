clear all




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

Id = lowpass(Id_raw(1:180), 0.5);
Deaths = lowpass(Deaths_raw(1:180), 0.5);  

% Plot of the data
figure()
plot(Dates(1:180), Id, Dates(1:180), Deaths)
legend('Infected','Deaths')

% Plot of raw data
figure
plot(Dates(1:180), Id_raw(1:180), Dates(1:180), Deaths_raw(1:180))
legend('Raw Infected', 'Raw Deaths')


