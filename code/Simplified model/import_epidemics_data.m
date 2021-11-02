function [Id,Deaths,Dates] = import_epidemics_data(filename)

table = readtable(filename);

Id_table = table(1:end,2);
Deaths_table = table(1:end,3);
Dates_table = table(1:end, 10);


Id = flipud(table2array(Id_table));
Deaths = flipud(table2array(Deaths_table));
Dates = flipud(table2array(Dates_table));
end