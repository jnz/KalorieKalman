function [table] = readhealthkitcsv(csvpath)

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 11);

% Specify range and delimiter
opts.DataLines = [3, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["type", "sourceName", "sourceVersion", "productType", "device", "startDate", "endDate", "unit", "value", "HKExt1", "HKExt2"];
opts.VariableTypes = ["string", "string", "string", "string", "string", "datetime", "datetime", "string", "double", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["type", "sourceName", "sourceVersion", "productType", "device", "unit", "HKExt1", "HKExt2"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["type", "sourceName", "sourceVersion", "productType", "device", "unit", "HKExt1", "HKExt2"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "startDate", "InputFormat", "yyyy-MM-dd HH:mm:ss Z", "TimeZone", 'Europe/Berlin');
opts = setvaropts(opts, "endDate", "InputFormat", "yyyy-MM-dd HH:mm:ss Z", "TimeZone", 'Europe/Berlin');
opts = setvaropts(opts, "value", "ThousandsSeparator", ",");

% Import the data
table = readtable(csvpath,opts);
% "yyyy-MM-dd HH:mm:ss Z"


end
