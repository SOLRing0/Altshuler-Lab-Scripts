function[invivo]=importinvivostr(filename, startRow)
%%   
delimiter = ',';
startRow = 2;

%% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
formatSpec = '%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Create output variables
Time = dataArray{:,1};
Position = dataArray{:,2};
% get subject ID
    idx=strfind(filename,' ');
Bird_ID=string(filename([(1:idx(1)-1)]));
%% save parameters
invivo.Bird=repelem(Bird_ID,length(Time));
invivo.Time=Time;
invivo.Position=Position;
