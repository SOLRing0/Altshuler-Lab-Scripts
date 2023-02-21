function[insitu] = importinsitustr(filename, startRow)
%%   
delimiter = ',';
startRow = 2;

%% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f) 
formatSpec = '%f%f%f%[^\n\r]';  
%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Create output variable
Rep_num = dataArray{:,1};
Time = dataArray{:,2};
Position = dataArray{:,3};
% get strain frequency
    stix=strfind(filename,'-');
stfreq=str2double(filename(stix(1)+1:stix(2)-1)); %'+1' excludes the preceding '-' and '-1' excludes the anteceding '-'

%% save parameters
insitu.Rep_num=Rep_num;
insitu.Time=Time;
insitu.Position=Position;
insitu.freq=repelem(stfreq,length(Time));