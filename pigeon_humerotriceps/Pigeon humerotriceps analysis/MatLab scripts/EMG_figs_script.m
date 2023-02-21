%% program to analyze all EMG data within a folder

clear
close all

DirName=uigetdir; % select 'Experiments' folder
Q=dir(DirName); 
cd(DirName);
FolderName=uigetdir; % select trial date folder within experiments folder
cd(FolderName);
PathName=uigetdir; % select treatment folder within trial date folder

% change directory
cd(PathName);
A=dir(fullfile(PathName,'*.ddf')); % saves data from all .ddf files into a structure called A 
scrsz = get(0,'ScreenSize');

for a=1:(length(A)) 
    FileName=A(a).name;
f(a)=extract_frequency(FileName,15,15);
    if f==8.6 
    numl=8290; 
elseif f==6.1 
    numl=11306;
elseif f==10.1 
    numl=7262;
    else 
      error('invalid frequency')
    end
[Sample,position,force,AI2,AI3,AI4,AI5,AI6,AI7,AO0,AO1,Stim] = importWLdata(FileName,20,numl);

% identify phase with file name
name = sprintf('%s', string(A(a).name));
fgr=figure();
set(fgr,'Name',name,'NumberTitle','off');
set(fgr,'Unit','normalized','Position',[0,0,1,1]);
plot(1:(length(A)));
hold off

% set x-axis end using frequency
    if f==8.6 
    xend=8000; 
elseif f==6.1 
    xend=11000;
elseif f==10.1 
    xend=7000;
    else 
      error('invalid frequency and x-limit')
    end
   
% plotting parameters 
% subplot EMG trace
subplot(2,1,1);
plot(Sample(1:xend),AI3(1:xend));
xlim([0,xend]);
set(gca,'box','off')
hold on

% subplot stimulus signal
subplot(2,1,2);
plot(Sample(1:xend),AI2(1:xend));
xlabel('Time (10^-^4 s)');
xlim([0,xend]);
set(gca,'box','off')
hold on

 end 

