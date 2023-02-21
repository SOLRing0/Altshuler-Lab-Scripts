%% program to analyze all data within a folder

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

% name figures with filenames
name = sprintf('%s',string(A(a).name));
fgr=figure();
set(fgr,'Name',name,'NumberTitle','off');
set(fgr,'Unit','normalized','Position',[0,0,1,1]);
plot(1:(length(A)));
hold off

GR=4;  % Choose gear ratio

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

% set plotting parameters
[AX, H1, H2]=plotyy(Sample(1:xend),((position(1:xend))/GR),Sample(1:xend),((force(1:xend))*GR));
set(H1,'linewidth',2);
set(H2,'linewidth',2);
xlabel('Time (10^-^4 s)');
xlim(AX(1), [0,xend]);
xlim(AX(2),[0,xend]);
set(get(AX(1),'Ylabel'),'String','\Delta Length (mm)'); 
ylim(AX(1),[-3,3]);
set(AX(1),'linewidth',1);
set(AX(2),'linewidth',1);
yticks(AX(1),[-3,0,3]);
ylim(AX(2), [0,40]);
yticks(AX(2), [0:5:40]);
set(get(AX(2),'Ylabel'),'String','Force (N)'); 
set(gca,'box','off')

 end 

