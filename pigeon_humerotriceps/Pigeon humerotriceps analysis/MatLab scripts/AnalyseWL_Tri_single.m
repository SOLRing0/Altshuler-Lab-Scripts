% code for running a single trial
clear


PathName=uigetdir; % select treatment folder within an experiment folder
cd(PathName);
% Calls file 
FileName=uigetfile('*.ddf');

%Program to analyze workloop data, calculating the work per cycle

% get trial frequency
f=extract_frequency(FileName,15,15);

% set line number for endRow using frequency setting
if f==8.6 
    numl=8290; 
elseif f==6.1 
    numl=11306;
elseif f==10.1 
    numl=7262;
    else 
      error('invalid frequency and endRow value')
end

[Sample,position,force,AI2,AI3,AI4,AI5,AI6,AI7,AO0,AO1,Stim] = importWLdata(FileName,20,numl);  

GR=4; % Choose gear ratio

%plot raw data
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/1.3 scrsz(4)/3])
subplot(1,4,[1 2])

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

[AX, H1, H2]=plotyy(Sample(1:xend),(position(1:xend))/GR,Sample(1:xend),(force(1:xend))*GR);
xlabel('Time (10^-^4 s)');
set(get(AX(1),'Ylabel'),'String','\Delta Length (mm)'); 
set(get(AX(2),'Ylabel'),'String','Force (N)'); 
set(AX(1),'YLim',[-2,2]);
set(AX(2),'YLim',[-0.5,40]); 

%% crop and stack data
% set startdex search parameters -- end point of search range
if f==8.6 
    endpt=1200; 
elseif f==6.1 
    endpt=1700;
elseif f==10.1 
    endpt=1200;
    else 
      error('invalid frequency and startdex parameter value')
end

startdex=find(position(1:endpt)==max(position(1:endpt)));
enddex=find(position(end-500:end)==max(position(end-500:end)))+length(position)-500;

% set EndLoop search parameters
% search range start
if f==8.6 
    ELs=1000; 
elseif f==6.1 
    ELs=1200;
elseif f==10.1 
    ELs=900;
    else 
      error('invalid frequency and EndLoop range start')
end
% search range end
if f==8.6 
    ELe=1500; 
elseif f==6.1 
    ELe=1700;
elseif f==10.1 
    ELe=1400;
    else 
      error('invalid frequency and EndLoop range end')
end

begLoop=startdex(1);
EndLoop=find(position(begLoop+ELs:begLoop+ELe)==max(position(begLoop+ELs:begLoop+ELe)))+begLoop+ELs;
loopLength=EndLoop-begLoop+1;

for i=1:5
    
    endLoop=loopLength+begLoop;

    Lengths(1:loopLength+1,i)=position(begLoop:endLoop)/GR;
    Forces(1:loopLength+1,i)=force(begLoop:endLoop)*GR;
    begLoop=endLoop;
end

subplot(1,4,3); hold on
plot(Lengths,Forces);
xlabel('\Delta Length (mm)'); ylabel('Force (N)');
xlim([-3,3]);

%% Calculate work

DS=1:round(length(Lengths)/2);
US=round(length(Lengths)/2):length(Lengths)-2;

for i=1:size(Lengths,2)
    uppercurve(1,i)=trapz(flipud(Lengths(DS,i)),flipud(Forces(DS,i)));
    lowercurve(1,i)=trapz(Lengths(US,i),Forces(US,i));
    work(1,i)=uppercurve(1,i)-lowercurve(1,i);
end

% Net work per cycle plot
subplot(1,4,4);
plot(1:5,work,'bo');
xlabel('cycle number'); ylabel('Work (mJ)');
xlim([0.5,5.5]);
ylim([-2.5,3.5]);

Work5=mean(work);
Work3=mean(work(3:5));
Work2=mean(work(4:5));
Work1=mean(work(end));

WorkTable=nan(1,4);
WorkTable(:,1)=Work1';
WorkTable(:,2)=Work2';
WorkTable(:,3)=Work3';
WorkTable(:,4)=Work5';

PowerTable=WorkTable;
PowerTable(1:4)=WorkTable(1:4)*f;
