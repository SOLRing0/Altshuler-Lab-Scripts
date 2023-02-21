% code for running a single trial
clear

% script to generate the sample work loop plots from one control cycle for 
% Figure 4.
DirName=uigetdir; % select 'Experiments' folder
cd(DirName);
PathName='2017-06-27'; % selects trial date that control work loop plots were generated from
cd(PathName);
% Calls control work loop file that plots were generated from
FileName='WL012.ddf';

%Program to analyze workloop data, calculating the work per cycle

% get trial frequency
f=extract_frequency(FileName,15,15);

% set line number for endRow using frequency setting
if f==8.6 
    numl=8290; 
    else 
      error('invalid frequency and endRow value for figure 4 plots')
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
    else
      error('invalid frequency and x-limit for figure 4 plots')
    end
% plot length and force traces on the same plot with separate y-axes
[AX, H1, H2]=plotyy(Sample(1:xend),(position(1:xend))/GR,Sample(1:xend),(force(1:xend))*GR);
xlabel('Time (10^-^4 s)'); % label x-axis
set(get(AX(1),'Ylabel'),'String','\Delta Length (mm)'); % label length change y-axis
set(get(AX(2),'Ylabel'),'String','Force (N)'); % label force y-axis
set(AX(1),'YLim',[-2.5,2.5]); % sets length change axis limits 
set(AX(2),'YLim',[0,25]); % sets Force axis limits

%% crop and stack data
% set startdex search parameters -- end point of search range
if f==8.6 
    endpt=1200; 
    else 
      error('invalid frequency and startdex parameter value for figure 4 plots')
end

startdex=find(position(1:endpt)==max(position(1:endpt)));
enddex=find(position(end-500:end)==max(position(end-500:end)))+length(position)-500;

% set EndLoop search parameters
% search range start
if f==8.6 
    ELs=1000; 
    else 
      error('invalid frequency and EndLoop range start for figure 4 plots')
end
% search range end
if f==8.6 
    ELe=1500; 
    else 
      error('invalid frequency and EndLoop range end for figure 4 plots')
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

% Work loop plots (force vs length) for cycles 2 through 5
subplot(1,4,3); hold on
plot(Lengths(:,(2:5)),Forces(:,(2:5)));
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

% Plot of net work per cycle for cycles 2 through 5
subplot(1,4,4);
plot(2:5,work(2:5),'bo');
xlabel('cycle number'); ylabel('Work (mJ)');
xlim([1,6]); % set x-axis limits--max=6 so that last cycle isn't cut off
ylim([-2,4]); % set y-axis limits

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
