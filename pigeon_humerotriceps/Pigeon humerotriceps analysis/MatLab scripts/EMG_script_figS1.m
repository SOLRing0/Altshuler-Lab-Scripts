%% Script to plot EMG data for Figure S1

clear
close all
%% Panel A

DirName=uigetdir; % select 'Experiments' folder
cd(DirName);
FolderName='2017-04-10'; % selects trial date EMGs in panel A were generated from
cd(FolderName);
PathName='8.6 Hz - 69pc cycle'; % select treatment folder within trial date folder
cd(PathName);
% Calls stimulus phase treatment that plots were generated from
FileName='-HT_GR4_varyPhase--30%B0_1.ddf';

% get trial frequency
f=extract_frequency(FileName,15,15);

% set line number for endRow using frequency setting
if f==8.6 
    numl=8290; 
    else 
      error('invalid frequency and endRow value for figure S1 plots')
end

[Sample,position,force,AI2,AI3,AI4,AI5,AI6,AI7,AO0,AO1,Stim] = importWLdata(FileName,20,numl);  

% plotting parameters 
% subplot stimulus signal
subplot(4,1,1);
plot(Sample(3600:5900),AI2(3600:5900),'-k');
xlabel('Time (10^-^4 s)');
xlim([3600,5900]);
ylim([3,4]);
axis off
set(gca,'box','off')
hold on

% subplot EMG trace
subplot(4,1,2);
plot(Sample(3600:5900),AI3(3600:5900)); %value range selects cycles 4 & 5
xlim([3600,5900]);
ylim([-13,11]);
axis off
set(gca,'box','off')
hold on

% plot y axis scale bar
plot([3600;3600], [-13;-2], '-k', [0;0], [0;0], '-k', 'LineWidth', 1)
% plot x axis scale bar
plot([0; 0], [0; 0], '-k', [3600; 4100], [-13; -13], '-k', 'LineWidth', 1)   
hold off
% add y scale bar label
text(3400,-7.5,'0.5 V');
% add x scale bar label
text(3750,-16, '50 ms');

%% Panel B

cd(DirName); % return to 'Experiments' folder
FolderName='2017-06-06'; % selects trial date EMGs in panel A were generated from
cd(FolderName);
PathName='8.6 Hz - 50pc cycle'; % select treatment folder within trial date folder
cd(PathName);
% Calls stimulus phase treatment that plots were generated from
FileName='-HT_GR4_varyPhase_8.6Hz_50pcstim-50%0_1.ddf';

% get trial frequency
f=extract_frequency(FileName,15,15);

% set line number for endRow using frequency setting
if f==8.6 
    numl=8290; 
    else 
      error('invalid frequency and endRow value for figure S1 plots')
end

[Sample,position,force,AI2,AI3,AI4,AI5,AI6,AI7,AO0,AO1,Stim] = importWLdata(FileName,20,numl);  

% plotting parameters 
% subplot stimulus signal
subplot(4,1,3);
plot(Sample(3350:5650),AI2(3350:5650),'-k');
xlabel('Time (10^-^4 s)');
xlim([3350,5650]);
ylim([3,4]);
axis off
set(gca,'box','off')
hold on

% subplot EMG trace
subplot(4,1,4);
plot(Sample(3350:5650),AI3(3350:5650),'-r'); %value range selects cycles 4 & 5
xlim([3350,5650]);
ylim([-13,11]);
axis off
set(gca,'box','off')
hold on