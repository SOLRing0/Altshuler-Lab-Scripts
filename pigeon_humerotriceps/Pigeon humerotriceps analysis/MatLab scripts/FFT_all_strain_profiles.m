%run an FFT analysis on muscle strain data and output the dominant frequency values
% plot strain profiles and FFT results 

clear all

%% set working directory
DirName=uigetdir; % select 'Experiments' folder ('Summarized data spreadsheets' folder for testing scripts without running all data processing scripts)
cd(DirName); % change working directory to experiments folder

min_PSD = 1.4180e-09; % max power of Tie  %6.0012e-10; pk 2 of Tie % 4.1655e-10; pk 3 of Tie
% determined from FFT averaged for each bird, where BirdTie
% had the lowest max power--cutoff chosen to include the main peak of
% BirdTie's FFT output
%% START in situ strain analysis %%
%% get all in situ strain data csv files
insitu_files= dir(fullfile(DirName, '*in situ*.csv'));

for a=1:length(insitu_files)
    filename=insitu_files(a).name;
    % get strain frequency
    stix(:,a)=strfind(filename,'-');
    stfreq(a)=str2double(filename(stix(1,a)+1:stix(2,a)-1)); %'+1' excludes the preceding '-' and '-1' excludes the anteceding '-'
end
%% loop to read in and stack strain data each trial of each experiment
for j=1:length(insitu_files)
    filename=insitu_files(j).name;
    insitu(j)= importinsitustr(filename, 2);  % startRow = 2 to skip column names
end

%%  loop to run FFT analysis on all trials and overplot FFT results from each experiment--sorted by frequency 
for k=1:length(stfreq)
    Fs = 10000;
    if unique(insitu(k).freq)==6.1
        
        for i=1:length(unique(insitu(k).Rep_num))
        Position = insitu(k).Position([find(insitu(k).Rep_num==i)]);
        [FFTdat(i), max_freq(i), max_power(i)]=FFTanalysis(Position,Fs,min_PSD);
        
        max_freq6(i)= max_freq(i);
        max_power6(i)= max_power(i);
        power6(:,i)=FFTdat(i).power;
        f6(i,:)=FFTdat(i).f;
        locs6(:,i)=FFTdat(i).locs;
        pks6(:,i)=FFTdat(i).pks;
        ffpow6(:,i)=FFTdat(i).ffpow;

    % plot strain from one trial and FFT results for all 6.1 Hz trials    
    % Time multiplied by 1000 to convert to milliseconds
    % minimum of time subtracted so that all strain profiles have same
    % start point
    subplot(7,1,3);
    plot(((insitu(k).Time([find(insitu(k).Rep_num==1)]))-min(insitu(k).Time([find(insitu(k).Rep_num==1)])))*1000,insitu(k).Position([find(insitu(k).Rep_num==1)]), '-g','LineWidth',1);
    xlim([0,1400]); 
    axis off
    ylim([-4,4]);
    box off
    hold on
    % plot FFT for 6.1 Hz
    subplot(7,1,7);
    plot(f6(i,:), power6(:,i), '-g','LineWidth',1);
    xlim([2,50]); % minimum set at 2 to exclude DC and noise if present
    ylim([0,6*10^-8]);
    box off
    hold on
        end
    elseif unique(insitu(k).freq)==8.6
        
        for i=1:length(unique(insitu(k).Rep_num))
        Position = insitu(k).Position([find(insitu(k).Rep_num==i)]);
        [FFTdat(i), max_freq(i), max_power(i)]=FFTanalysis(Position,Fs,min_PSD);
        
        max_freq8(i)= max_freq(i);
        max_power8(i)= max_power(i);
        power8(:,i)=FFTdat(i).power;
        f8(i,:)=FFTdat(i).f;
        locs8(:,i)=FFTdat(i).locs;
        pks8(:,i)=FFTdat(i).pks;
        ffpow8(:,i)=FFTdat(i).ffpow;
        
    % plot strain from one trial and FFT results for all 8.6 Hz trials
    % Time multiplied by 1000 to convert to milliseconds
    % minimum of time subtracted so that all strain profiles have same
    % start point
    subplot(7,1,4);
    plot(((insitu(k).Time([find(insitu(k).Rep_num==1)]))-min(insitu(k).Time([find(insitu(k).Rep_num==1)])))*1000,insitu(k).Position([find(insitu(j).Rep_num==1)]), '-k','LineWidth',1);
    xlim([0,1400]);
    axis off
    ylim([-4,4]);
    box off
    hold on
    % plot FFT for 8.6 Hz on all FFT plots
    subplot(7,1,6);
    plot(f8(i,:), power8(:,i), '-k','LineWidth',1);
    xlim([2,50]); % minimum set at 2 to exclude DC and noise if present
    box off
    hold on
    
    subplot(7,1,7);
    plot(f8(i,:), power8(:,i), '-k','LineWidth',1);
    xlim([2,50]); % minimum set at 2 to exclude DC and noise if present
    box off
    hold on
        end
    elseif unique(insitu(k).freq)==10.1
        
        for i=1:length(unique(insitu(k).Rep_num))
        Position = insitu(k).Position([find(insitu(k).Rep_num==i)]);
        [FFTdat(i), max_freq(i), max_power(i)]=FFTanalysis(Position,Fs,min_PSD);
        
        max_freq10(i)= max_freq(i);
        max_power10(i)= max_power(i);
        power10(:,i)=FFTdat(i).power;
        f10(i,:)=FFTdat(i).f;
        locs10(:,i)=FFTdat(i).locs;
        pks10(:,i)=FFTdat(i).pks;
        ffpow10(:,i)=FFTdat(i).ffpow;
    % plot strain from one trial and FFT results for all 10.1 Hz trials          
    % Time multiplied by 1000 to convert to milliseconds 
    % minimum of time subtracted so that all strain profiles have same
    % start point
    subplot(7,1,5);
    plot(((insitu(k).Time([find(insitu(k).Rep_num==1)]))-min(insitu(k).Time([find(insitu(k).Rep_num==1)])))*1000,insitu(k).Position([find(insitu(k).Rep_num==1)]), '-b','LineWidth',1);
    xlim([0,1400]);
    axis off
    plot([0; 0], [0; 0], '-k', [0; 100], [-4; -4], '-k', 'LineWidth', 1)   %x scale bar with labels at y=-4
    % add the scale bar labels
    text(25,-5,'100 ms','FontSize',8)
    ylim([-4,4]);
    box off
    hold on
    % Plot FFT for 10.1 Hz
    subplot(7,1,7);
    plot(f10(i,:), power10(:,i), '-b','LineWidth',1);
    % add axis labels to bottom panel
    xlim([2,50]); % minimum set at 2 to exclude DC and noise if present
    ylim([0,6*10^-8]);
    ylabel('Power');
    xlabel('Frequency (Hz)');
    box off
    hold on
        end
    else
        error('invalid in situ frequency')
    end
hold on


% total harmonic distortion in units of power for all 'peaks'
% harmonics of all 'fundamental' freqs
% locs (determined from minimum PSD cutoff) indicate that 6.1 Hz has 3 main 
% frequencies, 8.6 Hz has one, and 10.1 Hz has 2
       
    if unique(insitu(k).freq)==6.1
        % find 3 peaks 
        for i=1:length(unique(insitu(k).Rep_num))  
        % pk 1
        harmfreq6_1(:,i) = (f6(i,(rem(f6(i,:),pks6(1,i)) == 0 & f6(i,:)~=0)))';
        harmpow6_1(:,i) = power6((rem(f6(i,:),pks6(1,i)) == 0 & f6(i,:)~=0),i);
        % pk 2
        harmfreq6_2(:,i) = (f6(i,(rem(f6(i,:),pks6(2,i)) == 0 & f6(i,:)~=0)))';
        harmpow6_2(:,i) = power6((rem(f6(i,:),pks6(2,i)) == 0 & f6(i,:)~=0),i);
        % pk 3
        harmfreq6_3(:,i) = (f6(i,(rem(f6(i,:),pks6(3,i)) == 0 & f6(i,:)~=0)))';
        harmpow6_3(:,i) = power6((rem(f6(i,:),pks6(3,i)) == 0 & f6(i,:)~=0),i);
        % total harmonic distortion in units of power
        sum_pow6(i) = sum(harmpow6_1(:,i))+sum(harmpow6_2(:,i))+sum(harmpow6_3(:,i));
        pc6(:,i) = (ffpow6(:,i)/sum_pow6(i))*100;
        end
    elseif unique(insitu(k).freq)==8.6
        % find one peak
        for i=1:length(unique(insitu(k).Rep_num))  
        harmfreq8(:,i) = (f8(i,(rem(f8(i,:),pks8(1,i)) == 0 & f8(i,:)~=0)))';
        harmpow8(:,i) = power8((rem(f8(i,:),pks8(1,i)) == 0 & f8(i,:)~=0),i);
        % total harmonic distortion in units of power
        sum_pow8(i) = sum(harmpow8(:,i));
        pc8(:,i) = (ffpow8(:,i)/sum_pow8(i))*100;
        end
    elseif unique(insitu(k).freq)==10.1
        % find 2 peaks
        for i=1:length(unique(insitu(k).Rep_num))
        % pk 1
        harmfreq10_1(:,i) = (f10(i,(rem(f10(i,:),pks10(1,i)) == 0 & f10(i,:)~=0)))';
        harmpow10_1(:,i) = power10((rem(f10(i,:),pks10(1,i)) == 0 & f10(i,:)~=0),i);
        % pk 2
        harmfreq10_2(:,i) = (f10(i,(rem(f10(i,:),pks10(2,i)) == 0 & f10(i,:)~=0)))';
        harmpow10_2(:,i) = power10((rem(f10(i,:),pks10(2,i)) == 0 & f10(i,:)~=0),i);      
        % total harmonic distortion in units of power
        sum_pow10(i) = sum(harmpow10_1(:,i))+sum(harmpow10_2(:,i));
        pc10(:,i) = (ffpow10(:,i)/sum_pow10(i))*100;
        end
    else
        error('invalid strain frequency')
    end
        end

% calculate means of freqs and percentages across all trials for 6.1 Hz
 mean_pks6 = [mean(pks6(1,:)), mean(pks6(2,:)), mean(pks6(3,:))]';
 mean_pc6 = [mean(pc6(1,:)), mean(pc6(2,:)), mean(pc6(3,:))]';
 % calculate the weighted mean to get the fundamental frequency
 fund_freq6= repelem(sum(mean_pks6.*(mean_pc6/100)),length(mean_pc6))';
 % write a column to identify which strain trajectory the data are from
 strain_traj6 = repelem(string('6.1 Hz in situ'), length(mean_pc6));
% calculate means of freqs and percentages across all trials for 8.6 Hz
 mean_pks8 = (mean(pks8(1,:)))';
 mean_pc8 = (mean(pc8(1,:)))';
 % calculate the weighted mean to get the fundamental frequency
 fund_freq8= repelem(sum(mean_pks8.*(mean_pc8/100)),length(mean_pc8))';
 % write a column to identify which strain trajectory the data are from
 strain_traj8 = repelem(string('8.6 Hz in situ'), length(mean_pc8));
 
 % calculate means of freqs and percentages across all trials for 10.1 Hz
 mean_pks10 = [mean(pks10(1,:)), mean(pks10(2,:))]';
 mean_pc10 = [mean(pc10(1,:)), mean(pc10(2,:))]';
 % calculate the weighted mean to get the fundamental frequency
 fund_freq10= repelem(sum(mean_pks10.*(mean_pc10/100)),length(mean_pc10))';
 % write a column to identify which strain trajectory the data are from
 strain_traj10 = repelem(string('10.1 Hz in situ'), length(mean_pc10));

% create tables of the means for each frequency treatment
pc_f6 = table(strain_traj6',mean_pks6,mean_pc6,fund_freq6, 'VariableNames',{'strain_traj','freq','pc','fund_freq'});
pc_f8 = table(strain_traj8',mean_pks8,mean_pc8,fund_freq8, 'VariableNames',{'strain_traj','freq','pc','fund_freq'});
pc_f10 = table(strain_traj10',mean_pks10,mean_pc10,fund_freq10, 'VariableNames',{'strain_traj','freq','pc','fund_freq'});

% combine 3 in situ tables into one to later be combined with in vivo data
FFTinsitusummary = vertcat(pc_f6,pc_f8,pc_f10);

% write in situ summary table to csv
writetable(FFTinsitusummary, 'FFT analysis in situ summary data.csv');

     %% START in vivo strain analysis %%
%     %% get all in vivo strain data csv files
%     invivo_files=dir(fullfile(DirName, '*in vivo*.csv'));
% 
%     for a=1:length(invivo_files)
%         filename=invivo_files(a).name;
%         % get subject ID
%         idx(:,a)=strfind(filename,' ');
%         Bird_ID(a)=string(filename([(1:idx(1,a)-1)]));
%     end
% 
%     %% loop to read in strain data for in vivo
%     for j=1:length(invivo_files)
%         filename=invivo_files(j).name;
%         invivo(j)= importinvivostr(filename, 2);  % startRow = 2 to skip column names
%     end
%     % Used to find all frequencies with power above the minimum threshold
%     %%  loop to run FFT analysis on all in vivo trials and plot strain profiles and FFT results
%     for n=1:length(invivo)
%         if contains(unique(invivo(n).Bird), 'BirdTie') == 1
%             xt(n)=length(invivo(n).Position); % get the lengths of all trials from one bird
%         elseif contains(unique(invivo(n).Bird), 'BirdUuu') == 1
%             xu(n-5)=length(invivo(n).Position); % get the lengths of all trials from one bird
%         elseif contains(unique(invivo(n).Bird), 'BirdWww') == 1
%             xw(n-10)=length(invivo(n).Position); % get the lengths of all trials from one bird
%         elseif contains(unique(invivo(n).Bird), 'BirdXxx') == 1
%             xx(n-15)=length(invivo(n).Position); % get the lengths of all trials from one bird
%         elseif contains(unique(invivo(n).Bird), 'BirdYyy') == 1
%             xy(n-20)=length(invivo(n).Position); % get the lengths of all trials from one bird
%         end
%     end
% 
%     for k=1:length(invivo)
%         Fs=5000;
% 
%         if contains(unique(invivo(k).Bird),'BirdTie') == 1 % choose files from one Bird
%             L=min(xt); % use length of shortest trial for FFTs for this Bird
%             Position = (invivo(k).Position)*10; % convert from cm to mm
%             Time = invivo(k).Time; % Feed time vector to function to be processed for plotting strains
%             [FFTadat(k)]=FFTamplitude(Position,Fs,L,Time);
%             fbt=FFTadat(k).f;
%             ampbt(:,k)=FFTadat(k).amplitude;
%             positionbt(:,k)=FFTadat(k).position; % position adjusted to center around 0
%             timebt(:,k)=FFTadat(k).time;
% 
%         elseif  contains(unique(invivo(k).Bird),'BirdUuu') == 1 % choose files from one Bird
%             L=min(xu); % use length of shortest trial for FFTs for this Bird
%             Position = (invivo(k).Position)*10; % convert from cm to mm
%             Time = invivo(k).Time; % Feed time vector to function to be processed for plotting strains
%             [FFTadat(k)]=FFTamplitude(Position,Fs,L,Time);
%             fbu=FFTadat(k).f;
%             ampbu(:,k-5)=FFTadat(k).amplitude;
%             positionbu(:,k-5)=FFTadat(k).position; % position adjusted to center around 0
%             timebu(:,k-5)=FFTadat(k).time;
% 
%         elseif contains(unique(invivo(k).Bird),'BirdWww') == 1 % choose files from one Bird
%             L=min(xw); % use length of shortest trial for FFTs for this Bird
%             Position = (invivo(k).Position)*10; % convert from cm to mm
%             Time = invivo(k).Time; % Feed time vector to function to be processed for plotting strains
%             [FFTadat(k)]=FFTamplitude(Position,Fs,L,Time);
%             fbw=FFTadat(k).f;
%             ampbw(:,k-10)=FFTadat(k).amplitude;
%             positionbw(:,k-10)=FFTadat(k).position; % position adjusted to center around 0
%             timebw(:,k-10)=FFTadat(k).time;
% 
%         elseif  contains(unique(invivo(k).Bird),'BirdXxx') == 1 % how to choose files from one Bird?
%             L=min(xx); % use length of shortest trial for FFTs for this Bird
%             Position = (invivo(k).Position)*10; % convert from cm to mm
%             Time = invivo(k).Time; % Feed time vector to function to be processed for plotting strains
%             [FFTadat(k)]=FFTamplitude(Position,Fs,L,Time);
%             fbx=FFTadat(k).f;
%             ampbx(:,k-15)=FFTadat(k).amplitude;
%             positionbx(:,k-15)=FFTadat(k).position; % position adjusted to center around 0
%             timebx(:,k-15)=FFTadat(k).time;
% 
%         elseif  contains(unique(invivo(k).Bird),'BirdYyy') == 1 % choose files from one Bird
%             L=min(xy); % use length of shortest trial for FFTs for this Bird
%             Position = (invivo(k).Position)*10; % convert from cm to mm
%             Time = invivo(k).Time; % Feed time vector to function to be processed for plotting strains
%             [FFTadat(k)]=FFTamplitude(Position,Fs,L,Time);
%             fby=FFTadat(k).f;
%             ampby(:,k-20)=FFTadat(k).amplitude;
%             positionby(:,k-20)=FFTadat(k).position; % position adjusted to center around 0
%             timeby(:,k-20)=FFTadat(k).time;
% 
%         end
%     end
% 
%     for m=1:length(invivo)
%         if rem(m,5)==0
%             if contains(unique(invivo(m).Bird),'BirdTie') == 1
%                 L=min(xt); % use length of shortest trial for FFTs for this Bird
%                 mean_amp=mean(ampbt,2); % compute mean amplitude by row (for one bird across all trials)
%                 f=fbt;
%                 [FFTpdat, max_freq, max_power]=FFTpow(Fs,min_PSD,L,f,mean_amp);
% 
%                 max_freqbt= max_freq;
%                 max_powerbt= max_power;
%                 powerbt=FFTpdat.power;
%                 locsbt=FFTpdat.locs;
%                 pksbt=FFTpdat.pks([1:end]);
%                 ffpowbt=FFTpdat.ffpow([1:end]);
% 
%             elseif  contains(unique(invivo(m).Bird),'BirdUuu') == 1 % choose files from one Bird
%                 L=min(xu); % use length of shortest trial for FFTs for this Bird
%                 mean_amp=mean(ampbu,2); % compute mean amplitude by row (for one bird across all trials)
%                 f=fbu;
%                 [FFTpdat, max_freq, max_power]=FFTpow(Fs,min_PSD,L,f,mean_amp);
% 
%                 max_freqbu= max_freq;
%                 max_powerbu= max_power;
%                 powerbu=FFTpdat.power;
%                 locsbu=FFTpdat.locs;
%                 pksbu=FFTpdat.pks([1:end]);
%                 ffpowbu=FFTpdat.ffpow([1:end]);
% 
%             elseif contains(unique(invivo(m).Bird),'BirdWww') == 1 % choose files from one Bird
%                 L=min(xw); % use length of shortest trial for FFTs for this Bird
%                 mean_amp=mean(ampbw,2); % compute mean amplitude by row (for one bird across all trials)
%                 f=fbw;
%                 [FFTpdat, max_freq, max_power]=FFTpow(Fs,min_PSD,L,f,mean_amp);
% 
%                 max_freqbw= max_freq;
%                 max_powerbw= max_power;
%                 powerbw=FFTpdat.power;
%                 locsbw=FFTpdat.locs;
%                 pksbw=FFTpdat.pks([1:end]);
%                 ffpowbw=FFTpdat.ffpow([1:end]);
% 
%             elseif  contains(unique(invivo(m).Bird),'BirdXxx') == 1 % choose files from one Bird
%                 L=min(xx); % use length of shortest trial for FFTs for this Bird
%                 mean_amp=mean(ampbx,2); % compute mean amplitude by row (for one bird across all trials)
%                 f=fbx;
%                 [FFTpdat, max_freq, max_power]=FFTpow(Fs,min_PSD,L,f,mean_amp);
% 
%                 max_freqbx= max_freq;
%                 max_powerbx= max_power;
%                 powerbx=FFTpdat.power;
%                 locsbx=FFTpdat.locs;
%                 pksbx=FFTpdat.pks([1:end]);
%                 ffpowbx=FFTpdat.ffpow([1:end]);
% 
%             elseif  contains(unique(invivo(m).Bird),'BirdYyy') == 1 % choose files from one Bird
%                 L=min(xy); % use length of shortest trial for FFTs for this Bird
%                 mean_amp=mean(ampby,2); % compute mean amplitude by row (for one bird across all trials)
%                 f=fby;
%                 [FFTpdat, max_freq, max_power]=FFTpow(Fs,min_PSD,L,f,mean_amp);
% 
%                 max_freqby= max_freq;
%                 max_powerby= max_power;
%                 powerby=FFTpdat.power;
%                 locsby=FFTpdat.locs;
%                 pksby=FFTpdat.pks([1:end]);
%                 ffpowby=FFTpdat.ffpow([1:end]);
% 
%                 % plot strain from BirdYyy, trial 2 for figure 3
%                 subplot(7,1,1);
%                 plot(timeby(:,2),smooth(positionby(:,2),'sgolay'), '-m','LineWidth',1);
%                 xlim([0,1400]);
%                 axis off
%                 ylim([-4,4]);
%                 box off
%                 hold on
% 
%                 % plot strain from BirdUuu, trial 4 for figure 3
%                 subplot(7,1,2);
%                 plot(timebu(:,4),smooth(positionbu(:,4),'sgolay'), '-y','LineWidth',1);
%                 xlim([0,1400]);
%                 axis off
%                 ylim([-4,4]);
%                 box off
%                 hold on
% 
%                 % plot FFT results averaged across all trials for each Bird
%                 % BirdTie
%                 subplot(7,1,6);
%                 plot(fbt, powerbt, '-c','LineWidth',1);
%                 xlim([2,50]); % minimum set at 2 to exclude DC and noise if present
%                 ylim([0,6*10^-8]); % set all y-axes to the same limits
%                 box off;
%                 hold on
% 
%                 % BirdUuu
%                 subplot(7,1,6);
%                 plot(fbu, powerbu, '-y','LineWidth',1);
%                 xlim([2,50]); % minimum set at 2 to exclude DC and noise if present
%                 ylim([0,6*10^-8]); % set all y-axes to the same limits
%                 box off;
%                 hold on
% 
%                 % BirdWww
%                 subplot(7,1,6);
%                 plot(fbw, powerbw, '-r','LineWidth',1);
%                 xlim([2,50]); % minimum set at 2 to exclude DC and noise if present
%                 ylim([0,6*10^-8]); % set all y-axes to the same limits
%                 box off;
%                 hold on
% 
%                 % BirdXxx
%                 subplot(7,1,6);
%                 plot(fbx, powerbx, '-g','LineWidth',1);
%                 xlim([2,50]); % minimum set at 2 to exclude DC and noise if present
%                 ylim([0,6*10^-8]); % set all y-axes to the same limits
%                 box off;
%                 hold on
% 
%                 % BirdYyy
%                 subplot(7,1,6);
%                 plot(fby, powerby, '-m','LineWidth',1);
%                 xlim([2,50]); % minimum set at 2 to exclude DC and noise if present
%                 ylim([0,6*10^-8]); % set all y-axes to the same limits
%                 box off;
%                 hold on
% 
%             end
%         end
%     end
% 
%     %% total harmonic distortion in units of power for all 'peaks'
%     % harmonics of all 'fundamental' freqs
%     % locs (determined from minimum PSD cutoff) indicate that BirdTie has 1
%     % frequency, BirdUuu has 6 main frequencies, and BirdWww, BirdXxx and
%     % BirdYyy have 8
%     %% Find peaks for BirdTie
%     % pk 1
%     harmfreqbt_1 = fbt(rem(fbt,pksbt(1)) == 0 & fbt~=0);
%     harmpowbt_1 = powerbt(rem(fbt,pksbt(1)) == 0 & fbt~=0);
%     sum_powbt = sum(harmpowbt_1);
% 
%     % calculate the percentages of each contributing frequency
%     pcbt = (ffpowbt/sum_powbt)*100;
%     % calculate the weighted mean to get the fundamental frequency
%     fund_freq_bt= repelem(sum(pksbt.*(pcbt/100)),length(pcbt))';
%     % write a column to identify which strain trajectory the data are from
%     strain_trajbt = repelem(string('In vivo Bird 1 (Tie)'), length(pcbt));
%     % create a table of those percentages and frequencies
%     pc_fbt = table(strain_trajbt',pksbt',pcbt',fund_freq_bt, 'VariableNames',{'strain_traj','freq','pc','fund_freq'});
% 
%     %% Find peaks for BirdUuu
%     % pk 1
%     harmfreqbu_1 = fbu(rem(fbu,pksbu(1)) == 0 & fbu~=0);
%     harmpowbu_1 = powerbu(rem(fbu,pksbu(1)) == 0 & fbu~=0);
%     % pk2
%     harmfreqbu_2 = fbu(rem(fbu,pksbu(2)) == 0 & fbu~=0);
%     harmpowbu_2 = powerbu(rem(fbu,pksbu(2)) == 0 & fbu~=0);
%     % pk3
%     harmfreqbu_3 = fbu(rem(fbu,pksbu(3)) == 0 & fbu~=0);
%     harmpowbu_3 = powerbu(rem(fbu,pksbu(3)) == 0 & fbu~=0);
%     % pk4
%     harmfreqbu_4 = fbu(rem(fbu,pksbu(4)) == 0 & fbu~=0);
%     harmpowbu_4 = powerbu(rem(fbu,pksbu(4)) == 0 & fbu~=0);
%     % pk5
%     harmfreqbu_5 = fbu(rem(fbu,pksbu(5)) == 0 & fbu~=0);
%     harmpowbu_5 = powerbu(rem(fbu,pksbu(5)) == 0 & fbu~=0);
%     % pk6
%     harmfreqbu_6 = fbu(rem(fbu,pksbu(6)) == 0 & fbu~=0);
%     harmpowbu_6 = powerbu(rem(fbu,pksbu(6)) == 0 & fbu~=0);
%     sum_powbu = sum(harmpowbu_1) + sum(harmpowbu_2) + sum(harmpowbu_3) + sum(harmpowbu_4) + sum(harmpowbu_5) + sum(harmpowbu_6);
% 
%     % calculate the percentages of each contributing frequency
%     pcbu = (ffpowbu/sum_powbu)*100;
%     % calculate the weighted mean to get the fundamental frequency
%     fund_freq_bu= repelem(sum(pksbu.*(pcbu/100)),length(pcbu))';
%     % write a column to identify which strain trajectory the data are from
%     strain_trajbu = repelem(string('In vivo Bird 2 (Uuu)'), length(pcbu));
%     % create a table of those percentages and frequencies
%     pc_fbu = table(strain_trajbu',pksbu',pcbu',fund_freq_bu, 'VariableNames',{'strain_traj','freq','pc','fund_freq'});
% 
%     %% Find peaks for BirdWww
%     % pk 1
%     harmfreqbw_1 = fbw(rem(fbw,pksbw(1)) == 0 & fbw~=0);
%     harmpowbw_1 = powerbw(rem(fbw,pksbw(1)) == 0 & fbw~=0);
%     % pk2
%     harmfreqbw_2 = fbw(rem(fbw,pksbw(2)) == 0 & fbw~=0);
%     harmpowbw_2 = powerbw(rem(fbw,pksbw(2)) == 0 & fbw~=0);
%     % pk3
%     harmfreqbw_3 = fbw(rem(fbw,pksbw(3)) == 0 & fbw~=0);
%     harmpowbw_3 = powerbw(rem(fbw,pksbw(3)) == 0 & fbw~=0);
%     % pk4
%     harmfreqbw_4 = fbw(rem(fbw,pksbw(4)) == 0 & fbw~=0);
%     harmpowbw_4 = powerbw(rem(fbw,pksbw(4)) == 0 & fbw~=0);
%     % pk5
%     harmfreqbw_5 = fbw(rem(fbw,pksbw(5)) == 0 & fbw~=0);
%     harmpowbw_5 = powerbw(rem(fbw,pksbw(5)) == 0 & fbw~=0);
%     % pk6
%     harmfreqbw_6 = fbw(rem(fbw,pksbw(6)) == 0 & fbw~=0);
%     harmpowbw_6 = powerbw(rem(fbw,pksbw(6)) == 0 & fbw~=0);
%     % pk7
%     harmfreqbw_7 = fbw(rem(fbw,pksbw(7)) == 0 & fbw~=0);
%     harmpowbw_7 = powerbw(rem(fbw,pksbw(7)) == 0 & fbw~=0);
% 
%     % pk8
%     harmfreqbw_8 = fbw(rem(fbw,pksbw(8)) == 0 & fbw~=0);
%     harmpowbw_8 = powerbw(rem(fbw,pksbw(8)) == 0 & fbw~=0);
%     sum_powbw = sum(harmpowbw_1) + sum(harmpowbw_2) + sum(harmpowbw_3) + sum(harmpowbw_4) + sum(harmpowbw_5) + sum(harmpowbw_6) + sum(harmpowbw_7) + sum(harmpowbw_8);
% 
%     % calculate the percentages of each contributing frequency
%     pcbw = (ffpowbw/sum_powbw)*100;
%     % calculate the weighted mean to get the fundamental frequency
%     fund_freq_bw= repelem(sum(pksbw.*(pcbw/100)),length(pcbw))';
%     % write a column to identify which strain trajectory the data are from
%     strain_trajbw = repelem(string('In vivo Bird 3 (Www)'), length(pcbw));
%     % create a table of those percentages and frequencies
%     pc_fbw = table(strain_trajbw',pksbw',pcbw',fund_freq_bw, 'VariableNames',{'strain_traj','freq','pc','fund_freq'});
% 
%     %% Find peaks for BirdXxx
%     % pk 1
%     harmfreqbx_1 = fbx(rem(fbx,pksbx(1)) == 0 & fbx~=0);
%     harmpowbx_1 = powerbx(rem(fbx,pksbx(1)) == 0 & fbx~=0);
%     % pk2
%     harmfreqbx_2 = fbx(rem(fbx,pksbx(2)) == 0 & fbx~=0);
%     harmpowbx_2 = powerbx(rem(fbx,pksbx(2)) == 0 & fbx~=0);
%     % pk3
%     harmfreqbx_3 = fbx(rem(fbx,pksbx(3)) == 0 & fbx~=0);
%     harmpowbx_3 = powerbx(rem(fbx,pksbx(3)) == 0 & fbx~=0);
%     % pk4
%     harmfreqbx_4 = fbx(rem(fbx,pksbx(4)) == 0 & fbx~=0);
%     harmpowbx_4 = powerbx(rem(fbx,pksbx(4)) == 0 & fbx~=0);
%     % pk5
%     harmfreqbx_5 = fbx(rem(fbx,pksbx(5)) == 0 & fbx~=0);
%     harmpowbx_5 = powerbx(rem(fbx,pksbx(5)) == 0 & fbx~=0);
%     % pk6
%     harmfreqbx_6 = fbx(rem(fbx,pksbx(6)) == 0 & fbx~=0);
%     harmpowbx_6 = powerbx(rem(fbx,pksbx(6)) == 0 & fbx~=0);
%     % pk7
%     harmfreqbx_7 = fbx(rem(fbx,pksbx(7)) == 0 & fbx~=0);
%     harmpowbx_7 = powerbx(rem(fbx,pksbx(7)) == 0 & fbx~=0);
%     % pk8
%     harmfreqbx_8 = fbx(rem(fbx,pksbx(8)) == 0 & fbx~=0);
%     harmpowbx_8 = powerbx(rem(fbx,pksbx(8)) == 0 & fbx~=0);
%     sum_powbx = sum(harmpowbx_1) + sum(harmpowbx_2) + sum(harmpowbx_3) + sum(harmpowbx_4) + sum(harmpowbx_5) + sum(harmpowbx_6) + sum(harmpowbx_7) + sum(harmpowbx_8);
% 
%     % calculate the percentages of each contributing frequency
%     pcbx = (ffpowbx/sum_powbx)*100;
%     % calculate the weighted mean to get the fundamental frequency
%     fund_freq_bx= repelem(sum(pksbx.*(pcbx/100)),length(pcbx))';
%     % write a column to identify which strain trajectory the data are from
%     strain_trajbx = repelem(string('In vivo Bird 4 (Xxx)'), length(pcbx));
%     % create a table of those percentages and frequencies
%     pc_fbx = table(strain_trajbx',pksbx',pcbx',fund_freq_bx, 'VariableNames',{'strain_traj','freq','pc','fund_freq'});
% 
%     %% Find peaks for BirdYyy
%     % pk 1
%     harmfreqby_1 = fby(rem(fby,pksby(1)) == 0 & fby~=0);
%     harmpowby_1 = powerby(rem(fby,pksby(1)) == 0 & fby~=0);
%     % pk2
%     harmfreqby_2 = fby(rem(fby,pksby(2)) == 0 & fby~=0);
%     harmpowby_2 = powerby(rem(fby,pksby(2)) == 0 & fby~=0);
%     % pk3
%     harmfreqby_3 = fby(rem(fby,pksby(3)) == 0 & fby~=0);
%     harmpowby_3 = powerby(rem(fby,pksby(3)) == 0 & fby~=0);
%     % pk4
%     harmfreqby_4 = fby(rem(fby,pksby(4)) == 0 & fby~=0);
%     harmpowby_4 = powerby(rem(fby,pksby(4)) == 0 & fby~=0);
%     % pk5
%     harmfreqby_5 = fby(rem(fby,pksby(5)) == 0 & fby~=0);
%     harmpowby_5 = powerby(rem(fby,pksby(5)) == 0 & fby~=0);
%     % pk6
%     harmfreqby_6 = fby(rem(fby,pksby(6)) == 0 & fby~=0);
%     harmpowby_6 = powerby(rem(fby,pksby(6)) == 0 & fby~=0);
%     % pk7
%     harmfreqby_7 = fby(rem(fby,pksby(7)) == 0 & fby~=0);
%     harmpowby_7 = powerby(rem(fby,pksby(7)) == 0 & fby~=0);
%     % pk8
%     harmfreqby_8 = fby(rem(fby,pksby(8)) == 0 & fby~=0);
%     harmpowby_8 = powerby(rem(fby,pksby(8)) == 0 & fby~=0);
%     sum_powby = sum(harmpowby_1) + sum(harmpowby_2) + sum(harmpowby_3) + sum(harmpowby_4) + sum(harmpowby_5) + sum(harmpowby_6) + sum(harmpowby_7) + sum(harmpowby_8);
% 
%     % calculate the percentages of each contributing frequency
%     pcby = (ffpowby/sum_powby)*100;
%     % calculate the weighted mean to get the fundamental frequency
%     fund_freq_by= repelem(sum(pksby.*(pcby/100)),length(pcby))';
%     % write a column to identify which strain trajectory the data are from
%     strain_trajby = repelem(string('In vivo Bird 5 (Yyy)'), length(pcby));
%     % create a table of those percentages and frequencies
%     pc_fby = table(strain_trajby',pksby',pcby',fund_freq_by, 'VariableNames',{'strain_traj','freq','pc','fund_freq'});
%     % combine 5 tables into one to later be combined with in situ data and to one csv file
%     FFTinvivosummary = vertcat(pc_fbt,pc_fbu,pc_fbw,pc_fbx,pc_fby);
%     % combine in vivo and insitu tables into one to be written into one csv file
%     FFTsummary = vertcat(FFTinvivosummary,FFTinsitusummary);
%     % write in situ summary table to csv
%     writetable(FFTsummary, 'FFT analysis summary data.csv');     

