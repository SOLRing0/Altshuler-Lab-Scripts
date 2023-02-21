%% program to analyze all data within a folder

clear
close all

GR=4;% GR = inverted fraction (1/4) of hole distance from end of servo motor arm  

DirName=uigetdir; % select 'Experiments' folder
cd(DirName); % change working directory to experiments folder
% Get the number of directories (excludes files) within the experiment folder
% This gives the number of experiments--number of times to run analysis loop 
all_files = dir;
all_dir = all_files([all_files(:).isdir]);
Q = numel(all_dir);
 
delete('power_data.csv'); % delete power_data.csv so as not to append duplicates to existing data

for q = 3:Q % loop to move through different experiment folders
cd(DirName); % make sure to go back to experiments folder before loop starts
folderList=dir(pwd);    % gets a list of all the experiments folders
cd(folderList(q).name); % Sets next experiment folder as working directory
FolderName=dir(pwd);    % gets a list of treatment folders for a given experiment

   % sets conditions to determine # of treatments (and loops) per
   % experiment
    if folderList(q).name=='2017-04-05' % has two different treatments
        h=4; 
        elseif folderList(q).name=='2017-04-12' % has two different treatments
        h=4;
        else % the rest have only one treatment
        h=3;
    end
    
    for i = 3:h % loop to move through treatment folders within a experiment folder
    cd(FolderName(i).name);
    PathName=pwd;
    A=dir(fullfile(PathName,'*.ddf')); % saves data from all .ddf files into a structure called A 

    scrsz = get(0,'ScreenSize'); % makes plots fullscreen
    
        % loops through all trials for a given treatment
        for a=1:(length(A))
        FileName=A(a).name;
        % get trial frequency
        f(a)=extract_frequency(FileName,15,15); % extracts frequency from cell 15 because this is where the parameter settings are saved to each ddf file and frequency of the sine
                                                % wave is always located in cell 15
        pData(a)=AnalyseWLTri(FileName,GR,f);
        % get trial date 
        Date(a)=string(folderList(q).name);
        end
        
    % sort pData by ascending phase
    [sortedphase, phaseidx] = sort([pData.phase]);
    sorted_pData = pData(phaseidx);
    
     for j=1:(length(A)/2)
     
        % get phases from file name
        stimPhase(j)=sorted_pData(j).phase;
    
        Work5(j)=sorted_pData(j).work5;
        Work3(j)=sorted_pData(j).work3;
        Work2(j)=sorted_pData(j).work2;
        Work1(j)=sorted_pData(j).work1;
     end 
     
     for j=((length(A)/2)+1):length(A)
      
        % get phases from file name
        stimPhase(j)=sorted_pData(j).phase;
    
        Work5(j)=sorted_pData(j).work5;
        Work3(j)=sorted_pData(j).work3;
        Work2(j)=sorted_pData(j).work2;
        Work1(j)=sorted_pData(j).work1;    
     end
%% PLOTTING SECTION FOR VISUAL INSPECTION--UNCOMMENT TO VIEW PLOTS
%         % name figures with date and treatment and hold them onscreen
%         name = sprintf('%s',string(folderList(q).name),' - ', string(FolderName(i).name), ' - work loops & net work per cycle - negative phases');
%         fgr=figure();
%         set(fgr,'Name',name,'NumberTitle','off');
%         set(fgr,'Unit','normalized','Position',[0,0,1,1]);
%         hold off 
%     
%         for j=1:(length(A)/2)
%      
%         %   force plots (workloops)
%         subplot(2,length(A)/2,j);
%         plot(sorted_pData(j).Lengths,sorted_pData(j).Forces);
%         ylim([0,30]);
%         xlim([-3,3]);
%         title(stimPhase(j));
%     
% %       work plots
%         subplot(2,length(A)/2,j+(length(A)/2));
%         plot(1:1:5,sorted_pData(j).work,'bo');
%         ylim([-40,10]);
%         xlim([0,6]);
%       
%         end 
%         
%     subplot(2,length(A)/2,1);
%     ylabel('force (N)');
%     subplot(2,length(A)/2,1+(length(A)/2));
%     ylabel('work (mJ)');
% 
%         % name figures with date and treatment and hold them onscreen
%         name1 = sprintf('%s',string(folderList(q).name),' - ', string(FolderName(i).name), ' - work loops & net work per cycle - positive phases');
%         fgr=figure();
%         set(fgr,'Name',name1,'NumberTitle','off');
%         set(fgr,'Unit','normalized','Position',[0,0,1,1]);
%         hold off   
%         
%         for j=((length(A)/2)+1):length(A)
%     
%     %   force plots (workloops)
%         subplot(2,length(A)/2,j-(length(A)/2));
%         plot(sorted_pData(j).Lengths,sorted_pData(j).Forces);
%         ylim([0,30]);
%         xlim([-3,3]);
%         title(stimPhase(j));
%     
%     %   work plots    
%         subplot(2,length(A)/2,j);
%         plot(1:1:5,sorted_pData(j).work,'bo');
%         ylim([-90,0]);
%         xlim([0,6]);
%         end 
% 
%     subplot(2,length(A)/2,1);
%     ylabel('force (N)');
%     subplot(2,length(A)/2,1+(length(A)/2));
%     ylabel('work (mJ)');
% 
%     % Work vs Phase plot
%     % name figures with date and treatment and holds them onscreen
%     name2 = sprintf('%s',string(folderList(q).name),' - ', string(FolderName(i).name), ' - Work vs Phase');
%     fgr=figure();
%     set(fgr,'Name',name2,'NumberTitle','off');
%     set(fgr,'Unit','normalized','Position',[0,0,1,1]);
%     plot(i+2:(length(i+3)));
%     hold on 
%     
%     plot(stimPhase,Work5,'bo');
%     plot(stimPhase,Work3,'go');
%     plot(stimPhase,Work2,'ro');
%     plot(stimPhase,Work1,'mo');
%     plot(stimPhase,Work1,'m-');
% 
%     xlabel('phase (%)');ylabel('work (mJ)')
%     plot(xlim,[0,0],'k-'); 
%     legend('all 5','top 3','top 2','top 1');
%     % Power vs Phase plot
% 
%     % name figures with date and treatment and holds them onscreen
%     name2 = sprintf('%s',string(folderList(q).name),' - ', string(FolderName(i).name), ' - Power vs Phase');
%     fgr=figure();
%     set(fgr,'Name',name2,'NumberTitle','off');
%     set(fgr,'Unit','normalized','Position',[0,0,1,1]);
%     plot(i+3:(length(i+4)));
%     hold on 
%     
%     plot(stimPhase,Work5.*f,'bo');
%     plot(stimPhase,Work3.*f,'go');
%     plot(stimPhase,Work2.*f,'ro');
%     plot(stimPhase,Work1.*f,'mo');
%     plot(stimPhase,Work1.*f,'m-');
% 
%     xlabel('phase (%)');ylabel('power (mW)')
%     plot(xlim,[0,0],'k-'); 
%     legend('all 5','top 3','top 2','top 1');

    %% summary tables
        
    % Net work 
    WorkTable=zeros(length(stimPhase),6);
    WorkTable(:,1)=f;
    WorkTable(:,2)=stimPhase';
    WorkTable(:,3)=Work1';
    WorkTable(:,4)=Work2';
    WorkTable(:,5)=Work3';
    WorkTable(:,6)=Work5';

    % Power
    PowerTable=WorkTable;
    PowerTable(:,3:6)=WorkTable(:,3:6)*f(a);

    % Mean power per phase, from cycles 3-5--use for statistical analyses
    Power3Table=zeros(length(stimPhase),3);
    Power3Table(:,1)=WorkTable(:,1); % frequency
    Power3Table(:,2)=WorkTable(:,2); % stimulus phase
    Power3Table(:,3)=WorkTable(:,5)*f(a); % power data
        
    cd(DirName);
    % write data to csv and append new data after each loop
    dlmwrite('power_data.csv',Power3Table,'delimiter',',','-append');
    
    % clear variables to allow plotting and analysis from next experiment
    clearvars f pData sortedphase phaseidx sorted_pData stimPhase Date Work1 Work2 Work3 Work5 WorkTable PowerTable Power3Table
    
    cd(DirName);  % return to the top-level directory
    cd(folderList(q).name); % Return to experiment folder
    
    end
  cd(DirName);  % return to the top-level directory
  
end  

%% Join pigeon data and power data spreadsheets into one table

cd(DirName); % return to the top-level directory

% read in WL_pigeon_data and power_data csvs to join them into one
pigeondat=readtable('WL_pigeon_data.csv');
power3dat=readtable('power_data.csv');
% check that both have the same number of rows
height(pigeondat)==height(power3dat);
% Rename variables in power3dat to reflect their contents
power3dat.Properties.VariableNames={'Freq_Hz' 'Phase_pc' 'Power_3workLoops_mW'};
% Combine the two tables such that power data occupies columns 9-11 (empty
% columns of pigeon data spreadsheet
pigeonpower3dat=[pigeondat(:,1:8),power3dat,pigeondat(:,12)];

%% Write new csv from joined tables

cd(DirName); % return to the top-level directory
writetable(pigeonpower3dat,'pigeon HT power by freq and phase.csv');

fclose('all');