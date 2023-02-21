%% program to analyze all data within a folder

clear
close all

GR=4;% GR = inverted fraction (1/4) of hole distance from end of servo motor arm  

DirName=uigetdir; % select 'Manuscript experiments' folder
cd(DirName); % change working directory to experiments folder
% Get the number of directories (excludes files) within the experiment folder
% This gives the number of experiments--number of times to run analysis loop 
all_files = dir;
all_dir = all_files([all_files(:).isdir]);
Q = numel(all_dir);

for q = 3:Q % loop to move through different experiment folders
    cd(DirName);
    folderList = dir(pwd);  % gets a list of all the experiment folders
    cd(folderList(q).name); % Sets next experiment folder as working directory
    FolderName=dir(pwd);    % gets a list of treatment folders for a given experiment
 
    % sets conditions to determine # of treatments (and loops) per experiment
    if folderList(q).name=='2017-04-05' % has two different treatments
        h=4; 
        elseif folderList(q).name=='2017-04-12' % has two different treatments
        h=4;
        else % the rest have only one treatment
        h=3;
    end
    if folderList(q).name == '2017-04-05'
        subject='Bird 1';
    elseif folderList(q).name == '2017-04-10'
        subject='Bird 3';
    elseif folderList(q).name == '2017-04-12'
        subject='Bird 4';
    elseif folderList(q).name == '2017-04-19'
        subject='Bird 2';
    elseif folderList(q).name == '2017-06-06'
        subject='Bird 5';
    elseif folderList(q).name == '2017-06-11'
        subject='Bird 6';    
    elseif folderList(q).name == '2017-06-14'
        subject='Bird 7';
    elseif folderList(q).name == '2017-06-27'
        subject='Bird 8';
    elseif folderList(q).name == '2017-06-28'
        subject='Bird 9';  
    else
        error('invalid experiment date and subject ID')
    end    
    for i = 3:h % loop to move through treatment folders within a experiment folder
    cd(FolderName(i).name);
    PathName=pwd;
    A=dir(fullfile(PathName,'*.ddf')); % saves data from all .ddf files into a structure called A 
    scrsz = get(0,'ScreenSize');
    
    % name figures with date and treatment and holds them onscreen
    name = sprintf('%s',subject,' - ', string(FolderName(i).name), ' - Loops 3-5');
    fgr=figure();
    set(fgr,'Name',name,'NumberTitle','off');
    set(fgr,'Unit','normalized','Position',[0,0,1,1]);
    plot(1:(length(i)));
    hold off   
    
        % loops through all trials for a given treatment
        for a=1:(length(A)) 
        FileName=A(a).name;
        % get trial frequency
        f(a)=extract_frequency(FileName,15,15);
        pData(a)=AnalyseWLTri(FileName,GR,f); 
        end 
        
    % sort pData by ascending phase
    [sortedphase, phaseidx] = sort([pData.phase]);
    sorted_pData = pData(phaseidx);

        for j=1:(length(A))
        % get phases from file name
        stimPhase(j)=sorted_pData(j).phase;
    
    %   force plots (workloops 3-5)
        s(j) = subplot(3,length(A),j);
        plot(sorted_pData(j).Lengths(:,3:5),sorted_pData(j).Forces(:,3:5));
        set(gca,'ytick',[],'xtick',[],'layer','bottom','box','off') 
        % y min just below 0 prevents loop bottom cutoffs
        % y max must be 40 -- 10.1 Hz trials and for consistency
        ylim([-0.2,40]); 
        xlim([-3,3]); 
        title(stimPhase(j)); 
        axis off 
        end 

    subplot(3,length(A),1); 
    xticks([-3,0,3]) 
    hold on
    plot([-3; -3], [0; 10], '-k',  [0; 0], [0; 0], '-k', 'LineWidth', 2)    %y scale bar with labels at x=-3
    plot([0; 0], [0; 0], '-k', [-3; 3], [-0.2; -0.2], '-k', 'LineWidth', 1) %x scale bar with labels at y=-0.2
    hold off
    axis([[-3  3]    -0.2  40])
    text(-5,0,'0'); 
    text(-6.5,10,'10'); 
    text(-4,-2,'-3'); 
    text(2,-2,'3'); 
    
    % clear variables to allow plotting of work loops from next experiment
    clearvars f pData sortedphase phaseidx sorted_pData stimPhase
    
    cd(DirName);  % return to the top-level directory
    cd(folderList(q).name); % Return to experiment folder
    
    end
  cd(DirName);  % return to the top-level directory

 end