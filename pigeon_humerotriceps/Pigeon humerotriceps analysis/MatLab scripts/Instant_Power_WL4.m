%% program to analyze all data within a folder

clear
close all

GR=4;% GR = inverted fraction (1/4) of hole distance from end of servo motor arm  

DirName=uigetdir; % select 'Experiments' folder
cd(DirName); % change working directory to experiments folder
% read in WL_pigeon_data and power_data csvs to join them into one
pigeondat=readtable('WL_pigeon_data.csv');
% Get the mass of each HT muscle (g) sorted by Subject_ID
pigeonmassdat=unique([pigeondat(:,4),pigeondat(:,8)]);
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
    % reassign dates as corresponding subject # and assign muscle 
    % mass (which corresponds to subject # in pigeonmassdat table)
    % for normalizing power to mass
    if folderList(q).name == '2017-04-05'
        subject='Bird 1';
        HT_mass_kg=(pigeonmassdat.HuTri_mass_g(1))/1000;
    elseif folderList(q).name == '2017-04-10'
        subject='Bird 3';
        HT_mass_kg=(pigeonmassdat.HuTri_mass_g(3))/1000;
    elseif folderList(q).name == '2017-04-12'
        subject='Bird 4';
        HT_mass_kg=(pigeonmassdat.HuTri_mass_g(4))/1000;
    elseif folderList(q).name == '2017-04-19'
        subject='Bird 2';
        HT_mass_kg=(pigeonmassdat.HuTri_mass_g(2))/1000; 
    elseif folderList(q).name == '2017-06-06'
        subject='Bird 5';
        HT_mass_kg=(pigeonmassdat.HuTri_mass_g(5))/1000;
    elseif folderList(q).name == '2017-06-11'
        subject='Bird 6'; 
        HT_mass_kg=(pigeonmassdat.HuTri_mass_g(6))/1000;
    elseif folderList(q).name == '2017-06-14'
        subject='Bird 7';
        HT_mass_kg=(pigeonmassdat.HuTri_mass_g(7))/1000;
    elseif folderList(q).name == '2017-06-27'
        subject='Bird 8';
        HT_mass_kg=(pigeonmassdat.HuTri_mass_g(8))/1000;
    elseif folderList(q).name == '2017-06-28'
        subject='Bird 9'; 
        HT_mass_kg=(pigeonmassdat.HuTri_mass_g(9))/1000;
    else
        error('invalid experiment date, subject ID and muscle mass')
    end
       
    for i = 3:h % loop to move through treatment folders within a experiment folder
    cd(FolderName(i).name);
    PathName=pwd;
    A=dir(fullfile(PathName,'*.ddf')); % saves data from all .ddf files into a structure called A 
    scrsz = get(0,'ScreenSize');
    
    % name figures with date and treatment and holds them onscreen
    name = sprintf('%s',subject,' - ',string(FolderName(i).name), ' - Loop 4');
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
        % convert power to W/kg
        sorted_pData(j).MSPowerInstant=((sorted_pData(j).PowerInstant)/1000)/HT_mass_kg;
    
    %   Instantaneous power plots (4th workloop cycle)
        s(j) = subplot(3,length(A),j);
        plot(sorted_pData(j).MSPowerInstant(:,4));
        set(gca,'ytick',[],'xtick',[],'layer','bottom','box','off') 
        xlim([-1000,2000]);
        ylim([-4000,2000]);
        title(stimPhase(j)); 
        axis off
        hold on
        
        s(j) = subplot(3,length(A),j);
        plot([0,2000],[0,0],'k--');
        end 

    subplot(3,length(A),1); 
%     xticks([-3,0,3]) 
    hold on
    %y scale bar with labels at x=-50
    plot([-500; -500], [-2000; -1000], '-k',  [0; 0], [0; 0], '-k', 'LineWidth', 1)
    %x scale bar with labels at y=-2000
    plot([0; 0], [0; 0], '-k', [-500; 500], [-2000; -2000], '-k', 'LineWidth', 1)   
%     hold off
%     axis([[-3  3]    -0.2  40])
    % add the y scale bar label
    text(-2500,-1500,'1000'); 
    % add the x scale bar label
    text(-1000,-2300,'1000'); 
    % add zero line label
    text(-1000,0,'0');

    
    % clear variables to allow plotting of work loops from next experiment
    clearvars f pData sortedphase phaseidx sorted_pData stimPhase
    
    cd(DirName);  % return to the top-level directory
    cd(folderList(q).name); % Return to experiment folder
    
    end
  cd(DirName);  % return to the top-level directory

 end