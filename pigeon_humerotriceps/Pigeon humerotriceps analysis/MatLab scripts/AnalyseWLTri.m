function [pData]=AnalyseWLTri(FileName,GR,f)

% clear
%Program to analyze workloop data, calculating the work per cycle

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
%% get phase from file name
    ix=strfind(FileName,'-');
    stix=strfind(FileName,'%');
    phase=str2double(FileName(ix(2)+1:stix-1));
    
%% Calculate work

DS=1:round(length(Lengths)/2);
US=round(length(Lengths)/2):length(Lengths)-2;

for i=1:size(Lengths,2)
    uppercurve(1,i)=trapz(flipud(Lengths(DS,i)),flipud(Forces(DS,i)));
    lowercurve(1,i)=trapz(Lengths(US,i),Forces(US,i));
    work(1,i)=uppercurve(1,i)-lowercurve(1,i);
    Velocity(:,i)=-(Lengths(2:end,i)-Lengths(1:end-1,i))*10000;  %(t_f-t_i)*sample rate  mm/s
    InstantPower(:,i)=Forces(1:end-1,i).*Velocity(:,i);    %  mW
end

%% save parameters
pData.work=work;
pData.work5=mean(work);
pData.work3=mean(work(3:5));
pData.work2=mean(work(4:5));
pData.work1=mean(work(end));
pData.Lengths=Lengths;
pData.Forces=Forces;
pData.phase=phase;
pData.PowerInstant=InstantPower;
pData.Velocity=Velocity;


