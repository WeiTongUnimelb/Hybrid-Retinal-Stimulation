
clear all;
clc;
close all;
%% load recording data: Data2Save

load("Temporal Study.mat");


%% Determine the parameters for response detection
%%These parameters can be adjusted according to the recordings.


threshold=0; %% voltage of the threshold for action potential detection.
delta=[0 20e-3];%%20ms response latency
PeakDistance=2e-3;%%peak to peak distance is at least 2ms.
Response=[];
saveData={};
mode=1;
%% 
delay=1e-3;
fs=20e3; %%sampling frequency at 20kHz
A=Data2Save{1};
Data=A{1};
FreqAmp=A{2};E=FreqAmp(:,2);
Time=Data(:,1);
Recording=Data(:,2);%membrane potential
ElecStim=Data(:,4);
OptStim=Data(:,5);
figure();plot(ElecStim);
figure();plot(OptStim);
m1=max(Recording);
m2=min(Recording);
%% finding stimulation on time

ElecOpt=ElecStim/max(ElecStim(:)+1)+OptStim/max(OptStim(:)+1)/length(unique(E));
[pks,locs]=findpeaks(diff(ElecOpt),'MinPeakHeight',0,'MinPeakDistance',3e-3*fs);
ElecOpt=ElecStim+OptStim;

[pks1,locs1]=findpeaks(diff(ElecOpt),'MinPeakHeight',0,'MinPeakDistance',3e-3*fs);

locs=[locs locs1];
locs=min(locs');
locs=locs';
locs=locs+1; %%onset of all stimulation
ElecPk=ElecStim(locs+delay*fs); 
OptPk=OptStim(locs);
stimon=locs*1/fs; %%converted to 's'


%% Cutoff recording

[pks,locs]=findpeaks(diff(diff(stimon)),'MinPeakHeight',0.5);
locs=[locs+1;length(stimon)];
cutoff=stimon(locs)+0.9;
cutoff=floor(cutoff*fs);

%% Find Response
cuton=1;
sData=[];
Response=[];
for i=1:length(cutoff) 
    data=Data(cuton:cutoff(i),:);
    cuton=cutoff(i)+1;
    [sdata,mode,Freq] = FindTemporalResponse(data,threshold,delta,PeakDistance);
    sdata(:,2:6)=sdata;
    sdata(:,1)=num2cell(Freq);
    sData=[sData;sdata];
    Freq=round(Freq,1);
    Response(i,1)=Freq;
    Response(i,2)=cell2mat(sdata(1,2));%%electrical
    Response(i,3)=cell2mat(sdata(1,3));%%optical
    R=cell2mat(sdata(:,4));
    Response(i,4)=length(find(R))/length(R); %%%percentage of response
    yyaxis left;
    ylim([m2-0.01 m1+0.01]);
    pause(1);
    
end
%% summary

[~,idx] = sort(Response(:,1));
R=Response(idx,:);
[~,idx] = sort(R(:,2));
R=R(idx,:);
[~,idx] = sort(R(:,3));
R=R(idx,:);

E=unique(R(:,2));

figure();hold on;
index=find(R(:,2)==0);%optical only
RO=R(index,:);
h1=plot(RO(:,1),RO(:,4),'-ob');

index=find(R(:,3)==0); %electrical only
RE=R(index,:);
h2=plot(RE(:,1),RE(:,4),'-or');

index=find(R(:,3) & R(:,2) ); %hybrid
RH=R(index,:);
h3=plot(RH(:,1),RH(:,4),'-ok');

legend([h1 h2 h3],{'optical only','electrical only','hybrid'})

set(gca,'FontSize',20);
xlabel('Frequency (Hz)');
ylabel('Response Efficacy (100%)');


%% final result sData: 
% Column 1 Freq, Column 2 electrical stimulation, 
% Column 3 Optical stimulation, Column 4 num of responses, Colum 5 AP
% amplitudes, Column 6 AP response latency
clearvars -except sData Data2Save;

