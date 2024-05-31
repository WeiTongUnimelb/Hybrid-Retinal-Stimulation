function [sData,mode,Freq] = FindResponse(Data,threshold,delta,PeakDistance)
%%Find responses that evoked by the stimulation
%%Return SData: ElecAmp,OptAmp;Response Number; Response Amplitude, Response latency(from stim onset);
%%Return mode: mode=1 hybrid, 2 elec only, 3 opt only
%%Return Freq: stimulation frequency

Time=Data(:,1);
Recording=Data(:,2);
ElecStim=Data(:,4);
OptStim=Data(:,5);

sData={};

dt=Time(2)-Time(1);
fs=20e3;
dt=1/fs;
%%% below updated on 2023 Dec 21

mode=2;%%1 hybrid, 2 elec only, 0 opt only

if max(OptStim)>0 & max(ElecStim)==0 %%optical only
    mode=0;
    [stimpks,opton]=findpeaks(diff(OptStim));
    optpks=OptStim(opton+1);
    stimon=opton+1; %%opt only
elseif max(OptStim)==0 %%electrical only
    mode=2;
    [stimpks,elecon]=findpeaks(diff(ElecStim));
    [elecpks,elecoff]=findpeaks(-diff(ElecStim));%this is finding the electrical stimulation off times.for removing elecstim artefact
    elecpks=ElecStim(elecoff);
    elecoff=elecoff*dt;
    stimon=elecon+1;
else %%hybrid stimulation
    mode=1;

    [stimpks,opton]=findpeaks(diff(OptStim));
    optpks=OptStim(opton+1);
    gap=(opton(2)-opton(1))/2;
    stimon=opton+1;
    elecon=stimon;
    elecoff=elecon;
    opton=opton+1;
    for i=1:length(opton)
        k1(i)=max(1,floor(opton(i)-gap));
        k2(i)=min(length(ElecStim),floor(opton(i)+gap));
        k=find(ElecStim(k1(i):k2(i)))+k1(i);
        elecpks(i)=max(ElecStim(k1(i):k2(i)));
        if k
            elecon(i)=min(k);
            elecoff(i)=max(k);
            stimon(i)=min(elecon(i),opton(i));
        end
    end
    elecoff=elecoff*dt;
end


stimon=stimon*dt;

%%find stim freq
Freq=1/(stimon(2)-stimon(1));

%% plot stimulation and recording

figure();
hold on;
yyaxis right;
if mode
    plot(Time,ElecStim,'k'); %%plot elec stimulation for elec only or hybrid
else
    plot(Time,OptStim,'k'); %%plot optic stimulation for opt only
end
yyaxis left;
plot(Time,Recording,'r');
%% find peaks higher than threshold

fs=1/dt;
pks=[];locs=[];
if max(Recording)>=threshold
    [pks,locs]=findpeaks(Recording,fs,'MinPeakHeight',threshold,'MinPeakDistance',PeakDistance);%this is finding the action potentials distance between peaks >delta 15 ms
end
%% find response

for i=1:length(stimon)
    %%
    if mode %% for elec only and hybrid
        k1=find(locs>stimon(i) & locs<elecoff(i));%%find the response during stimulation
        k2=find(locs>elecoff(i)+delta(1) & locs<stimon(i)+delta(2));%%remove artefact (within delta(1) when stim is off) and find the response within delta(2) when stim is off
        k=[k1;k2];

    else
        k=find(locs>stimon(i) & locs<stimon(i)+delta(2));
    end

    pksi=pks(k);
    locsi=locs(k);
    if mode==1 %%hybrid
        sData{i,1}=elecpks(i);
        sData{i,2}=optpks(i);
    elseif mode==2 %%eleconly
        sData{i,1}=elecpks(i);
        sData{i,2}=0;
    else %%optonly
        sData{i,1}=0;
        sData{i,2}=optpks(i);
    end
    sData{i,3}=length(k);
    sData{i,4}=pksi';
    sData{i,5}=(locsi-stimon(i))';
    plot(locsi,pksi,'o');
end

drawnow;
if mode

    [~,idx] = sort(elecpks);
    sData = sData(idx,:);
end


end

