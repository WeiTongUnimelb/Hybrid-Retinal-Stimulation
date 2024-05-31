

clc;close all;clear all;
%%

load ("Threshold Study.mat");

%% define parameters: These parameters need to be adjusted according to every recording

threshold=0;
delta=[5e-4 20e-3];%%+0.5ms is electrical artefact latency or set it to -0.5ms to remove the function; 15ms is response latency, count from stim onset
PeakDistance=2e-3;%%peak to peak distance is at least 5ms.
pulses=10;
Num=length(Data2Save);
TrialNum=1:Num;

%% find Response

saveData=[];
for n=1:length(TrialNum)
    i=TrialNum(n);
    A=Data2Save{i};
    Data=A{1};
    [sData,mode,Freq]=FindThresholdResponse(Data,threshold,delta,PeakDistance);
    saveData=[saveData;sData];
end

R=saveData(:,1:3);
R=cell2mat(R);
elecpks=R(:,1);
optpks=R(:,2);
[~,idx] = sort(optpks);
saveData = saveData(idx,:);

%% find electrical and optical stimulation conditions

E=[];O=[]; Response=[]; %%optical and electrical stimulation amplitudes

E=unique(elecpks);O=unique(optpks);
R=R(idx,:);

for i=1:length(E)
    for j=1:length(O)
        k=j+(i-1)*length(O);
        Response(k,1)=E(i);
        Response(k,2)=O(j);
        index=find(R(:,1)==E(i) & R(:,2)==O(j));
        index=find(R(index,3));
        Response(k,3)=length(index);
    end
end

%% find electrical thresholds

if length(E) > 1

    ft = fittype( [num2str(pulses) '/(1+exp(-b*(x-c)))'], 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';

    b=E(floor(length(E)/2));
    opts.StartPoint = [pulses/b,b];

    TE=zeros(length(O),1);

    figure();
    j=0;%%for plotting plots
    for i=1:length(O)

        if j<13
        else
            j=0;  figure();
        end
        j=j+1;

        if length(O)>1
            subplot(4,3,j);

        end

        hold on;

        index=find(Response(:,2)==O(i));
        data=Response(index,:);
        x=data(:,1);y=data(:,3); %%x:current,y:response
        plot(x,y,'o');ylim([0 pulses]);

        if max(y)>7 & min(y)<5  %%change max and min number if needed
            [fitresult, gof] = fit( x, y, ft, opts );
            plot(fitresult);
            TE(i)=fitresult.c;
            if TE(i)<0
                TE(i)=0;
            end
        end
        title(['Laser Power: ' num2str(O(i))]);

    end

    if length(O)==1
        fprintf('Electrical Threshold: %f pA\n',TE);
    else
        figure();plot(O,TE,'o');
        xlabel('Optical Power(mV)');
        ylabel('Electrical Current(pA)');
    end
end

%% find optical thresholds

if length(O) > 1

    ft = fittype( [num2str(pulses) '/(1+exp(-b*(x-c)))'], 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';

    b=O(floor(length(O)/2));
    opts.StartPoint = [pulses/b,b];

    TO=zeros(length(E),1);

    figure();
    j=0;
    for i=1:length(E)

        if j<13
        else
            j=0;  figure();
        end
        j=j+1;

        if length(E)>1
            subplot(4,3,j);
        end

        hold on;

        index=find(Response(:,1)==E(i));
        data=Response(index,:);
        x=data(:,2);y=data(:,3);
        plot(x,y,'o');ylim([0 pulses]);

        if max(y)>7 & min(y)<5 %%change max and min number if needed
            [fitresult, gof] = fit( x, y, ft, opts );
            plot(fitresult);
            TO(i)=fitresult.c;
            if TO(i)<0
                TO(i)=0;
            end
        end
        title(['Electrical Current: ' num2str(E(i))]);

    end


    if length(E)==1
        fprintf('Optical Threshold: %f \n',TO);
    else
        figure();plot(E,TO,'o');
        ylabel('Optical Power(mV)');
        xlabel('Electrical Current(pA)');
    end
end


%% Color Response

cm=parula;
figure();hold on;
for i=1:length(Response)
    x=Response(i,1);
    y=Response(i,2);
    z=Response(i,3);
    c=ceil((z+1)/11*256);
    c=cm(c,:);
    plot(x,y,'o','MarkerFaceColor',c,'MarkerSize',20,'MarkerEdgeColor',[0.8 0.8 0.8]);
end

set(gca,'FontSize',20)
xlabel('Current pA');
ylabel('Laser Power mV')
xlabel('Electric Current pA');
box on
colorbar
colorbar('Ticks',[0 0.2 0.4 0.6 0.8 1],...
    'TickLabels',{'0','2','4','6','8','10'})

%% final result saveData:
% Column 1 Freq, Column 2 electrical stimulation,
% Column 3 Optical stimulation, Column 4 num of responses, Colum 5 AP
% amplitudes, Column 6 AP response latency

saveData(:,2:6)=saveData;
saveData(:,1)=num2cell(Freq);
clearvars -except saveData Data2Save;
