%% 

% clc;clearvars -except S;
clc;close all;
clear all;

%% Select folder where files are located
disp('Select folder where folders are located')
selpath = uigetdir;


%% for each folder
findfolders=dir(selpath);
dirFlags=[findfolders.isdir];
subfolders =findfolders(dirFlags);
for sf = 3:length(subfolders)
    close all;
    currDir=subfolders(sf).name;%get name of the subfolder
    filepath=strcat(selpath,'/',currDir);
    file=dir([filepath '/*.oif']);
    filezip=dir([selpath '/*.zip']);
    numfile=length(file);
    stimulated_cell_location={};
    unstimulated_cell_location={};
    all_cell_locations={};
    spikenum_all={};


    for fileno=1:1:numfile
        %get file name
        f=file(fileno).name;
        fname=strcat(filepath,'/', f);


        %get associated ROI file
        froi=filezip(fileno).name;
        ROIName=strcat(selpath,'/',froi);     

        %current=100;

        %SaveFolder='C:/Users/ebrunton/Documents/datas/calcium imaging/Results/';
        dt=0.12933;%0.21820;
        fs=1/dt;
        dx=1;
        % dx=317.952/256;
        %  dx=dx*2;
        Filter=[2 1 -1 -2];
        window=4;
        %%%%background time duration (in second) 
        BackStart=1;
        BackEnd=2;
        %%%%Stimulation Interpulse Duration (in second)
        Interpulse=2.1;
        %%%%Stimulation Start Time (in second)
        StartTime=0.5;
        %%%%%Stimulation Repeats
        Repeats=10;
        %%%Stimulation Phase
        Phase='0.033ms';
        phase=0.033;
        W=256;H=256;
        ft = fittype( '10/(1+exp(-a*(x-b)))+0', 'independent', 'x', 'dependent', 'y' );    
        sROI=ReadImageJROI(ROIName);
        cells=length(sROI);
        [X Y]=meshgrid(1:W,1:H);
        CELLS={};corx=[];cory=[];
         for i=1:cells
                ellipse=[];
                v=sROI{i}.vnRectBounds;
                Y0=v(1)+v(3);
                X0=v(2)+v(4);
                X0=X0/2;Y0=Y0/2;
                l=abs(v(1)-Y0);
                w=abs(v(2)-X0);
                corx=[corx;X0];cory=[cory;Y0];
                ellipse=((X-X0)/l).^2+((Y-Y0)/w).^2<=1;
                CELLS=[CELLS;find(ellipse)];

         end      
        %for step=1:length(current) %looking at the image calculating the intensity
        %of the image
        %  fname=sprintf('%s%d %d.oif',Folder,Image(step),current(step))

           A = bfopen(fname); 
           A=A{1};A=A(:,1);
           frames=length(A);

            I=A{1};%figure();hold on;imshow(I,[min(I(:)),max(I(:))]);ang=0:0.01:2*pi;
        %     for i=1:cells
        %         v=sROI{i}.vnRectBounds;
        %         Y0=v(1)+v(3);
        %         X0=v(2)+v(4);
        %         X0=X0/2;Y0=Y0/2;
        %         l=abs(v(1)-Y0);
        %         w=abs(v(2)-X0);
        %         ellipse=((X-X0)/l).^2+((Y-Y0)/w).^2<=1;
        %         x=X0+w*cos(ang);y=Y0+l*sin(ang);
        % %         plot(x,y,'r');
        %     end  
        F=[];
        for k = 1:frames
                I=A{k};
                %F_ALL(k)=mean(I,'all');% mean intensity over whole image
            for i=1:cells
                F(k,i)=mean(I(CELLS{i}));
                
            end
            F_ALL(k)=mean(F(k,:));% mean intensity over ALL CELLS
        end

        figure();
        plot(F_ALL)
        figure();
        deltaF=(F-mean(F))./mean(F);
        plot(deltaF)
        %prompt user to select start time
        prompt = "What is the start time of stimulus? (integer) ";
        stimulus_start = input(prompt);

        T=dt:dt:dt*frames; 
        filtered=zeros(frames,cells);
        bs=floor(BackStart/dt);
        be=floor(BackEnd/dt);
        SNR=zeros(1,cells);
        mean_signal=zeros(1,cells);
        stimulated=zeros(1,cells);
        
        for i=1:cells
                    pk=[];
                    loc=[];
                    f=F(:,i);
                    f=(f-mean(f))/mean(f); 
                    %f=highpass(f,0.3,fs);
                    F1=conv(f,Filter);%removes DC component
                    F1=F1(2:end-2);
                    %calculate welch PSDF estimate
                    [pxx,f2]=pwelch(f,[],[],[],fs);%change to f to get raw signal
                    %calculate SNR at 0.5 Hz/ 2 Hz
                    SNR(i)=pxx(round((1/2.1)/(max(f2)/length(f2))))/mean(pxx(round(2/(max(f2)/length(f2))):end));%divide time by the frequency interval???
                    mean_signal(i)=mean(f);
                    time=1:1:length(F1);
                    time=time*dt;
                    F2=F1;
                    std_sig=std((F1(30/dt:end)));
                   threshold=mean(F1(30/dt:end))+1.645*std_sig;%1.4*max(F2(floor(32/dt):end));
                   timebtwpeaks=floor(2/dt); %2.1 seconds between stimuli
                   peakwidth=floor(0.6/dt);
                    if  SNR(i)<5 ||stimulus_start == 0%
                      spikenum=0;
                    else
                    [pk, loc]=findpeaks(F1,'MinPeakHeight',threshold,'MinPeakDistance',timebtwpeaks);%,'MinPeakWidth',peakwidth);
                     
                    %filtering data
                    if isempty(loc)
                    
                    else
                    firstspike=loc(1);
                    indices = find((loc<stimulus_start)|loc>(stimulus_start+round(21/dt)));
                    loc(indices)=[];
                    pk(indices)=[];
                    end
                    

        %             spikeloc=[];
        %             spike=[];
                     spikenum=length(pk);
        % %            amp=0;
        %              for k=1:Repeats
        %                  id=find(loc>(Interpulse*(k-1)+StartTime)/dt+1&loc<=(Interpulse*(k-1)+StartTime)/dt+window);
        %                  if(~isempty(id))
        %                      spikeloc(1,end+1:end+length(id)) = loc(id)';
        %                      spike(1,end+1:end+length(id))=pk(id)';
        %                      spikenum=spikenum+1;   
        %              end
                    end  
                   if i<17
                   figure(2+fileno*8);subplot(4,4,i);plot(time,F2,loc*dt,pk,'o');
                   figure(2+fileno*8);subplot(4,4,i);plot(f2,10*log10(pxx));%plot welch power spectrum
                   title(SNR(i));
                   elseif i <33
                   figure(3+fileno*8);subplot(4,4,i-16);plot(time,F2,loc*dt,pk,'o'); 
                   %figure(3+fileno*8);subplot(4,4,i-16);plot(f2,10*log10(pxx));%plot welch power spectrum
                   title(SNR(i));
                   elseif i<49
                   figure(4+fileno*8);subplot(4,4,i-32);plot(time,F2,loc*dt,pk,'o');
                   %figure(4+fileno*8);subplot(4,4,i-32);plot(f2,10*log10(pxx));
                   title(SNR(i));
                   elseif i<65        
                   figure(5+fileno*8);subplot(4,4,i-48);plot(time,F2,loc*dt,pk,'o');
                   %figure(5+fileno*8);subplot(4,4,i-48);plot(f2,10*log10(pxx));
                   title(SNR(i));
                   elseif i<81
                   figure(6+fileno*8);subplot(4,4,i-64);plot(time,F2,loc*dt,pk,'o'); 
                   %figure(6+fileno*8);subplot(4,4,i-64);plot(f2,10*log10(pxx));
                   title(SNR(i));
                   elseif i<97
                   figure(7+fileno*8);subplot(4,4,i-80);plot(time,F2,loc*dt,pk,'o'); 
                   %figure(7+fileno*8);subplot(4,4,i-80);plot(f2,10*log10(pxx));
                   title(SNR(i));
                   else
                       
                    end

                    stimulated(i)=spikenum;         

        end
             %  here we find position of stimulated and unstimulated cells

            %figure(1);hold on;
            %index=find(SNR(:)<=10);%not stimulated
            index=find(stimulated(:)<5);
           % plot(corx(index),cory(index),'o','MarkerEdgeColor',[0.8 0.8 0.8],'MarkerFaceColor','white','MarkerSize',6);
            %index2=find(SNR(:)>10);%stimulated
            index2=find(stimulated(:)>=5);
           % plot(corx(index2),cory(index2),'o','MarkerEdgeColor',[0.8 0.8 0.8],'MarkerFaceColor','red','MarkerSize',8);     
            %axis image;
            %set(gca,'Ydir','reverse');
            %title(['Phase: ' Phase ' Current: ' '\muA']);
            %set(gca,'FontSize',20);
        %save stimulated and unstimulated cells
    stimulated_cell_location{fileno}=[corx(index2),cory(index2)];
    unstimulated_cell_location{fileno}=[corx(index),cory(index)];
    all_cell_locations{fileno}=[corx,cory];
    SNR(SNR > 20) = 20;
    SNR_allcells{fileno}=SNR;
    SNR_stimmed{fileno}=SNR(index2);
    SNR_unstimmed{fileno}=SNR(index);
    spikenum_all{fileno}=stimulated;
    %end files for loop
    end


    %show mosaic image
    figure()
    imshow(strcat(selpath,'/Loc0.jpg'))
    hold on
    
    %plot stimulated cells %plot adjusted cell positions
%    
%     for count=1:1:numfile
%         plot(unstimulated_cell_location{1,count}(:,1),unstimulated_cell_location{1,count}(:,2),'o','MarkerEdgeColor',[0.8 0.8 0.8],'MarkerFaceColor','white','MarkerSize',6);
%         plot(stimulated_cell_location{1,count}(:,1),stimulated_cell_location{1,count}(:,2),'o','MarkerEdgeColor',[0.8 0.8 0.8],'MarkerFaceColor','red','MarkerSize',6);
%         hold on
%  
%     end
    %remove duplicates
     unstimulated_cell_location={cat(1,unstimulated_cell_location{:})};
     unstimulated_cell_location=cell2mat(unstimulated_cell_location);
     %[unstimulated_cell_location,UA, UB]=uniquetol(unstimulated_cell_location,0.01,'ByRows',true);
     stimulated_cell_location={cat(1,stimulated_cell_location{:})};
     stimulated_cell_location=cell2mat(stimulated_cell_location);
     %[stimulated_cell_location, SA, SB]=uniquetol(stimulated_cell_location,0.01,'ByRows',true);
     
    plot(unstimulated_cell_location(:,1),unstimulated_cell_location(:,2),'o','MarkerEdgeColor',[0.8 0.8 0.8],'MarkerFaceColor','white','MarkerSize',6);
    plot(stimulated_cell_location(:,1),stimulated_cell_location(:,2),'o','MarkerEdgeColor',[0.8 0.8 0.8],'MarkerFaceColor','red','MarkerSize',6);
%plot(stimulated_cell_location(:,1),stimulated_cell_location(:,2),'o','MarkerEdgeColor',[0.8 0.8 0.8],'MarkerFaceColor','red','MarkerSize',6);
%     

    
    
    %savefigure
    saveas(gcf,strcat(filepath,'stim_notstim'),'jpg')
    
    
    
    
    
    figure()
    
    hold on
    
    
        all_cell_locations={cat(1,all_cell_locations{:})};
        all_cell_locations=cell2mat(all_cell_locations);    
        %[all_cell_locations,IA,IC]=uniquetol(all_cell_locations,0.01,'ByRows',true);
        spikenum_all={cat(2,spikenum_all{:})};
        spikenum_all=cell2mat(spikenum_all);
        %spikenum_all=spikenum_all(IA);
        scatter(all_cell_locations(:,1),all_cell_locations(:,2),10,spikenum_all,'filled');
        hold on
 
    colormap(jet)
    colorbar

    
    %savefigure
    saveas(gcf,strcat(filepath,'_spikenum'),'jpg')
    %savecelldata
    save(filepath,'unstimulated_cell_location','stimulated_cell_location','all_cell_locations','SNR_allcells','spikenum_all')
end



