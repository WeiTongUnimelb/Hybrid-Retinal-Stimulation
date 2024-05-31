% clc;clearvars -except S;
clc;close all;
clear all;

%% Select folder where files are located
disp('Select folder where folders are located')
selpath = uigetdir;
file=dir([selpath '/*.mat']);
numfile=length(file);
category={numfile};
cellsactive=[];
distance_all=[];
g = [];
Stim_CR = [];
Stim_notCR = [];

    %Calculate distance from electrode

    distance_elec={numfile};
%4/05/23  elec_x=201; and elec_y=505; 
%9/05/23  elec_x=364; and elec_y=432;
%29/05/23 elec_x=152; and elec_y=139; 
%elec co-ordinates - needs to be flipped
    elec_x=152*1.24;
    elec_y=(256-139)*1.24;
load(strcat(selpath,'/', 'chrimsonrcells.mat'))
chrimsonrcells=chrimsonrcells*1.24;
area=zeros(1,numfile);
map=[ 0.1 0.1 0.1; 1 0.2 0; 1 0.6 0; 1 1 0];%define colourmap
for fileno=1:1:numfile-1
      
    f=file(fileno).name;
    load(strcat(selpath,'/', file(fileno).name))
      
      unstimulated_cell_location=1.24*unstimulated_cell_location;
      stimulated_cell_location=1.24*stimulated_cell_location;
      distance_combined=stimulated_cell_location;
      chrimsonrcellsstimmed=intersect(stimulated_cell_location,chrimsonrcells,'rows');

    figure()
    %imshow(strcat(selpath,'/MosaicJ.jpg'))
    hold on
    %256 pixels = 317.952 microns
    %1 pixel =1.242 microns
    %30 x 30 µm square will be 24.15 pixels 
    plot(unstimulated_cell_location(:,1),unstimulated_cell_location(:,2),'o','MarkerEdgeColor',[0.8 0.8 0.8],'MarkerFaceColor','white','MarkerSize',8);

    plot(stimulated_cell_location(:,1),stimulated_cell_location(:,2),'o','MarkerEdgeColor',[0.8 0.8 0.8],'MarkerFaceColor','red','MarkerSize',8);
    plot(chrimsonrcells(:,1),chrimsonrcells(:,2),'x','MarkerEdgeColor',[0.8 0.8 0.8],'MarkerFaceColor','blue','MarkerSize',8);
    plot(chrimsonrcellsstimmed(:,1),chrimsonrcellsstimmed(:,2),'o','MarkerEdgeColor',[0.8 0.8 0.8],'MarkerFaceColor','blue','MarkerSize',8);
    axis equal
    title(file(fileno).name)
    %find area stimulated 
    figure()
    if stimulated_cell_location~=0
    %plot the stimulated cell location
    hist3(stimulated_cell_location,'Ctrs',{12:24:336 12:24:336},'CdataMode','auto')
    hold on
    plot(elec_x,elec_y,'o','MarkerEdgeColor',[0.8 0.8 0.8],'MarkerFaceColor','white','MarkerSize',8)
    axis equal
    
    caxis manual
    caxis([0 3]);
    cmap=colormap(map);
    colorbar
    view(2) 
    
    [N,edges]=hist3(stimulated_cell_location,'Ctrs',{12:24:336 12:24:336},'CdataMode','auto');
    area(fileno)=nnz(N)*30e-6*30e-6*1e6;%(mm^2)
    title(file(fileno).name)
    
    %
    distance_elec_x=abs(distance_combined(:,1)-elec_x);
    distance_elec_y=abs(distance_combined(:,2)-elec_y);
    distance=sqrt(distance_elec_x.^2+distance_elec_y.^2);
    end
    %distance={cat(1,distance_elec{:})};
    %distance=cell2mat(distance);
    TF = isempty(distance_combined);
    if TF == 1
        distance=0;  
    end    
    distance_all=[distance_all;distance];
    g = [g; fileno*ones(size(distance))];
    cellsactive=[cellsactive;length(distance_combined)];
    OpticalStim=regexp(file(fileno).name,'\d*','Match','Once');
    ElectricalStim=regexp(file(fileno).name,'\d*','Match');
    ElectricalStim=ElectricalStim(2);
    category{fileno}=strcat(OpticalStim,'/',ElectricalStim);
    Stim_CR=[Stim_CR;length(chrimsonrcellsstimmed)];
    Stim_notCR=[Stim_notCR;(length(distance_combined)-length(chrimsonrcellsstimmed))];
    %show mosaic image
    
end


figure()
boxplot(distance_all,g)
xticklabels(cat(1,category{:}))
ylabel('Distance (pixels')
figure()
bar(cellsactive)
xticklabels(cat(1,category{:}))
ylabel('Number of stimulated cells')
figure()
bar(area)
xticklabels(cat(1,category{:}))
ylabel('Area Active')
% figure()
% bar(cellsactive./area')
% xticklabels(cat(1,category{:}))
% ylabel('Cells/area')
figure()
bar(Stim_CR)
xticklabels(cat(1,category{:}))
ylabel('Cells expressing chrimson R')
figure()
bar(Stim_notCR)
xticklabels(cat(1,category{:}))
ylabel('Cells not expressing chrimson R')
figure()
bar((Stim_CR./cellsactive))
xticklabels(cat(1,category{:}))
ylabel('Cells with chrimsonR/ total active cells')

