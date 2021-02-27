%% animation for dilute polymer

clear all; clc; close all;

dlg_title='User Inputs';
prompt={'Enter the topology';...
        'Enter the com calculation mode';...
        'Enter the configuration file';...
        'Enter the com file';...
        'Enter total number of snapshots';...
        'Enter number of chains';...
        'Enter the index of the chain to be animated';...
        'Enter total number of segments';...
        'Enter number of segments in backbone (for comb)';...
        'Enter number of segments in arms (for comb)';...
        'Enter if outputing Hencky strain is desired';...
        'Enter frequency of dumping data';...
        'Enter the Weissenberg number';
        };
num_lines=1;
defaultans={'linear';'1';'./R.flow.dat';'./CoM.flow.dat';...
    '200';'120';'16';'16';'10';'10';'0';'0.01';'2'};

Inputs=inputdlg(prompt,dlg_title,num_lines,defaultans);

% specifying inputs

tplgy=Inputs{1};
com=str2num(Inputs{2});
R_file=Inputs{3};
com_file=Inputs{4};

ntime=str2num(Inputs{5}); % number of snapshots
nchain=str2num(Inputs{6}); % number of chains
ichain=str2num(Inputs{7}); % chain index to be animated
nseg=str2num(Inputs{8}); % number of segments
nseg_bb=nseg;
if (strcmp(tplgy,'comb'))
  nseg_bb=str2num(Inputs{9}); % number of segments along the backbone
  nseg_ar=str2num(Inputs{10}); % number of segments in the arm
end
strshow=str2num(Inputs{11});
frq=str2double(Inputs{12});
Wi=str2double(Inputs{13});

% ----

hstar=0.061; % the radius of the beads

A=load(R_file,' ','-ascii');
B=load(com_file,' ','-ascii');

if (strshow)
    hTxtCoords = text(75,12, sprintf('(%5.2f)',0.0), ... 
        'Color',[0.2 0.2 0.2], 'FontSize',14, 'EraseMode','normal', ... 
        'HorizontalAlignment','left', 'VerticalAlignment','top'); 
end

x = zeros(ntime,nseg+1);
y = zeros(ntime,nseg+1);
z = zeros(ntime,nseg+1);

xcm = zeros(ntime);
ycm = zeros(ntime);
zcm = zeros(ntime);

xchr = zeros(ntime);
ychr = zeros(ntime);
zchr = zeros(ntime);

time = 1:1:ntime;   % Time data

for itime=1: ntime

    osc=(itime-1)*nchain;
    osb=(itime-1)*nchain*(nseg+1) + (ichain-1)*(nseg+1);
    
    if (com)
        xcm(itime)=B(osc+ichain,1);
        ycm(itime)=B(osc+ichain,2);
        zcm(itime)=B(osc+ichain,3);
    else
        xcm(itime)=0.0;
        ycm(itime)=0.0;
        zcm(itime)=0.0;
    end       
    
    for ibead=1: nseg+1
        x(itime,ibead)=A(osb+ibead,1)+xcm(itime);
        y(itime,ibead)=A(osb+ibead,2)+ycm(itime);
        z(itime,ibead)=A(osb+ibead,3)+zcm(itime);
    end
    
end

% initialization of the figure

figure('units','normalized','outerposition',[0 0 1.0 1.0]);hold on;
set(gca,'FontName','Times New Roman','FontSize',18,'FontWeight',...
    'Bold','LineWidth',2,'TickLength',[0.015 0.015]);
set(gcf,'Renderer','opengl','DoubleBuffer','off');

h = zeros(nseg+1);
[K u0] = sphere(50);

color_b=zeros(nseg+1,3);
iarm=1;

for ibead = 1: nseg+1    
    u=hstar*sqrt(pi)*u0;
    u(:,1)=u0(:,1)+x(1,ibead);
    u(:,2)=u0(:,2)+y(1,ibead);
    u(:,3)=u0(:,3)+z(1,ibead);
    h(ibead)=trisurf(K,u(:,1),u(:,2),u(:,3));

    if (strcmp(tplgy,'comb'))
        if (ibead <= nseg_bb+1)
            color_b(ibead,:)=[50,205,50]/255;
        else
            color_b(ibead,:)=[255,0,0]/255;
        end
    elseif (strcmp(tplgy,'linear'))
        color_b(ibead,:)=[50,205,50]/255;
    end
    set(h(ibead),'FaceColor',color_b(ibead,:),'EdgeColor','none',...
        'FaceLighting','gouraud');
    axis equal;
   
end

if (strshow)
    hTxtCoords = text(75,12, sprintf('(%5.2f)',0.0), ... 
        'Color',[0.2 0.2 0.2], 'FontSize',14, 'EraseMode','normal', ... 
        'HorizontalAlignment','left', 'VerticalAlignment','top');
end
    
for iseg=1: nseg
    sx=[x(1,iseg),x(1,iseg+1)];
    sy=[y(1,iseg),y(1,iseg+1)];
    sz=[z(1,iseg),z(1,iseg+1)];
    hs(iseg) = plot3(sx,sy,sz,'linewidth',2,'color',[50,205,50]/255);
end

% CoM
hcm = line('XData',xcm(1),'YData',ycm(1),'ZData',zcm(1),'EraseMode','normal',  ... 
    'Color','c', 'Marker','.', 'MarkerSize',30);
set(hcm,'EraseMode','normal');

grid on;
alpha(1.0);
alphamap('rampdown'); 
camlight(45,45); 
lighting gouraud
view(0,90);

limitvaluex=180;
limitvaluey=25;
limitvaluez=20;
xlim([-limitvaluex,limitvaluex]);
ylim([-limitvaluey,limitvaluey]);
zlim([-limitvaluez,limitvaluez]);

xlabel('x-axis','FontName','Times New Roman',...
    'FontSize',26,'FontWeight','Bold'); % change if necessary
ylabel('y-axis','FontName','Times New Roman',...
    'FontSize',26,'FontWeight','Bold') % change if necessary
zlabel('z-axis','FontName','Times New Roman',...
    'FontSize',26,'FontWeight','Bold') % change if necessary

% prepare video output

useVideoWriter = ~verLessThan('matlab','7.11'); 
if useVideoWriter 
    vid = VideoWriter('./');
    vidObj.Quality = 100;
    vid.FrameRate = 20;
    open(vid);
else
    vid = avifile('./lin/ch84.avi','fps',30, 'quality',100,...
        'compression','None'); 
end

% Animation Loop
outputfolder=fullfile(cd,'frames');
i = 2;
while i<=ntime
    
    for ibead = 1: nseg+1
        u(:,1)=u0(:,1)+x(i,ibead);
        u(:,2)=u0(:,2)+y(i,ibead);
        u(:,3)=u0(:,3)+z(i,ibead);
        set(h(ibead),'Vertices',u);
    end
    
    for iseg=1: nseg
        sx=[x(i,iseg),x(i,iseg+1)];
        sy=[y(i,iseg),y(i,iseg+1)];
        sz=[z(i,iseg),z(i,iseg+1)];
        set(hs(iseg),'XData',sx,'YData',sy,'ZData',sz);
    end 
    
    set(hcm,'XData',xcm(i),'YData',ycm(i),'ZData',zcm(i));

    if (strshow)
        set(hTxtCoords, 'Position',[75,15], ...
        'String',sprintf('(%5.2f)',i*frq*Wi));
    end
    
    if ~ishandle(h(1,1)), break; end        %# if you close the figure
    
    % capture frame 
    if useVideoWriter 
        writeVideo(vid,getframe(gcf)); 
%         outputbase=sprintf('%3.3d.eps',i);
%         outputfull=fullfile(outputfolder,outputbase);
%         if (mod(i,50)==0)      
%             print(gcf,outputfull,'-depsc2','-r600')
%         end
    else
        vid = addframe(vid, getframe(gcf)); 
    end 
    i = i+1;
end

% close and save video output 
if useVideoWriter 
    close(vid);
else 
    vid = close(vid);
end 

% open AVI file using system default player 
%winopen('vid.avi');
