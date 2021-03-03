%% Reading the files for the animation

clear all; clc; close all;

A = importdata('Rb_elong.dat',' ');
B = importdata('CoM_elong.dat',' ');
C = importdata('BoxConfig_elong.dat',' ');

%% Selecting the color map to be used for the facecolor of the beads

myMap=get(groot,'defaultAxesColorOrder');

%% The main animation code

% ViewBox= 30; % Viewing Box
L=[30.2 30.2 30.2]; % Simulation Box

% Set the flow type here: Equil, PSF, PEF
FlowType='PEF';

% --- Changing the outer box specification details ---
% az=0; % azimuth view degrees
% el=90; % elevation view degrees
% view(az,el)
% xmin=0; xmax=ViewBox; ymin=0; ymax=ViewBox; zmin=0; zmax=ViewBox;
% axis([xmin xmax ymin ymax zmin zmax])
% axis on % may be off
% grid off % may be off
% box on % may be off
% camlight;


% Position data for Beads
nrun=10; % Total number of runs (independent box simulations)
irun=1; % The number of run currently visualized
nchain=20; % Number of chain
nseg=9; % Number of segments
nbead=nseg+1; % Number of beads
ntotseg=nchain*nseg; % Total number of segments in the box
ntotbead=nbead*nchain; % Total number of beads in the box
ntime=100; % Total number of dumps in each run (ndmp in simulations)
hstar=0.05; % h*

x = zeros(ntime,nchain,nbead);
y = zeros(ntime,nchain,nbead);
z = zeros(ntime,nchain,nbead);

xcm = zeros(ntime,nchain);
ycm = zeros(ntime,nchain);
zcm = zeros(ntime,nchain);

delrxL = zeros(ntime);

L1=zeros(ntime,2);
L2=zeros(ntime,2);

time = 1:1:ntime;   % Time data

osb=(irun-1)*ntotbead*ntime;
osc=(irun-1)*nchain*ntime;

for itime=1: ntime  
    for ichain=1: nchain
        for ibead = 1: nbead
            x(itime,ichain,ibead)=A(osb+ntotbead*(itime-1)+(ichain-1)*nbead+ibead,1);
            y(itime,ichain,ibead)=A(osb+ntotbead*(itime-1)+(ichain-1)*nbead+ibead,2);
            z(itime,ichain,ibead)=A(osb+ntotbead*(itime-1)+(ichain-1)*nbead+ibead,3);            
        end
        xcm(itime,ichain)=B(osc+nchain*(itime-1)+ichain,1); 
        ycm(itime,ichain)=B(osc+nchain*(itime-1)+ichain,2);
        zcm(itime,ichain)=B(osc+nchain*(itime-1)+ichain,3);
    end
    switch FlowType
        case 'PSF'
            delrxL(itime)=C(itime);
        case 'PEF'
            L1(itime,1)=C(itime,1);L1(itime,2)=C(itime,2);
            L2(itime,1)=C(itime,3);L2(itime,2)=C(itime,4);
    end
end

% h is for beads and hs is for springs:

% Initialization of the figure
figure('units','normalized','outerposition',[0 0 1.0 1.0]);hold on;
set(gca,'FontName','Times New Roman','FontSize',18,'FontWeight',...
    'Bold','LineWidth',4,'TickLength',[0.015 0.015]);
set(gcf,'Renderer','opengl','DoubleBuffer','off');

h=zeros(nchain,nbead);
hcm=zeros(nchain);
[K u0]=sphere_wp(30);
% Beads initial position
for ichain=1: nchain
    
    % --- This part manually assign colors for the beads. change as needed ---
    color=[1 1 0; ...
           1 0 1; ...
           0 1 1; ...
           1 0 0; ...
           0 1 0; ...
           0 0 1; ...
           0.1 0.5 1; ...
           0.2 0.4 0.3; ...
           0.1 0.2 0.7; ...
           0.5 0.1 0.1; ...
           0.1 0.8 0.3; ...
           0.4 0.4 0.4; ...
           0.6 0.5 1; ...
           1 0.3 0.1; ...
           0.2 1 0.6; ...
           0.7 0.9 0.9; ...
           0.5 1 1; ...
           0.2 0.2 0.9; ...
           0.9 0.1 0.1; ...
           0.3 0.3 0.3; ...
           0.9 0.9 0.9; ...
           0.7 0.1 0.9; ...
           0.6 0.2 0.2; ...
           0.5 0.5 0.7; ...
           0.8 0.1 0.9; ...
           0.1 0.1 0.1; ...
           0.6 0.5 0.9; ...
           1 1 0.9; ...
           0.2 1 0.2; ...
           0.5 0.7 1;...
           0.1 0.3 0.3; ...
           0.2 0.1 0.5; ...
           0.1 0.1 0.8; ...
           0.5 0.2 0.7; ...
           0.3 0.3 0.5; ...
           0.1 0.2 0.2; ...
           0.6 0.3 0.8; ...
           0.9 0.2 0.2; ...
           0.1 1 0.9; ...
           0.5 0.2 0.1];
    for ibead=1: nbead
        u=hstar*sqrt(pi)*u0;
        u(:,1)=u0(:,1)+x(1,ichain,ibead);
        u(:,2)=u0(:,2)+y(1,ichain,ibead);
        u(:,3)=u0(:,3)+z(1,ichain,ibead);
        % Handle for all beads of all chains
        h(ichain,ibead)=trisurf(K,u(:,1),u(:,2),u(:,3));
        set(h(ichain,ibead),'FaceColor',color(ichain,:),'EdgeColor','none','FaceLighting','gouraud');
        axis equal;
        % --- Change the light position as needed ---
        %light('Position',[1 -1 0],'Style','infinite');
        %light('Position',[10 10 10],'Style','infinite'); 
    end
    for iseg=1: nseg
        sx=[x(1,ichain,iseg),x(1,ichain,iseg+1)];
        sy=[y(1,ichain,iseg),y(1,ichain,iseg+1)];
        sz=[z(1,ichain,iseg),z(1,ichain,iseg+1)];
        % Handle for all segments of all chains
        hs(ichain,iseg) = plot3(sx,sy,sz,'linewidth',4,'color', color(ichain,:));
    end
    % CoM
    hcm(ichain) = line('XData',xcm(1,ichain),'YData',ycm(1,ichain),'ZData',zcm(1,ichain),...
        'EraseMode','normal','Color',color(ichain,:), 'Marker','*', 'MarkerSize',10);
    set(hcm(ichain),'EraseMode','normal');
end

% Deformed Box

% h(1) = axes('Position',[0.2 0.2 0.6 0.6]); 
% vert = [-L(1)/2 -L(2)/2 -L(3)/2; L(1)/2 -L(2)/2 -L(3)/2; -L(1)/2 L(2)/2 -L(3)/2; L(1)/2 L(2)/2 -L(3)/2 ; ... 
%         -L(1)/2 -L(2)/2 L(3)/2; L(1)/2 -L(2)/2 L(3)/2; -L(1)/2 L(2)/2 L(3)/2; L(1)/2 L(2)/2 L(3)/2]; 
vert = [0 0 0; L(1) 0 0; 0 L(2) 0; L(1) L(2) 0 ; ... 
        0 0 L(3); L(1) 0 L(3); 0 L(2) L(3); L(1) L(2) L(3)]; 
fac = [1 3 4 2; ... 
    1 3 7 5; ... 
    1 2 6 5; ... 
    2 4 8 6; ... 
    3 4 8 7; ... 
    5 6 8 7];
hr=patch('Faces',fac,'Vertices',vert,'FaceColor',myMap(1,:),'EdgeColor',myMap(1,:),'LineStyle','-','LineWidth',4);  % patch function 
% set(hr,'FaceLighting','phong','EdgeLighting','phong');
set(hr,'EraseMode','normal','FaceAlpha',0.1);
% light('Position',[1 3 2]);
% light('Position',[-3 -1 3]);
% material shiny;
% alpha(1.0);
alphamap('rampdown'); 
camlight(25,45);
camlight('headlight');
lighting gouraud 
view(25,45);


% --- Setting for your outer box names (change if needed) ---
% xlabel('x','FontName','Times New Roman',...
%     'FontSize',26,'FontWeight','Bold'); % change if necessary
% xlabh = get(gca,'xlabel');
% set(xlabh,'Position',[0 -4 0]);
% ylabel('y','FontName','Times New Roman',...
%     'FontSize',26,'FontWeight','Bold') % change if necessary
% ylabh = get(gca,'ylabel');
% set(ylabh,'Position',[66 20 0]);
% zlabel('z','FontName','Times New Roman',...
%     'FontSize',26,'FontWeight','Bold') % change if necessary

set(gca,'xtick',[],'ytick',[],'ztick',[])
box on
% set(gca,'visible','off')

% --- Setting for your outer box (change if needed) ---
limitvaluex=60;
limitvaluey=55;
limitvaluez=35;
xlim([-40,60]);
ylim([0,40]);
zlim([0,limitvaluez]);

% prepare video output 
useVideoWriter = ~verLessThan('matlab','7.11'); 
if useVideoWriter 
    vid = VideoWriter('vid1.avi');
    vidObj.Quality = 100; 
    vid.FrameRate = 20; 
    open(vid); 
else
    vid = avifile('vid1.avi','fps',20,'quality',100,'compression','None'); 
end

% Animation Loop
i = 2;
while i<=ntime
    
    for ichain=1: nchain
        for ibead=1: nbead
            u(:,1)=u0(:,1)+x(i,ichain,ibead);
            u(:,2)=u0(:,2)+y(i,ichain,ibead);
            u(:,3)=u0(:,3)+z(i,ichain,ibead);
            set(h(ichain,ibead),'Vertices',u);
            %set(h(ibead),'XData',x(i,ibead),'YData',y(i,ibead),'ZData',z(i,ibead));
            %drawnow;
%             pause(1/50);
        end
        for iseg=1: nseg
            sx=[x(i,ichain,iseg),x(i,ichain,iseg+1)];
            sy=[y(i,ichain,iseg),y(i,ichain,iseg+1)];
            sz=[z(i,ichain,iseg),z(i,ichain,iseg+1)];
            set(hs(ichain,iseg),'XData',sx,'YData',sy,'ZData',sz);
            %drawnow;
            %pause(1/50);
        end        
        set(hcm(ichain),'XData',xcm(i,ichain),'YData',ycm(i,ichain),'ZData',zcm(i,ichain));
        %drawnow;pause(1/50);
        %set(hTxtcom,'XData',xc(i),'YData',yc(i),'ZData',zc(i));
    end
    switch FlowType
        case 'PSF'
            vert = [0 0 0; L(1) 0 0; delrxL(i) L(2) 0; delrxL(i)+L(1) L(2) 0; 0 0 L(3); L(1) 0 L(3); ...
                delrxL(i) L(2) L(3); delrxL(i)+L(1) L(2) L(3)];
        case 'PEF'
            vert = [0 0 0; L1(i,1) L1(i,2) 0; L2(i,1) L2(i,2) 0; L2(i,1)+L1(i,1) L2(i,2)+L1(i,2) 0 ; ... 
              0 0 L(3); L1(i,1) L1(i,2) L(3); L2(i,1) L2(i,2) L(3); L2(i,1)+L1(i,1) L2(i,2)+L1(i,2) L(3)];
    end
    set(hr,'Vertices',vert);
    
     % if you close the figure
    if ~ishandle(h(1,1)), break; end       
    
    %# capture frame 
    if useVideoWriter 
        writeVideo(vid,getframe(gcf)); 
    else 
        vid = addframe(vid, getframe(gcf)); 
    end 
    i = i+1;
    
    % Printing the last frame in the home directory
    if (i==ntime)
        hold off;
        print('~/semi.png','-dpng','-r600');
    end
    
    
end

%# close and save video output 
if useVideoWriter 
    close(vid); 
else 
    vid = close(vid); 
end 

%# open AVI file using system default player 
winopen('vid1.avi');
