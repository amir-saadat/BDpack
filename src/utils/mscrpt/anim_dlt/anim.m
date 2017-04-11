%% animation for dilute polymer

% getting user input:

dlg_title='User Inputs';
prompt={'Enter the topology';...
        'Enter the com calculation mode';...
        'Enter the configuration file';...
        'Enter the com file';...
        };
num_lines=1;
defaultans={'linear';'0';'R.flow.dat';'CoM.flow.dat'};

Inputs=inputdlg(prompt,dlg_title,num_lines,defaultans);

% if com is desired 1=true, 0=false
%com=1;

A = load('./R.flow.dat',' ','-ascii');
if (com)
    B = load('./CoM.flow.dat',' ','-ascii');
end

% Position data for Beads
nchain=120;                 % number of chains
ichain=1;                   % the chain index which is desired to be animated
nseg=16;                    % number of segments
if (strcmp(tplgy,'comb'))
  nseg_bb=16;               % number of segments along the backbone
  nseg_ar=30;               % number of segments in the arm
end
ntime=400;                  % number of snapshots saved for animation
hstar=0.061;                % diameter of the beads (optional)

x = zeros(ntime,nseg+1);
y = zeros(ntime,nseg+1);
z = zeros(ntime,nseg+1);

xcm = zeros(ntime);
ycm = zeros(ntime);
zcm = zeros(ntime);

% time data
time = 1:1:ntime;

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

% h is for beads and hs is for springs:
% figure(1);
% set(gcf,'Renderer','zbuffer','DoubleBuffer','on');
%initialization of the figure
figure('units','normalized','outerposition',[0 0 1.0 1.0]);hold on;
set(gca,'FontName','Times New Roman','FontSize',18,'FontWeight',...
    'Bold','LineWidth',2,'TickLength',[0.015 0.015]);
set(gcf,'Renderer','opengl','DoubleBuffer','off');

h = zeros(nseg+1);
[K u0] = sphere(50);
color=[1 1 0; ...
       1 0 1; ...
       0 1 1; ...
       1 0 0; ...
       0 1 0; ...
       1 1 1; ...
       0 0 0; ...
       0 0 1; ...
       0.1 0.5 1; ...
       0.5 0.7 1];
% beads
for ibead = 1: nseg+1    
    u=hstar*sqrt(pi)*u0;
    u(:,1)=u0(:,1)+x(1,ibead);
    u(:,2)=u0(:,2)+y(1,ibead);
    u(:,3)=u0(:,3)+z(1,ibead);
    h(ibead)=trisurf(K,u(:,1),u(:,2),u(:,3));

%     if (ibead <= nseg_bb+1)
%         color_b(ibead,:)=[1 1 0];
%     else
%         if (ibead-nseg_bb-1-(iarm-1)*nseg_ar <= nseg_ar)
%             color_b(ibead,:)=[1 0 0];
%             if (ibead == nseg_bb+1+iarm*nseg_ar)
%                 iarm=iarm+1;
%             end,
%         end,
%     end,

    set(h(ibead),'FaceColor',[1 1 0],'EdgeColor','none',...
        'FaceLighting','gouraud');
    axis equal;  
end

% segments
% for iseg=1: nseg
%     sx=[x(1,iseg),x(1,iseg+1)];
%     sy=[y(1,iseg),y(1,iseg+1)];
%     sz=[z(1,iseg),z(1,iseg+1)];
%     hs(iseg) = plot3(sx,sy,sz,'linewidth',4,'color', [1 0 0]);
% end

% CoM
hcm = line('XData',xcm(1),'YData',ycm(1),'ZData',zcm(1),...
    'EraseMode','normal','Color','c','Marker','.','MarkerSize',30);
set(hcm,'EraseMode','normal');

% Springs initial position
x1 = zeros(ntime,nseg);x2 = zeros(ntime,nseg);
y1 = zeros(ntime,nseg);y2 = zeros(ntime,nseg);
z1 = zeros(ntime,nseg);z2 = zeros(ntime,nseg);
hs = zeros(nseg);

for itime = 1 : ntime
    for iseg = 1: nseg
        x1(itime,iseg) = x(itime,iseg);x2(itime,iseg) = x(itime,iseg+1);
        y1(itime,iseg) = y(itime,iseg);y2(itime,iseg) = y(itime,iseg+1);
        z1(itime,iseg) = z(itime,iseg);z2(itime,iseg) = z(itime,iseg+1);
    end
end

% for iseg = 1: nseg 
%     
%     dx=x2(1,iseg)-x1(1,iseg); % extent of this spring along x-axis
%     dy=y2(1,iseg)-y1(1,iseg); % extent of this spring along y-axis
%     dz=z2(1,iseg)-z1(1,iseg); % extent of this spring along z-axis
%     ds=sqrt(d3x^2+dy^2+dz^2); % needed length of this spring
% 
%     % rotate the master spring towards the needed orientation
%     p=[xs' ys' zs'*ds]*vrrotvec2mat(vrrotvec([dx dy dz],[0 0 1]));
% 
%     sx=p(:,1)+x1(1,iseg); % shift x-coordinates
%     sy=p(:,2)+y1(1,iseg); % shift y-coordinates
%     sz=p(:,3)+z1(1,iseg); % shift z-coordinates
% 
% %     hs(iseg) = plot3(sx,sy,sz,'linewidth',2,'color', 'r');   
% %     
% %     set(hs(iseg),'EraseMode','normal');
%     
% end

% grid on;
alpha(1.0);
alphamap('rampdown'); 
camlight(45,45); 
lighting gouraud
view(0,90);
linear
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
while i<=ntime
useVideoWriter = ~verLessThan('matlab','7.11'); 
if useVideoWriter 
    vid = VideoWriter('./ch1.avi');
    vidObj.Quality = 100;
    vid.FrameRate = 20;
    open(vid);
else
    vid = avifile('./ch1.avi','fps',30, 'quality',100,'compression','None'); 
end

% Animation Loop
i = 2;
    
    % beads
    for ibead = 1: nseg+1
        u(:,1)=u0(:,1)+x(i,ibead);
        u(:,2)=u0(:,2)+y(i,ibead);
        u(:,3)=u0(:,3)+z(i,ibead);
        set(h(ibead),'Vertices',u);
    end
    % segments
%     for iseg=1: nseg
%         sx=[x(i,iseg),x(i,iseg+1)];
%         sy=[y(i,iseg),y(i,iseg+1)];
%         sz=[z(i,iseg),z(i,iseg+1)];
%         set(hs(iseg),'Vertices',[sx,sy,sz]);
%     end       
    % com
    set(hcm,'XData',xcm(i),'YData',ycm(i),'ZData',zcm(i));
    
    if ~ishandle(h(1,1)), break; end        % if you close the figure
    
    % capture frame
    if useVideoWriter 
        writeVideo(vid,getframe(gcf)); 
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
winopen('vid.avi');
