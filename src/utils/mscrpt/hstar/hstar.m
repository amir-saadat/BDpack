% clc;clear;

%% The program to find h*

Nb=4; % number of beads
Ree=13.85*10^-3;%(sqrt(5253)*0.001803); % end-to-end distance
Lc=79.4*10^-3;%9.47; % contour length
ls=Lc/(Nb-1); % segment length
dc=0.5*10^-3; % chain diameter

% The factor to convert hstar to a
afctr=sqrt(pi/3)*Ree/sqrt(Nb-1);
% afctr=sqrt(pi)*0.0637;

% initial guess for hstar
hstar_0=0.01;
a_0=hstar_0*afctr;


hstar_1=0.2;
a_1=hstar_1*afctr;

dhstar=hstar_1-hstar_0;

count=0;
while (abs(dhstar) > 1e-5)
    
    f_0=hsfunc(Nb,Lc,ls,dc,a_0);
    f_1=hsfunc(Nb,Lc,ls,dc,a_1);
      
    hstar_new=(hstar_0*f_1-hstar_1*f_0)/(f_1-f_0);

    dhstar=hstar_new-hstar_1;
    hstar_0=hstar_1;
    a_0=hstar_0*afctr;
    hstar_1=hstar_new;
    a_1=hstar_1*afctr;

    count=count+1;
end

display('Number of iterations:');display(count);
display('h*:');display(hstar_new);
