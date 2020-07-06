function fval=hsfunc(Nb,Lc,ls,dc,a)

% Constructing the diffusion tensor
% The terms are normalized by kBT/eta_s
D=zeros(Nb,Nb);
x=zeros(Nb);

for i=1: Nb
    for j=1: Nb
        
        if (i == j)
            D(i,j)=1/(6*pi*a);
        else
            D(i,j)=1/(6*pi)*3/(4*abs(i-j)*ls)*...
                ( (1+2*a^2/(3*(i-j)^2*ls^2)) + (1-2*a^2/((i-j)^2*ls^2)) );
        end
        
    end % j
    x(i)=-Lc/2+(i-1)*ls;
end % i

Fb=-D^-1*x;

% The spring force on the center of the chain
fspr=0;
for i=floor(Nb/2)+1: Nb
    fspr=fspr+Fb(i);
end

% f=fspr+fdrag
fval=fspr+1/8*2*pi/log(Lc/dc)*Lc^2;

end 