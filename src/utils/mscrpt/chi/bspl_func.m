function bspl_val=bspl_func(m,K,p)

prefact=complex(cos(2*pi*(p-1)*m/K),sin(2*pi*(p-1)*m/K));

sum_spl=0.0;
for k=0: p-2,
    
    sum_spl=sum_spl+W_func(k+1,p)*complex(cos(2*pi*m*k/K),sin(2*pi*m*k/K)) ;
    
end,

bspl_val=prefact/sum_spl;