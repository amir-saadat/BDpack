function b_val=bfunc(Nk,Nb,method)

  b_val=fsolve(@(b) b_chi2func(Nk,Nb,method,b),1);

end

function func=b_chi2func(Nk,Nb,method,b)

%   term=Nk*3*(Nb+1)/(Nb*(Nb-1));
  term=Nk*3/(Nb-1);
  if (strcmp(method,'ILC'))
      func=b/chi_ILC(b)^2-term;
  end
  
  if (strcmp(method,'RWS'))
      N_ks=Nk/(Nb-1);
      func=b/chi_RWS(b,N_ks)^2-term;
  end  
  
  if (strcmp(method,'FENE'))
      func=b/chi_FENE(b)^2-term;
  end
  
  if (strcmp(method,'WLC'))
      func=b/chi_WLC(b)^2-term;
  end
  
  if (strcmp(method,'WLC_UD'))
      N_ks=Nk/(Nb-1);
      func=b/chi_WLC_UD(b,N_ks)^2-term;
  end
  
end