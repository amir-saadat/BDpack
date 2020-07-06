function chi=chi_WLC(b)

  chi=sqrt( b/3 * ( integral(@(q) integrand_nom_WLC(q,b),0,1,'AbsTol',1e-15) / ...
                  ( integral(@(q) integrand_denom_WLC(q,b),0,1,'AbsTol',1e-15)) ) );  

end

function func_nom=integrand_nom_WLC(q,b)

  func_nom=q.^4.*exp( -b/6*(2*q.^2+1./(1-q)-q) );
  
end

function func_denom=integrand_denom_WLC(q,b)
  
  func_denom=q.^2.*exp( -b/6*(2*q.^2+1./(1-q)-q) );
    
end