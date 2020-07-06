function chi=chi_FENE(b)

  chi=sqrt( b/3 * ( quadl(@(q) integrand_nom_FENE(q,b),0,1,1e-15) / ...
                  ( quadl(@(q) integrand_denom_FENE(q,b),0,1,1e-15)) ) );  

end

function func_nom=integrand_nom_FENE(q,b)

  func_nom=q.^4.*( (1-q.^2).^(b/2) );
  
end

function func_denom=integrand_denom_FENE(q,b)
  
  func_denom=q.^2.*( (1-q.^2).^(b/2) );
    
end