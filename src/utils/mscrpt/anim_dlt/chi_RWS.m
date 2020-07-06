
function chi=chi_RWS(b,N_ks)
  
  C=3-10/(3*N_ks)+10/(27*N_ks^2);
  D=1+2/(3*N_ks)+10/(27*N_ks^2);
  
  chi=sqrt( b/3 * ( quadl(@(q) integrand_nom_RWS(q,b,C,D),0,1,1e-15) / ...
                  ( quadl(@(q) integrand_denom_RWS(q,b,C,D),0,1,1e-15)) ) );  

end

function func_nom=integrand_nom_RWS(q,b,C,D)

  func_nom=q.^4.*( (1-q.^2).^(b*( C/6-D/6 )) .* exp(-b*D/6*q.^2) );
  
end

function func_denom=integrand_denom_RWS(q,b,C,D)
  
  func_denom=q.^2.*( (1-q.^2).^(b*( C/6-D/6 )) .* exp(-b*D/6*q.^2) );
    
end