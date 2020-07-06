function chi=chi_WLC_UD(b,N_ks)
  
  A=3/32-3/(8*N_ks)-3/(2*N_ks^2);
  B=(13/32+0.4086/N_ks-14.79/(4*N_ks^2))/(1-4.225/(2*N_ks)+4.87/(4*N_ks^2));
  
  chi=sqrt( b/3 * ( integral(@(q) integrand_nom_WLC_UD(q,b,A,B,N_ks),0,1,'AbsTol',1e-15) / ...
                  ( integral(@(q) integrand_denom_WLC_UD(q,b,A,B,N_ks),0,1,'AbsTol',1e-15)) ) );  

end

function func_nom=integrand_nom_WLC_UD(q,b,A,B,N_ks)

  func_nom=q.^4.*( (1-q.^2).^(-7*b/(6*N_ks)).*exp( -2*b/3*( 1./(2-2*q.^2)+A/2*q.^2-B/4*q.^2.*(q.^2-2) )) );
  
end

function func_denom=integrand_denom_WLC_UD(q,b,A,B,N_ks)
  
  func_denom=q.^2.*( (1-q.^2).^(-7*b/(6*N_ks)).*exp( -2*b/3*( 1./(2-2*q.^2)+A/2*q.^2-B/4*q.^2.*(q.^2-2) )) );
    
end