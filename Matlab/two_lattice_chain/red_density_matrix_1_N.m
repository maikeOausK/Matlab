function Con_1_N = red_density_matrix_1_N(Psi,N) 

    up   = [1;0];
    down = [0;1];
    plus = 1/sqrt(2) *[1;1];
    one     = eye(2);
    sigma_y = [0,-i;i,0];
    
    Rho_red = Psi*Psi';
    Rho_r_1_N = Rho_red;

for cnt = N-1:-1:2
  
    
    psi_J_1 = kron(up,one);
    psi_J_2 = kron(down,one);
    
    for cnt2 = 1:cnt-1 
     
        psi_J_1 = kron(one,psi_J_1);
        psi_J_2 = kron(one,psi_J_2);
                
    end
   
      Rho_r_1_N = psi_J_1' * Rho_r_1_N * psi_J_1 + psi_J_2' * Rho_r_1_N * psi_J_2;

end

sigsig = kron(sigma_y,sigma_y);

rho_tilt = sigsig*Rho_r_1_N*sigsig;

EV = sort(eig(Rho_r_1_N*rho_tilt),'descend');

Con_1_N = max(0,(sqrt(EV(1))-sqrt(EV(2))-sqrt(EV(3))-sqrt(EV(4))));
% 
% R = sqrtm(sqrtm(Rho_r_1_N)*rho_tilt*sqrtm(Rho_r_1_N));
% EV_R = sort(eig(R),'descend');
% Con_1_N  = (EV_R(1)-EV_R(2)-EV_R(3)-EV_R(4));

Con_1_N = real(Con_1_N);
if Con_1_N <= 0
    Con_1_N = 0;
end
    
end