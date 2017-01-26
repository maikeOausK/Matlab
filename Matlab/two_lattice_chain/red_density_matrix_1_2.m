%Calculation of the density matrix
function Con_1_2 = red_density_matrix_1_2(Psi,N)
    Rho_red = Psi*Psi';
    %-----------------------------------------------------------------------
    % reduced density matrix rho_{N-1,N}
    Rho_r = Rho_red;
    % Single particle operators and states

    up   = [1;0];
    down = [0;1];
    plus = 1/sqrt(2) *[1;1];
    one     = eye(2);
    sigma_y = [0,-i;i,0];
    for cnt = N-1:-1:2


 
      
    psi_J_1 = up;
    psi_J_2 =down;
   
    
    for cnt2 = 1:cnt
     
        psi_J_1 = kron(one,psi_J_1);
        psi_J_2 = kron(one,psi_J_2);
                
    end
  
          Rho_r = psi_J_1' * Rho_r * psi_J_1 + psi_J_2' * Rho_r * psi_J_2;

    end




    sigsig = kron(sigma_y,sigma_y);

    rho_tilt = sigsig*Rho_r*sigsig;

    EV = sort(eig(Rho_r*rho_tilt),'descend');

    Con_1_2 = max(0,sqrt(EV(1))-sqrt(EV(2))-sqrt(EV(3))-sqrt(EV(4)));

%     R = sqrtm(sqrtm(Rho_r)*rho_tilt*sqrtm(Rho_r));
%     EV_R = sort(eig(R),'descend');
%     Con_1_2  = (EV_R(1)-EV_R(2)-EV_R(3)-EV_R(4));
    
    Con_1_2 = real(Con_1_2);
    if Con_1_2 <= 0
        Con_1_2 = 0;
    end
end