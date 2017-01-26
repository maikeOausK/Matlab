%Calculation of the density matrix
function Con_Nm1_N = red_density_matrix_Nm1_N(Psi,N)
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

    psi_J_1 = one;
    psi_J_2 = one;
   
    
    for cnt2 = 1:cnt -1
     
        psi_J_1 = kron(one,psi_J_1);
        psi_J_2 = kron(one,psi_J_2);
                
    end
      psi_J_1 = kron(up,psi_J_1);
       psi_J_2 = kron(down,psi_J_2);

          Rho_r = psi_J_1' * Rho_r * psi_J_1 + psi_J_2' * Rho_r * psi_J_2;

    end




    sigsig = kron(sigma_y,sigma_y);

    rho_tilt = sigsig*Rho_r*sigsig;

    EV = sort(eig(Rho_r*rho_tilt),'descend');

    Con_Nm1_N = max(0,sqrt(EV(1))-sqrt(EV(2))-sqrt(EV(3))-sqrt(EV(4)));

   % R = sqrtm(sqrtm(Rho_r)*rho_tilt*sqrtm(Rho_r));
   % EV_R = sort(eig(R),'descend');
   % Con_Nm1_N  = (EV_R(1)-EV_R(2)-EV_R(3)-EV_R(4));
    
    Con_Nm1_N = real(Con_Nm1_N);
    if Con_Nm1_N <= 0
        Con_Nm1_N = 0;
    end
end

