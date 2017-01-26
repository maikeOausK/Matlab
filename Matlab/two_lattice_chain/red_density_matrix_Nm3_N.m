function Con_Nm3_N = red_density_matrix_Nm3_N(Psi,N) 
% % reduced density matrix rho_{N-3,N}
% Rho_r_N-3_N = Rho_red;
    up   = [1;0];
    down = [0;1];
    plus = 1/sqrt(2) *[1;1];
    one     = eye(2);
    sigma_y = [0,-i;i,0];
    
    Rho_red = Psi*Psi';
    
        Rho_r_Nm3_N = Rho_red;
    
        for cnt = N-1:-1:4
    
    
            psi_J_1 = one;
            psi_J_2 = one;
    
            for cnt2 = 1:cnt-1
    
                psi_J_1 = kron(one,psi_J_1);
                psi_J_2 = kron(one,psi_J_2);
    
            end
    
               psi_J_1 = kron(up,psi_J_1);
               psi_J_2 = kron(down,psi_J_2);
               Rho_r_Nm3_N = psi_J_1' * Rho_r_Nm3_N * psi_J_1 + psi_J_2' * Rho_r_Nm3_N * psi_J_2;
    
        end
        
    psi_J_1 = kron(one,one);
    psi_J_1 = kron(up,psi_J_1);
    psi_J_1 = kron(one,psi_J_1);
    psi_J_2 = kron(one,one);
    psi_J_2 = kron(down,psi_J_2);
    psi_J_2 = kron(one,psi_J_2);
    Rho_r_Nm3_N = psi_J_1' * Rho_r_Nm3_N * psi_J_1 + psi_J_2' * Rho_r_Nm3_N * psi_J_2;
    
    
    psi_J_1 = kron(up,one);
    psi_J_1 = kron(one,psi_J_1);
    
    psi_J_2 = kron(down,one);
    psi_J_2 = kron(one,psi_J_2);
    Rho_r_Nm3_N = psi_J_1' * Rho_r_Nm3_N * psi_J_1 + psi_J_2' * Rho_r_Nm3_N * psi_J_2;
    
    
    sigsig = kron(sigma_y,sigma_y);
    
    rho_tilt = sigsig*Rho_r_Nm3_N*sigsig;
    
    EV = sort(eig(Rho_r_Nm3_N*rho_tilt),'descend');
    
   Con_Nm3_N = max(0,sqrt(EV(1))-sqrt(EV(2))-sqrt(EV(3))-sqrt(EV(4)));
%     
%     R = sqrtm(sqrtm(Rho_r_Nm3_N)*rho_tilt*sqrtm(Rho_r_Nm3_N));
%     EV_R = sort(eig(R),'descend');
%     Con_Nm3_N  = max(0,EV_R(1)-EV_R(2)-EV_R(3)-EV_R(4));

    Con_Nm3_N = real(Con_Nm3_N);
    if Con_Nm3_N <= 0
        Con_Nm3_N = 0;
    end
end