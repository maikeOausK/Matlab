function Con_3_6 = red_density_matrix_3_6(Psi,N) 

    up   = [1;0];
    down = [0;1];
    plus = 1/sqrt(2) *[1;1];
    one     = eye(2);
    sigma_y = [0,-i;i,0];
    Rho_red = Psi*Psi';
       
    Rho_r_3_6 = Rho_red;

    for cnt = N-1:-1:N-3


        psi_J_1 = kron(one,up);
        psi_J_2 = kron(one,down);

        for cnt2 = 1:cnt-1

            psi_J_1 = kron(one,psi_J_1);
            psi_J_2 = kron(one,psi_J_2);

        end

       
           Rho_r_3_6 = psi_J_1' * Rho_r_3_6 * psi_J_1 + psi_J_2' * Rho_r_3_6 * psi_J_2;

    end
    
    for cnt = 5:-1:4
        psi_J_1 = kron(up,one);
        psi_J_2 = kron(down,one);

        for cnt2 = 1:cnt-1

            psi_J_1 = kron(one,psi_J_1);
            psi_J_2 = kron(one,psi_J_2);

        end

       
           Rho_r_3_6 = psi_J_1' * Rho_r_3_6 * psi_J_1 + psi_J_2' * Rho_r_3_6 * psi_J_2;
        
    end
    
       for cnt = 3:-1 :2
        psi_J_1 =one;
        psi_J_2 = one;

        for cnt2 = 1:cnt-1

            psi_J_1 = kron(one,psi_J_1);
            psi_J_2 = kron(one,psi_J_2);

        end

            psi_J_1 = kron(up,psi_J_1);
            psi_J_2 = kron(down,psi_J_2);

       
           Rho_r_3_6 = psi_J_1' * Rho_r_3_6 * psi_J_1 + psi_J_2' * Rho_r_3_6 * psi_J_2;
        
    end
    
    
   
   sigsig = kron(sigma_y,sigma_y);

   rho_tilt = sigsig*Rho_r_3_6*sigsig;

   EV = sort(eig(Rho_r_3_6*rho_tilt),'descend');

  Con_3_6 = max(0,sqrt(EV(1))-sqrt(EV(2))-sqrt(EV(3))-sqrt(EV(4)));
% 
%     R = sqrtm(sqrtm(Rho_r_2_4)*rho_tilt*sqrtm(Rho_r_2_4));
%     EV_R = sort(eig(R),'descend');
%     Con_2_4  = max(0,EV_R(1)-EV_R(2)-EV_R(3)-EV_R(4));


Con_3_6 = real(Con_3_6);
if Con_3_6 <= 0
    Con_3_6= 0;
end
    
end
  