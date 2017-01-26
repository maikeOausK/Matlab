%
% 08.11.20216
%-------------
%
% Simulation of a spin chain consisting of two alternating
% lattices.
% The system is described the Hamiltonian
%
% H = Omega * \sum_k phi* sigma_x^k (1-n_{k-1})  (1-n_{k+1})
%
% Omega: Rabi frequency
%
% The length of the pulse is parameterised
%
% 1st Pulse: Omega*t = phi
% 2nd Pulse Omega*t = chi
% 3rd pulse: Omega*t = phi

clear all;
format long;
% Construction of the spin chain


% Single particle operators and states

up   = [1;0];
down = [0;1];
plus = 1/sqrt(2) *[1;1];

sigma_x = [0,1;1,0];
sigma_y = [0,-i;i,0];

num     = [1,0; 0 ,0];
one     = eye(2);

    phi = [0:0.1:2*pi];
    chi = [0:0.1:2*pi];


Concurrence_5_4(length(chi),1:length(phi)) = 0;
Concurrence_5_4_measured(length(chi),1:length(phi)) = 0;



% Construction of the Hamiltonian and the unitary operator
%---------------------------------


 % Construction of all number operators n_k
    
     
    a = 1.0; % spacing between the spins (same disctance between the two chains)
    N = 5;
    x = a*(1:N); % spin chain

    % Initial State 
    init = down;
    for cnt = 1:N-1
          init   = sparse(kron(down,init)); % initial state vector
    end


 num_op_m = sparse(N*2^N,2^N); % all N-particle number operators n_k are stored in this matrix
 help_vec_n(1:2,1:2*N) = 0;    % helping vector to produce the N particle sigma_x operator
 
 sig_y_k = sparse(N*2^N,2^N); % all N-particle Sigma_y are stored in this matrix
 help_vec_s(1:2,1:2*N) = 0;  % helping vector to produce the N particle sigma_x operator
      for cnt = 1:N

              %Construction of number operators     
              for cnt2 = 1:N

                   if cnt2 == cnt
                          help_vec_n(:,2*cnt2-1:2*cnt2) = num;
                          continue
                   end
                   help_vec_n(:,2*cnt2-1:2*cnt2) = one;

              end

              n_k = sparse(help_vec_n(:,2*N-1:2*N));

              for cnt3 = N-1:-1:1

                  n_k =sparse(kron(help_vec_n(:,2*cnt3-1:2*cnt3),n_k));

              end
              num_op_m((cnt-1)*2^N+1:cnt*2^N,:) = n_k;

              %Many particle sigma_x_k operator
                for cnt2 = 1:N
                     if cnt2 == cnt
                          help_vec_s(:,2*cnt2-1:2*cnt2) = sigma_y;
                          continue
                     end
                     help_vec_s(:,2*cnt2-1:2*cnt2) = one;

                end

                s_yk = help_vec_s(:,2*N-1:2*N);

                for cnt3 = N-1:-1:1

                  s_yk =sparse(kron(help_vec_s(:,2*cnt3-1:2*cnt3),s_yk));

                end
                sig_y_k ((cnt-1)*2^N+1:cnt*2^N,:) = s_yk;
      end


    % 1st step: pi/2 pulse on all odd atoms
    % --------- 
    %          U(Omega*t = phi)
    %          t = pi/(4*Omega)
    for cnt_chi = 1:length(chi)

        for counter = 1:length(phi)

           % Rabi frequency

            Hamil_odd = sparse(2^N,2^N);
            for cnt = 1:2:N 
              

                if (cnt >1 && cnt <N   )
                    
                       Hamil_odd = Hamil_odd + sig_y_k((cnt-1)*2^N+1:cnt*2^N,:)...
                             *(eye(2^N)- num_op_m((cnt-2)*2^N+1:(cnt-1)*2^N,:))...
                             *(eye(2^N)- num_op_m((cnt)*2^N+1:(cnt+1)*2^N,:));
                elseif cnt ==1
                       Hamil_odd = Hamil_odd + sig_y_k((cnt-1)*2^N+1:cnt*2^N,:)...
                         *(eye(2^N)- num_op_m((cnt)*2^N+1:(cnt+1)*2^N,:));
                elseif cnt == N
                         Hamil_odd = Hamil_odd + sig_y_k((cnt-1)*2^N+1:cnt*2^N,:)...
                             *(eye(2^N)- num_op_m((cnt-2)*2^N+1:(cnt-1)*2^N,:));...

                end

            end


            Unitary_phi = sparse(expm(-i*phi(counter)*Hamil_odd));         
            Psi = Unitary_phi *init;
            
            
%           %%%%%%%%%%%%%%%%%%%%No interaction Hamiltonian%%%%%%%%%%%%%%%%%%%
%              Hamil_odd_no_int = sparse(2^N,2^N);
%             for cnt = 1:2:N 
%               
% 
%                 if (cnt >1 && cnt <N   )
%                     
%                        Hamil_odd_no_int = Hamil_odd_no_int + sig_y_k((cnt-1)*2^N+1:cnt*2^N,:);
%                      
%                 elseif cnt ==1
%                        Hamil_odd_no_int = Hamil_odd_no_int + sig_y_k((cnt-1)*2^N+1:cnt*2^N,:);
%                         
%                 elseif cnt == N
%                          Hamil_odd_no_int = Hamil_odd_no_int + sig_y_k((cnt-1)*2^N+1:cnt*2^N,:);
%                          
% 
%                 end
% 
%             end
%             
%             
%             
%             Unitary_phi_no_int = sparse(expm(-i*phi(counter)*Hamil_odd_no_int));
%             Psi2 = Unitary_phi_no_int *init;
%             
            
            
            % 2nd step: pi pulse on all even atoms
            % --------- 
            %          U(Omega*t = phi)
            %          t = pi/(2*Omega)
            
        Hamil_even =  sparse(2^N,2^N);
        for cnt = 2:2:N

             if (cnt >1 && cnt <N   )
                       Hamil_even = Hamil_even + sig_y_k((cnt-1)*2^N+1:cnt*2^N,:)...
                             *(eye(2^N)- num_op_m((cnt-2)*2^N+1:(cnt-1)*2^N,:))...
                             *(eye(2^N)- num_op_m((cnt)*2^N+1:(cnt+1)*2^N,:));
                elseif cnt ==1
                       Hamil_even = Hamil_even + sig_y_k((cnt-1)*2^N+1:cnt*2^N,:)...
                        *(eye(2^N)- num_op_m((cnt)*2^N+1:(cnt+1)*2^N,:));
                elseif cnt == N
                         Hamil_even = Hamil_even + sig_y_k((cnt-1)*2^N+1:cnt*2^N,:)...
                             *(eye(2^N)- num_op_m((cnt-2)*2^N+1:(cnt-1)*2^N,:));...

                end
        end
            
    %         
            Unitary_chi = sparse(expm(-i*chi(cnt_chi)*Hamil_even));
            Psi = Unitary_chi * Psi;

      
              % Measurement in basis M_i = 1/sqrt(2) *(|0> +- |1>)
              M  = one; %Measurement
            
              for cnter = 1:N-2
                  M = kron(M,one);

              end
              ONE = kron(up,M)  ;
              ZERO = kron(down,M);
              
              M = 1/sqrt(2)*(exp(i*pi/2)*ONE + ZERO);
              Psi_new = M' * Psi/(sqrt(Psi'*M*M'*Psi));

    %          %reduced density matrix rho_{N-1,N}
             Con_Nm1_N = red_density_matrix_Nm1_N(Psi,N);
             Concurrence_5_4(cnt_chi,counter) = Con_Nm1_N;

             
            %reduced density matrix rho(5_4) after the measurement 
             Rho_red = Psi_new*Psi_new';
            %-----------------------------------------------------------------------
            % reduced density matrix rho_{N-1,N}
            Rho_r = Rho_red;
    
            for cnt = N-2:-1:2

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

        

            Con_Nm1_N = real(Con_Nm1_N);
            if Con_Nm1_N <= 0
                Con_Nm1_N = 0;
            end
          Concurrence_5_4_measured(cnt_chi,counter) = Con_Nm1_N;
        end

    end



%   figure
% surf(phi,chi,Concurrence_5_4)
% xlabel('\phi')
% ylabel('\chi')
% zlabel('C_{54}')
% xlim([0 2*pi])
% ylim([0 2*pi])
% colormap jet
% colorbar
% 
% 
% figure
% contourf(phi,chi,Concurrence_5_4)
% xlabel('\phi')
% ylabel('\chi')
% title('C_{54}')
% colormap jet
% colorbar


  figure
surf(phi,chi,Concurrence_5_4_measured)
xlabel('\phi')
ylabel('\chi')
zlabel('C_{54}')
xlim([0 2*pi])
ylim([0 2*pi])
colormap jet
colorbar


figure
contourf(phi,chi,Concurrence_5_4_measured)
xlabel('\phi')
ylabel('\chi')
title('C_{54}')
colormap jet
colorbar




