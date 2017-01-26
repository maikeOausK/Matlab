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


Concurrence_5_4_ni(length(chi),1:length(phi)) = 0;
Concurrence_5_3_ni(length(chi),1:length(phi)) = 0;

Concurrence_5_4(length(chi),1:length(phi)) = 0;
Concurrence_5_3(length(chi),1:length(phi)) = 0;
Concurrence_5_2(length(chi),1:length(phi)) = 0;
Concurrence_5_1(length(chi),1:length(phi)) = 0;
Concurrence_1_2(length(chi),1:length(phi)) = 0;

Concurrence_7_4(length(chi),1:length(phi)) = 0;

 % Concurrence_2_4(1,1:length(Omega_vec)) = 0;

% Construction of the Hamiltonian and the unitary operator
%---------------------------------


 % Construction of all number operators n_k
    
     
    a = 1.0; % spacing between the spins (same disctance between the two chains)
    N = 7;
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

      
            % 3rd step: pi/2 pulse on all odd atoms
            % --------- 
            %          U(Omega*t = pi/4)
            %          t = pi/(4*Omega)


        %   Psi = Unitary_phi* Psi;

    %          %reduced density matrix rho_{N-1,N}
             Con_Nm1_N = red_density_matrix_Nm1_N(Psi,N);
             Concurrence_5_4(cnt_chi,counter) = Con_Nm1_N;

    % 
             % % reduced density matrix rho_{N-2,N}
            Con_Nm2_N = red_density_matrix_Nm2_N(Psi,N) ;
            Concurrence_5_3(cnt_chi,counter) =   Con_Nm2_N;

                   
            
 
             Con_2_4 = red_density_matrix_2_4(Psi,N);
              Concurrence_2_4(cnt_chi,counter) =  Con_2_4;
  
             % % reduced density matrix rho_{2,N}
            Con_2_N = red_density_matrix_2_N(Psi,N); 
            Concurrence_5_2(cnt_chi,counter) =  Con_2_N;
    
            % % % reduced density matrix rho_{1,N}
             Con_1_N = red_density_matrix_1_N(Psi,N) ;
             Concurrence_5_1(cnt_chi,counter) =  Con_1_N;
         
             %reduced density matrix rho_{1,2}
             Con_1_2 = red_density_matrix_1_2(Psi,N);
             Concurrence_1_2(cnt_chi,counter) =   Con_1_2;
         end

    end
%  figure
% surf(phi,chi,abs(Concurrence_5_4-Concurrence_5_4_ni)./Concurrence_5_4)
% xlabel('\phi')
% ylabel('\chi')
% zlabel('C_{54}')
% xlim([0 pi])
% ylim([0 2*pi])
% colormap jet
% colorbar
% 
% 
% figure
% contourf(phi,chi,abs(Concurrence_5_4-Concurrence_5_4_ni)./Concurrence_5_4)
% xlabel('\phi')
% ylabel('\chi')
% title('C_{54}')
% colormap jet
% colorbar



%     figure
% surf(phi,chi,abs(Concurrence_5_3-Concurrence_5_3_ni)./Concurrence_5_3_ni)
% xlabel('\phi')
% ylabel('\chi')
% zlabel('C_{53}')
% xlim([0 pi])
% ylim([0 2*pi])
% colormap jet
% colorbar
% 
% figure
% contourf(phi,chi,abs(Concurrence_5_3-Concurrence_5_3_ni)./Concurrence_5_3_ni)
% xlabel('\phi')
% ylabel('\chi')
% title('C_{53}')
% colormap jet
%  colorbar
 
 
  figure
surf(phi,chi,Concurrence_5_4)
xlabel('\phi')
ylabel('\chi')
zlabel('C_{76}')
xlim([0 2*pi])
ylim([0 2*pi])
colormap jet
colorbar


figure
contourf(phi,chi,Concurrence_5_4)
xlabel('\phi')
ylabel('\chi')
title('C_{76}')
colormap jet
colorbar

    figure
surf(phi,chi,Concurrence_1_2)
xlabel('\phi')
ylabel('\chi')
zlabel('C_{12}')
xlim([0 2*pi])
ylim([0 2*pi])
colormap jet
colorbar
figure
contourf(phi,chi,Concurrence_1_2)
xlabel('\phi')
ylabel('\chi')
title('C_{12}')
colormap jet
colorbar


% figure
% contourf(phi,chi,abs(Concurrence_1_2-Concurrence_5_4)./Concurrence_5_4)
% xlabel('\phi')    figure
surf(phi,chi,Concurrence_2_4)
xlabel('\phi')
ylabel('\chi')
zlabel('C_{24}')
xlim([0 2*pi])
ylim([0 2*pi])
colormap jet
colorbar

figure
contourf(phi,chi,Concurrence_2_4)
xlabel('\phi')
ylabel('\chi')
title('C_{24}')
colormap jet
colorbar
% ylabel('\chi')
% title('relative error')
% colormap jet
% colorbar
% 
%     figure
% surf(phi,chi,abs(Concurrence_1_2-Concurrence_5_4)./Concurrence_5_4)
% xlabel('\phi')
% ylabel('\chi')
% zlabel('relative error ')
% xlim([0 pi])
% ylim([0 2*pi])
% colormap jet
% colorbar





    figure
surf(phi,chi,Concurrence_5_3)
xlabel('\phi')
ylabel('\chi')
zlabel('C_{75}')
xlim([0 2*pi])
ylim([0 2*pi])
colormap jet
colorbar

figure
contourf(phi,chi,Concurrence_5_3)
xlabel('\phi')
ylabel('\chi')
title('C_{75}')
colormap jet
 colorbar
 
 
    figure    figure
surf(phi,chi,Concurrence_2_4)
xlabel('\phi')
ylabel('\chi')
zlabel('C_{24}')
xlim([0 2*pi])
ylim([0 2*pi])
colormap jet
colorbar

figure
contourf(phi,chi,Concurrence_2_4)
xlabel('\phi')
ylabel('\chi')
title('C_{24}')
colormap jet
colorbar
surf(phi,chi,Concurrence_7_4)
xlabel('\phi')
ylabel('\chi')
zlabel('C_{74}')
xlim([0 2*pi])
ylim([0 2*pi])
colormap jet
colorbar
 
 
figure
contourf(phi,chi,Concurrence_7_4)
xlabel('\phi')
ylabel('\chi')
title('C_{74}')
colormap jet
 colorbar

 
% 


    figure
surf(phi,chi,Concurrence_5_2)
xlabel('\phi')
ylabel('\chi')
zlabel('C_{72}')
xlim([0 2*pi])
ylim([0 2*pi])
colormap jet
colorbar

figure
contourf(phi,chi,Concurrence_5_2)
xlabel('\phi')
ylabel('\chi')
title('C_{72}')
colormap jet
colorbar


    figure
surf(phi,chi,Concurrence_5_1)
xlabel('\phi')
ylabel('\chi')
zlabel('C_{71}')
xlim([0 2*pi])
ylim([0 2*pi])
colormap jet
colorbar

figure
contourf(phi,chi,Concurrence_5_1)
xlabel('\phi')
ylabel('\chi')
title('C_{71}')
colormap jet
colorbar



    figure
surf(phi,chi,Concurrence_2_4)
xlabel('\phi')
ylabel('\chi')
zlabel('C_{24}')
xlim([0 2*pi])
ylim([0 2*pi])
colormap jet
colorbar

figure
contourf(phi,chi,Concurrence_2_4)
xlabel('\phi')
ylabel('\chi')
title('C_{24}')
colormap jet
colorbar
