%
% 08.11.20216
%-------------
%
% Simulation of a spin chain consisting of two alternating
% lattices.
% The system is described by a realistic Hamiltonian
%
% H = Omega * \sum_k sigma_x^k + V_o * \sum_{m>k } n_m*n_k/(m-k)^gamma
%
% Omega: Rabi frequency
% V_o: interaction strength
% ( zero detuning)

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







% Parameters  

V_0   = 20;
gamma = 6;  


Omega_vec=[0.001:0.01:1 1.1:0.1:20];

Concurrence_Nm1_N(1,1:length(Omega_vec)) = 0;
Concurrence_Nm2_N(1,1:length(Omega_vec)) = 0;
Concurrence_Nm3_N(1,1:length(Omega_vec)) = 0;
Concurrence_2_N(1,1:length(Omega_vec)) = 0;
Concurrence_1_N(1,1:length(Omega_vec)) = 0;
Concurrence_1_2(1,1:length(Omega_vec)) = 0;


 % Concurrence_2_4(1,1:length(Omega_vec)) = 0;

% Construction of the Hamiltonian and the unitary operator
%---------------------------------


 % Construction of all number operators n_k
 
 integer = 1;
 for N = 5
     
     
    a = 1.0; % spacing between the spins (same disctance between the two chains)
    %N = 7;   % number of atoms in the chain (odd number) 
    x = a*(1:N); % spin chain

    % Initial State 
    init = down;
    for cnt = 1:N-1
          init   = sparse(kron(down,init)); % initial state vector
    end

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
    %          U(Omega*t = pi/4)
    %          t = pi/(4*Omega)

    for counter = 1:length(Omega_vec)

       % Rabi frequency
       Omega = Omega_vec(counter);

        t = pi/(4*Omega);
        Hamil_odd = sparse(2^N,2^N);
        V_int_o = sparse(2^N,2^N); %interaction Hamiltonian acting on the odd atoms

        for cnt = 1:2:N 


            Hamil_odd = Hamil_odd + Omega* sig_y_k((cnt-1)*2^N+1:cnt*2^N,:);
            for cnt2 = 1:N

                  if cnt2 == cnt
                       continue
                  end

                  V_int_o = sparse(V_int_o + V_0/(abs(x(cnt)-x(cnt2))^gamma) * ...
                         num_op_m((cnt-1)*2^N+1:cnt*2^N,:)*num_op_m((cnt2-1)*2^N+1:cnt2*2^N,:)); 
            end

        end

        Hamil_odd = Hamil_odd + V_int_o;

        Unitary_pi_2 = sparse(expm(-i*t*Hamil_odd));

        Psi = Unitary_pi_2 *init;

        % 2nd step: pi pulse on all odd atoms
        % --------- 
        %          U(Omega*t = pi/2)
        %          t = pi/(2*Omega)

        t2 = pi/(2*Omega);
        Hamil_even =  sparse(2^N,2^N);
        V_int_e    = sparse(2^N,2^N) ;%interaction Hamiltonian acting on the even atoms

        for cnt = 2:2:N


            Hamil_even = Hamil_even + Omega* sig_y_k((cnt-1)*2^N+1:cnt*2^N,:);
            for cnt2 = 1:N

                  if cnt2 == cnt
                       continue
                  end

                  V_int_e = sparse(V_int_o + V_0/(abs(x(cnt)-x(cnt2))^gamma) * ...
                         num_op_m((cnt-1)*2^N+1:cnt*2^N,:)*num_op_m((cnt2-1)*2^N+1:cnt2*2^N,:)); 
            end

        end
        Hamil_even = Hamil_even + V_int_e;
        Unitary_pi = sparse(expm(-i*t2*Hamil_even));
        Psi = Unitary_pi * Psi;


        % 3rd step: pi/2 pulse on all odd atoms
        % --------- 
        %          U(Omega*t = pi/4)
        %          t = pi/(4*Omega)
% 
%          Unit_pi_2_dagger =  ctranspose(Unitary_pi_2); %;
%         Psi = Unit_pi_2_dagger* Psi;


%          %reduced density matrix rho_{N-1,N}
%          Con_Nm1_N = red_density_matrix_Nm1_N(Psi,N);
%          Concurrence_Nm1_N(integer,counter) =   Con_Nm1_N;
%          
%          %reduced density matrix rho_{1,2}
%          Con_1_2 = red_density_matrix_1_2(Psi,N);
%          Concurrence_1_2(integer,counter) =   Con_1_2;
% 
%          % % reduced density matrix rho_{N-2,N}
%         Con_Nm2_N = red_density_matrix_Nm2_N(Psi,N) ;
%         Concurrence_Nm2_N(integer,counter) =   Con_Nm2_N;
% 
%         % reduced density matrix rho_{N-3,N}
%         if N >3
%         % 
%            Con_Nm3_N =red_density_matrix_Nm3_N(Psi,N);
%            Concurrence_Nm3_N(integer,counter) =  Con_Nm3_N;
%         % 
%         end
% 
%         if N >=5
%          Con_2_4 = red_density_matrix_2_4(Psi,N);
%           Concurrence_2_4(counter) =  Con_2_4;
%         end
%          % % reduced density matrix rho_{2,N}
%         Con_2_N = red_density_matrix_2_N(Psi,N); 
%         Concurrence_2_N(integer,counter) =  Con_2_N;
% 
%         % % % reduced density matrix rho_{1,N}
%          Con_1_N = red_density_matrix_1_N(Psi,N) ;
%          Concurrence_1_N(integer,counter) =  Con_1_N;
% %      
%      end
    integer = integer + 1;
 end
% 
%     figure
%   plot(Omega_vec ,Concurrence_Nm1_N(1,:),Omega_vec, Concurrence_1_2(1,:),Omega_vec ,Concurrence_Nm1_N(2,:),Omega_vec, Concurrence_1_2(2,:),...
%       Omega_vec ,Concurrence_Nm1_N(3,:),Omega_vec, Concurrence_1_2(3,:))
%   legend('C_{N-1,N} N = 5','C_{12},N = 5','C_{N-1,N}, N = 7','C_{12}, N = 7','C_{N-1,N}, N = 9','C_{12}, N = 9')
    
%  plot(Omega_vec,Concurrence_2_4)
%   legend('C_{2,4}')
%   xlabel('Omega')
%   ylabel('Concurrence')
% %   title('N = 7')
% 

%   figure
%   plot(Omega_vec ,Concurrence_Nm1_N(1,:),'+r',Omega_vec ,Concurrence_Nm1_N(2,:),'db', Omega_vec ,Concurrence_Nm1_N(3,:),'c')
%   legend('N = 5','N = 7','N = 9')
%   xlabel('Omega')
%   ylabel('Concurrence')
%   title('C_{N-1,N}')
%   
%   figure
%   plot(Omega_vec ,Concurrence_Nm2_N(1,:),'+r',Omega_vec ,Concurrence_Nm2_N(2,:),'db', Omega_vec ,Concurrence_Nm2_N(3,:),'c')
%   legend('N = 5','N = 7','N = 9')
%   xlabel('Omega')
%   ylabel('Concurrence')
%   title('C_{N-2,N}')
%   
%   figure
%   plot(Omega_vec ,Concurrence_Nm3_N(1,:),'+r',Omega_vec ,Concurrence_Nm3_N(2,:),'db', Omega_vec ,Concurrence_Nm3_N(3,:),'c')
%   legend('N = 5','N = 7','N = 9')
%   xlabel('Omega')
%   ylabel('Concurrence')
%   title('C_{N-3,N}')
%   
%   figure
%   plot(Omega_vec ,Concurrence_2_N(1,:),'+r',Omega_vec ,Concurrence_2_N(2,:),'db', Omega_vec ,Concurrence_2_N(3,:),'c')
%   legend('N = 5','N = 7','N = 9')
%   xlabel('Omega')
%   ylabel('Concurrence')
%   title('C_{2,N}')
% 
%   figure
%   plot(Omega_vec ,Concurrence_1_N(1,:),'+r',Omega_vec ,Concurrence_1_N(1,:),'db', Omega_vec ,Concurrence_1_N(1,:),'c')
%   legend('N = 5','N = 7','N = 9')
%   xlabel('Omega')
%   ylabel('Concurrence')
%   title('C_{1,N}')
%   
% if N <=3
%     plot(Omega_vec ,Concurrence_Nm1_N,'g', Omega_vec,Concurrence_Nm2_N,'b',...
%     Omega_vec,Concurrence_2_N,'k',Omega_vec,Concurrence_1_N,'c')
%     legend('C_{N-1,N}','C_{N-2,N}','C_{2,N}','C_{1,N}')
%     xlabel('Omega')
%     ylabel('Concurrence')
% 
% % % elseif N>5
%   figure(1)
%     plot(Omega_vec ,Concurrence_Nm1_N,'g', Omega_vec,Concurrence_Nm2_N,'b',...
%     Omega_vec,Concurrence_2_N,'k',Omega_vec,Concurrence_1_N,'c',Omega_vec,Concurrence_2_4,'r')
% hold on;
% plot(Omega_vec,Concurrence_1_2,'Color',[0.3 0.6 0.2])
%     legend('C_{N-1,N}','C_{N-2,N}','C_{2,N}','C_{1,N}','C_{2,4}','C_{1,2}')
%     xlabel('Rabi frequency \Omega')
%     ylabel('Concurrence')
%     

% end
% 
