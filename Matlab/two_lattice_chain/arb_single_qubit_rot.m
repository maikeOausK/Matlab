%
% 08.11.20216
%-------------
% PERFORM AN ARBITRARY SINGLE QUBIT ROTATION
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
sigma_z = [1,0;0,-1];

num     = [1,0; 0 ,0];
one     = eye(2);

    phi = pi/4; 
    chi = pi/2;    
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

        end
       
        
    end
    
    
    % Perform arbitrary rotation on the last qubit of Psi
    %
    % Perform measurement in the bases
    % B_j(delta_j) = 1/sqrt(2){|0> + exp(i delta_j) |1> , |0> - exp(i delta_j)} |1>}
    % measurement outcome: s_j in {0,1
    % s_j = 0 -> qubit j projected onto 1st state of B_j(delta_j)
    %
    % ALGORITHM
    %
    % 1. measure qubit 1 in B_1(0)
    % 2. measure qubit 2 in B_2((-1)^{s_1+1} alpha)
    % 3. measure qubit 2 in B_2((-1)^{s_2} beta)
    % 4. measure qubit 2 in B_2((-1)^{s_1+s_3} gamma)
     
     
    % 1. measure qubit 1 in B_1(0)
     
     delta = 0;
     % determin the measurement basis state
     s1 = round(rand); %creates a random number 0 or 1
     
     if s1 == 0
       
         B = (down +exp(i*delta)*up)/sqrt(2);
      
     elseif s1 == 1
        
         B = (down -exp(i*delta)*up)/sqrt(2) ;
         
     end
     Measure = one;
     Measure = kron(Measure,one);
     Measure = kron(Measure,one);
     Measure = kron(Measure,one);
     Measure = kron(B,Measure);
     
     Psi = Measure' * Psi/(sqrt(Psi'*Measure*Measure'*Psi));
     
     % 2. measure qubit 2 in B_2((-1)^(s_1+1)beta)
     
     delta = ((-1)^(s1+1))*pi/2;
     
     % determin the measurement basis state
     s2 = round(rand); %creates a random number 0 or 1
     
     if s2 == 0
       
         B = (down +exp(i*delta)*up)/sqrt(2);
      
     elseif s2 == 1
        
         B = (down -exp(i*delta)*up)/sqrt(2) ;
         
     end
     Measure = one;
     Measure = kron(Measure,one);
     Measure = kron(Measure,one);
     Measure = kron(B,Measure);
     
     Psi = Measure' * Psi/(sqrt(Psi'*Measure*Measure'*Psi));
     
     
    % 3. measure qubit 3 in B_3((-1)^s_2*gamma)
     
     delta = ((-1)^s2)*pi/4;
     
     % determin the measurement basis state
     s3 = round(rand); %creates a random number 0 or 1
     
     if s3 == 0
       
         B = (down +exp(i*delta)*up)/sqrt(2) ;
      
     elseif s3 == 1
        
         B = (down -exp(i*delta)*up)/sqrt(2) ;
         
     end
     Measure = one;
     Measure = kron(Measure,one);
     Measure = kron(B,Measure);
     
     Psi = Measure' * Psi/(sqrt(Psi'*Measure*Measure'*Psi));
     
     
     % 4. measure qubit 3 in B_3((-1)^{s_1+s_3}*gamma)
     
     delta = ((-1)^(s1+s3))*pi/3;
     
     % determin the measurement basis state
     s4 = round(rand); %creates a random number 0 or 1
     
     if s4 == 0
       
         B = (down +exp(i*delta)*up)/sqrt(2);
      
     elseif s4 == 1
        
         B = (down -exp(i*delta)*up)/sqrt(2) ;
         
     end
     Measure = one;
     Measure = kron(B,Measure);
     
     Psi = Measure' * Psi/(sqrt(Psi'*Measure*Measure'*Psi));

     Psi = sigma_x^(s2+s4)*sigma_z^(s1+s3)*Psi ;

  
%     EINS = Psi(1)      
%     
%     Z1 = abs(EINS)
%     if real()