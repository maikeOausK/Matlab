
% Construction of a GHZ state
% GHZ = |+-+-+-...> +|-+-+-+...>
% |+> = spin up
% |-> = spin down
% The initial state: |Psi> = |++++...>
%
% The system is treated as an effective 2 state system
% with the states |R> and |+>  
%(it is assumed that the (pi,-) pulse works perfectly and all occupation
% of the |R> state is transferred to the |-> state => the (pi,-) is neglected
% and in the end, the occupation of |R> can be interpreted as the 
% occupation of the |-> state)

% Hamiltonian of the system:
% H = sum_k Omega_k * sigma_x^(k) + sum_k Delta_k * n_k + sum_k,m V_{km} * n_k * n_m

% Omega_k : Rabi frequency of particle k
% Delta_k : Detuning of the rydberg state of particle k
% V_{km}  : interaction potential
% V_km = V_0 * 1/|k-m|^gamma
%Approximation
% Delta_k = 0 -> no Detuning
% H = sum_k Omega_k * sigma_x^(k) + sum_k,m V_{km} * n_k * n_m

clear all;
format long;

a = 1.0;     % distance between the atoms
    % Definition of the states
ryd  = [1;0]; % Rydberg state
up   = [0;1]; % up state

% single particle operators
num     = ryd*ryd'; % number operator n = |R><R|
sigma_x = [0, 1;1,0]; %
one = eye(2);% identity matrix

% many particle state
gamma = 6; % exponent of interaction force


Omega_vec=[0.001:0.01:1 1.1:0.1:20];
Fidelity(1:3,1:length(Omega_vec)) = 0;
integer = 1;



for N = 4:2:8

    x = a*(1:N); % position of the atoms
    
    num_op_m=sparse(N*2^N,2^N);  % all N particle number operators n_k are stored in this matrix
    help_vec_n(1:2,1:2*N) = 0;  % helping vector to produce the N particle sigma_x operator



    O_V(1:length(Omega_vec))       = 0;
    GS = up;
    % Construction of the many particle GS
    for cnt =1:N-1
              GS = kron(GS,up);  %construction of many particle product state
    end

    % Construction of all number operators n_k
    for cnt = 1:N
            % Construction of number operators     
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

    end
    
       %Construction of the correct GHZ state

        GHZ_l = up;
        GHZ_r = ryd;

        for cnt = 1:N-1
          if mod(cnt,2) == 0
              GHZ_l = kron(up,GHZ_l);
              GHZ_r = kron(ryd,GHZ_r);

          else
              GHZ_l = kron(ryd,GHZ_l);
              GHZ_r = kron(up,GHZ_r);
          end
        end

        GHZ = (1/sqrt(2)) * (GHZ_l+GHZ_r);

    for counter = 1:length(Omega_vec) % Measuring the fidelity



        % Rabi frequency for the two transitions
        Omega = Omega_vec(counter);
        % Interaction potential
        V_0 = 20; % nearest neighbor interaction

        O_V (counter) = Omega/V_0;

        Psi= GS;

        for k_atom = 1:N

            %Many particle sigma_x_k operator
            help_vec_s(1:2,1:2*N) = 0; %helping vector to produce the N particle sigma_x operator


            for cnt2 = 1:N
                 if cnt2 == k_atom
                      help_vec_s(:,2*cnt2-1:2*cnt2) = sigma_x;
                      continue
                 end
                 help_vec_s(:,2*cnt2-1:2*cnt2) = one;

            end

            s_xk = help_vec_s(:,2*N-1:2*N);

            for cnt3 = N-1:-1:1

              s_xk =sparse(kron(help_vec_s(:,2*cnt3-1:2*cnt3),s_xk));

            end


            % Interaction Potential
            V_int = sparse(2^N,2^N) ;
            Hamil_k = sparse(2^N,2^N);
            for cnt = 1:N

                if cnt == k_atom 
                    continue
                end

                 V_int = sparse(V_int + V_0/(abs(x(k_atom)-x(cnt))^gamma) * ...
                 num_op_m((k_atom-1)*2^N+1:k_atom*2^N,:)*num_op_m((cnt-1)*2^N+1:cnt*2^N,:)); %* num_all others       
            end
           
            t = pi/(2*Omega);
            Hamil_k = sparse(Omega * s_xk + V_int);
            
            if k_atom == 1
                t = pi/(4*Omega);
            end

            Psi = expm(-1i*t*Hamil_k)*Psi;

        end


        Fidelity(integer,counter) = abs(Psi'*GHZ);
    end
    integer = integer+1;
end
plot( O_V,Fidelity(1,:),'xg', O_V,Fidelity(2,:),'+m', O_V,Fidelity(3,:))
grid on;
grid minor;
legend('L=4 (2)','L=6 (2)','L=8 (2)')
xlabel('\Omega/V_0')
ylabel('Fidelity F')

% %--------------------------------------------------------------------------
% GHZ State as a three state scheme
%
%

Fidelity_3_state(1:3,1:length(Omega_vec)) = 0;
integer = 1;
%
for N = 4:2:8
    % Definition of the chain

    x = a*(1:N); % position of the atoms
   
    %Construction of the correct GHZ state

    GHZ_l = up;
    GHZ_r = ryd;

    for cnt = 1:N-1
          if mod(cnt,2) == 0
              GHZ_l = kron(up,GHZ_l);
              GHZ_r = kron(ryd,GHZ_r);

          else
              GHZ_l = kron(ryd,GHZ_l);
              GHZ_r = kron(up,GHZ_r);
          end
    end

     GHZ = (1/sqrt(2)) * (GHZ_l+GHZ_r);
    
    for counter = 1:length(Omega_vec)% Measuring the fidelity
        
        % Rabi frequency for the two transitions
        Omega = Omega_vec(counter);
        % Interaction potential
        O_V (counter) = Omega/V_0;
        GS = up;
        % Construction of the many particle GS
        for cnt =1:N-1
                  GS = kron(GS,up);  %construction of many particle product state
        end


        s_x1   = kron(sigma_x,one); %sigma_x1
        s_x2   = kron(one,sigma_x); %sigma_x2

        n_k    = kron(num,one); %n_k
        nk_1   = kron(one,num); %n_k+1

        V_int_12 =  sparse(V_0/((abs(a))^gamma) * n_k*nk_1); %* num_all others       
        V_int_21 =  sparse(V_0/((abs(a))^gamma) * nk_1*n_k);
        %

        % pi/2 pulse to produce superposition state -> acts on 1st atom
        Hamil_12 = sparse(Omega*s_x1 + V_int_12);
        t        = pi/(4*Omega);
        Unit_12  = sparse(expm(-i*t*Hamil_12)); % unitary operator
        if N > 2
            unt = one;
            for cnt = 1:N-3
              unt = sparse(kron(one,unt));
            end
            Unit_12 = sparse(kron(Unit_12,unt));
        end

        Psi = Unit_12*GS;

        % Series of pi pulses
        Hamil    = Omega*s_x2 + V_int_21;

        for k_atom = 2:N
          unt = one;
          t      = pi/(2*Omega);
          Unit   = expm(-i*t*Hamil); % unitary operator  

          if k_atom ~= N
              for cnt = N-1:-1:2

                  if cnt ==k_atom 
                     unt = sparse(kron(Unit,unt)) ;
                     continue
                  end
                  unt = sparse(kron(one,unt));

              end
              Psi = unt * Psi;
          elseif k_atom ==N 

              unt = sparse(kron(one,Unit));
              for cnt =  N-2:-1:2
                  unt = sparse(kron(one,unt));
              end
              Psi= unt*Psi;
          end

        end
    Fidelity_3_state(integer,counter) = abs(Psi'*GHZ);
    end
       integer = integer+1;
end
plot( O_V,Fidelity(1,:),'xg', O_V,Fidelity(2,:),'+m', O_V,Fidelity(3,:),'cd',O_V,Fidelity_3_state(1,:),'b',...
    O_V,Fidelity_3_state(2,:),'k',O_V,Fidelity_3_state(3,:),'r')
grid on;
grid minor;
legend('L=4 (2)','L=6 (2)','L=8 (2)','L=4 (3)','L=6 (3)', 'L=8 (3)')
xlabel('\Omega/V_0')
ylabel('Fidelity F')