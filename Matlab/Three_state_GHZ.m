% GHZ State as a three state scheme
%
%
%
clear all;

M = 1000;  % number of fidelity measurements
Fidelity_3_state(1:2,1:M) = 0;
integer = 1;
for N = 4:2:6
    % Definitopm of the chain
    %N =4  ;   % number of atoms
    %L = 4.7;       % length of chain
    a = 1.175;     % distance between the atoms
    x = a*(1:N); % position of the atoms

    % V_0 = 1;
    % Omega = 1;
    % Definition of the states
    ryd  = [1;0]; % Rydberg state
    up   = [0;1]; % up state
    %down = [0;0;1]; % down state

    % single particle operators
    num     = ryd*ryd'; % number operator n = |R><R|
    sigma_x = [0, 1;1,0]; %
    one = eye(2);% identity matrix

    % many particle state
    gamma = 6; % exponent of interaction force

    
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
    
    for counter = 1:M % Measuring the fidelity
        
        % Rabi frequency for the two transitions
        Omega = 0.0000001*counter;
        % Interaction potential
        V_0 = 0.0001; % nearest neighbor interaction

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

        V_int_12 =  V_0/((abs(a))^gamma) * n_k*nk_1; %* num_all others       
        V_int_21 =  V_0/((abs(a))^gamma) * nk_1*n_k;
        %

        % pi/2 pulse to produce superposition state -> acts on 1st atom
        Hamil_12 = Omega*s_x1 + V_int_12;
        t        = pi/(4*Omega);
        Unit_12  = expm(-i*t*Hamil_12); % unitary operator
        if N > 2
            unt = one;
            for cnt = 1:N-3
              unt = kron(one,unt);
            end
            Unit_12 = kron(Unit_12,unt);
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
                     unt = kron(Unit,unt) ;
                     continue
                  end
                  unt = kron(one,unt);

              end
              Psi = unt * Psi;
          elseif k_atom ==N 

              unt = kron(one,Unit);
              for cnt =  N-2:-1:2
                  unt = kron(one,unt);
              end
              Psi= unt*Psi;
          end

        end
    Fidelity_3_state(integer,counter) = abs(Psi'*GHZ);
    end
       integer = integer+1;
end
figure(2)
plot(O_V,Fidelity_3_state(1,:), O_V,Fidelity_3_state(2,:))
grid on;
legend('L=4','L=6')
xlabel('\Omega/V_0')
ylabel('Fidelity F')