%14.11.16
%
% 1D spin chain consisting of two different latticed 
% Pulse sequence
% --------------
%                 1st pi/2 pulse on odd atoms -> U(pi/2)
%                 2nd pi pulse on even atoms
%                 3rd U^d(pi/2) (U dagger)
%
% Identification of the final state for N = 3 
%
clear all;
format long

N  = 3; % total number of atoms


up = [1;0]; % spin up
do = [0;1]; % spin down

plus = 1/sqrt(2) * [up+do];
minus = 1/sqrt(2) * [-up+do];

sigma_y = [0, -i; i 0];
num     = [1,0; 0 ,0];
one     = eye(2);
hadamard = 1/sqrt(2) * [-1,1;1,1];


init = do;


for cnt = 1:N-1
      init   = sparse(kron(do,init)); % initial state vector
end


       Hada_N =  hadamard;
       for cnt = 1: N-1
             if mod(cnt,2) == 0 % cnt ungerade
                 
                Hada_N = sparse(kron(hadamard,Hada_N)); 
            
             else
               
                 Hada_N = sparse(kron(one,Hada_N));
             end
       end

Psi_final = Hada_N*init;



sig_2 = kron ([1,1;-1,1],one);
sig_2 = sparse(kron(one,sig_2));

num_1 = kron(one,one);
num_1 = sparse(kron(num,num_1));


num_3 = kron(one,num);
num_3 = sparse(kron(one,num_3));

one_N = kron(one,one);
one_N = sparse(kron(one_N,one));


G = (one_N - sig_2 * (one_N - num_1)*(one_N-num_3)); % acting on atom j-1,j,j+1
if N > 3
    for k_atom = 2:2:N-1
       
        G_N = one; %reset G_N

        if k_atom == N-1 % for G_{N-1}untitled
            
           G_N = sparse(kron(G_N,G));
           
           for cnt = 1:k_atom -3
               G_N = sparse(kron(one,G_N));
           end
         
           
        else % all other G_j
            
            for cnt = N-3:-1:1
                
                if cnt == k_atom -1
                  
                    G_N = sparse(kron(G,G_N));
                    continue
                end
                G_N = sparse(kron(one,G_N));  
            end
         
            
        end

        Psi_final = G_N * Psi_final;  
       
    end
else
    G_N = G;
    Psi_final = G * Psi_final;
end



Psi_final = Hada_N*Psi_final;

%Calculation of the density matrix
Rho = Psi_final*Psi_final';
%-----------------------------------------------------------------------
% reduced density matrix rho_{N-1,N}
Rho_A = Rho;

psi_J_1 = kron(one,up);
psi_J_1 = kron(one,psi_J_1);
psi_J_2 = kron(one,do);
psi_J_2 = kron(one,psi_J_2);

Rho_A = psi_J_1' * Rho_A * psi_J_1 + psi_J_2' * Rho_A * psi_J_2;

psi_J_1 = kron(one,up);
psi_J_2 = kron(one,do);
Rho_A = psi_J_1' * Rho_A * psi_J_1 + psi_J_2' * Rho_A * psi_J_2;

S_A = det(Rho_A)


Rho_B = Rho;
psi_J_1 = kron(one,one);
psi_J_1 = kron(up,psi_J_1);
psi_J_2 = kron(one,one);
psi_J_2 = kron(do,psi_J_2);

Rho_B = psi_J_1' * Rho_B * psi_J_1 + psi_J_2' * Rho_B * psi_J_2;

psi_J_1 = kron(one,up);
psi_J_2 = kron(one,do);
Rho_B = psi_J_1' * Rho_B * psi_J_1 + psi_J_2' * Rho_B * psi_J_2;

S_B = det(Rho_B)

Rho_C = Rho;
psi_J_1 = kron(one,one);
psi_J_1 = kron(up,psi_J_1);
psi_J_2 = kron(one,one);
psi_J_2 = kron(do,psi_J_2);

Rho_C = psi_J_1' * Rho_C * psi_J_1 + psi_J_2' * Rho_C * psi_J_2;

psi_J_1 = kron(up,one);
psi_J_2 = kron(do,one);
Rho_C = psi_J_1' * Rho_C * psi_J_1 + psi_J_2' * Rho_C * psi_J_2;

S_C = det(Rho_C)




