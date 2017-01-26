%
% 18.01.2017
%-------------
% Teleportation of a single qubit state 
% a |0> + b|1> over an array of N qubits 
% by just applying a sequence of pi pulses
% interaction is restricted to nearest neighbor interaction

clear all;
format short;



N   = 5; % total number of atoms


% Single particle operators and states

v0  = [0;1];
v1  = [1;0];



sigma_x = [0,1;1,0];
sigma_y = [0,i;-i,0];
sigma_z = [-1,0;0,1];


b = 1;
a = 1;
Psi = [b;a];
Psi = Psi/norm(Psi);

num     = [1,0; 0 ,0];
one     = eye(2);

% 3 particle operators
sig_2 = kron (sigma_y,one); 
sig_2 = sparse(kron(one,sig_2)); % 3 particle sigma operator 
                                 % 1 x sigma x 1
                                 
num_1 = kron(one,one);
num_1 = sparse(kron(num,num_1));


num_3 = kron(one,num);
num_3 = sparse(kron(one,num_3));

one_N = kron(one,one);
one_N = sparse(kron(one_N,one));

init = v0;
sigma_x_N = kron(one,sigma_x);

for cnt = 1:N-2
     init  = sparse((kron(init,v0))); % initial state vector
     sigma_x_N = kron(one,sigma_x_N);
end
init = sparse(kron(Psi,init));




G_c2 = sparse(kron(sigma_y,one) * (eye(4)- kron(one,num)));

G_c1 = sparse(kron(one,sigma_y) * (eye(4)- kron(num,one)));

G_3 = sparse( sig_2 * (one_N - num_1)*(one_N-num_3)); % acting on atom j-1,j,j+1



if N > 3

 
    for cnt =1:N-1 % determination of Toffoli gates
        Gate = one;
        G_2 = one;
        
     
         
        if cnt == N-1  % last Toffoli
         
            Gate = G_3;
            G_2 = G_c1;
            for cnt2 = 1:N-3
                Gate = sparse(kron(one,Gate));
           
                G_2 = sparse(kron(one,G_2));
 
            end
            
            G_2 = sparse(kron(one,G_2));
    
            %Seq =Seq+GateN+G_2;
            Psif = sparse(expm(-i*(pi/2)*G_2)*Psif);
            Psif = sparse(expm(-i*(pi/2)*GateN)*Psif);
%             
% 
%                 
       elseif cnt == 1  % 1st Toffoli (correct)
              
            for cnt2 = 1:N-4
                Gate = sparse( kron(Gate,one));
                G_2 = sparse(kron(G_2,one));
            end
            Gate = sparse(kron(G_3,Gate));
  
            G_2 = sparse(kron(G_2,one));
            G_2 = sparse(kron(G_c2,G_2));

            
            Psif = sparse(expm(-i*(pi/2)*Gate)*init);
            Psif = sparse(expm(-i*(pi/2)*G_2)*Psif);
            GateN = Gate;
           
        
        elseif cnt == N-2
            
            Gate = G_3;
            for cnt2 = 1:N-3
                Gate = sparse(kron(one,Gate));
            end
             
             Psif = sparse(expm(-i*(pi/2)*Gate)*Psif);
             Psif = sparse(expm(-i*(pi/2)*GateN)*Psif);
             GateN = Gate;
     
        else
                   
            
         for cnt2 = N-3:-1:1 %needed for N>4
   
             if cnt2 == cnt  %FILL IN THE MISSING PART
            

                 Gate = sparse(kron(G_3,Gate));
            
                 continue
             end
             Gate = sparse(kron(one,Gate));
           
              
         end
             Psif = sparse(expm(-i*(pi/2)*Gate)*Psif);
             Psif = sparse(expm(-i*(pi/2)*GateN)*Psif);
    
         GateN = Gate;

        

        
        end
    
    end
   
 end

M00 = one;
for cnt = 1:N-1
    M00 = sparse(kron(v0,M00));
end

Psif = (M00' *Psif);

((sigma_x)^(N-1))*Psif;


