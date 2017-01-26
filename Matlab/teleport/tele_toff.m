%
% 18.01.2017
%-------------
% Teleportation of a single qubit state 
% a |0> + b|1> over an array of N qubits 
% by just applying a sequence of pi pulses

clear all;
format long;



N   = 15; % total number of atoms
% auxiliary variable to determine nr of applied toffolis

% Single particle operators and states

v0  = [0;1];
v1  = [1;0];



sigma_x = [0,1;1,0];
sigma_y = [0,i;-i,0];
sigma_1 = [-1,0;0,1];

Hada    = [-1,1;1,1]/sqrt(2);

b = 1;
a = 1;
Psi = [b;a];
Psi = Psi/norm(Psi);

num     = [1,0; 0 ,0];
one     = eye(2);




init = v0;
for cnt = 1:N-2
     init  = sparse((kron(init,v0))); % initial state vector
end
init = sparse(kron(Psi,init));



% TOFFOLI GATE 
% 1st and 3rd qubit -> control qubits
% 2nd qubit target qubit
toff = [1,0,0,0,0,0,0,0; 0,1,0,0,0,0,0,0; 0,0,1,0,0,0,0,0; 0,0,0,1,0,0,0,0; 0,0,0,0,1,0,0,0;...
        0,0,0,0,0,0,0,1; 0,0,0,0,0,0,1,0; 0,0,0,0,0,1,0,0]; 
    
% CNOT GATE    
CNOTc2 = [1,0,0,0; 0,0,0,1; 0,0,1,0; 0,1,0,0]; % 2nd qubit control qubit
    
CNOTc1 = [1,0,0,0; 0,1,0,0; 0,0,0,1; 0,0,1,0]; % 1st qubit control qubit


GateN = 0;
Gate1 = 0;

if N > 3
    
    for cnt =1:N-1 % determination of Toffoli gates
        Gate = one;
        CNOT = one;
        
     
         
        if cnt == N-1  % last Toffoli
         
            Gate = toff;
            CNOT = CNOTc1;
            for cnt2 = 1:N-3
                Gate = sparse(kron(one,Gate));
           
                CNOT = sparse(kron(one,CNOT));
 
            end
            
            CNOT = sparse(kron(one,CNOT));
            Seq =GateN*CNOT*Seq;
          
%                     
%             
% 
%                 
        elseif cnt == 1  % 1st Toffoli (correct)
              
            for cnt2 = 1:N-4
                Gate = sparse( kron(Gate,one));
                CNOT = sparse(kron(CNOT,one));
            end
            Gate = sparse(kron(toff,Gate));
  
            CNOT = sparse(kron(CNOT,one));
            CNOT = sparse(kron(CNOTc2,CNOT));
            Seq = CNOT*Gate;
            GateN = Gate;
           
%         
        elseif cnt == N-2
            
            Gate = toff;
            for cnt2 = 1:N-3
                Gate = sparse(kron(one,Gate));
            end
             Seq = GateN*Gate*Seq;
             
             GateN = Gate;
     
        else
                   
            
         for cnt2 = N-3:-1:1 %needed for N>4
      
             if cnt2 == cnt  %FILL IN THE MISSING PART
            

                 Gate = sparse(kron(toff,Gate));
            
                 continue
             end
             Gate = sparse(kron(one,Gate));
           
              
         end
          
         Seq = GateN*Gate*Seq;
         GateN = Gate;

        

        
         end
    
    end
end

Psif = (Seq*init);

M00 = one;
for cnt = 1:N-1
    M00 = sparse(kron(v0,M00));
end

Psif = (M00' *Psif);

((sigma_x)^(N-1))*Psif;

