% 28.10.2016
%
% 1D spin chain consisting of two different latticed 
% Pulse sequence
% --------------
%                 1st pi/2 pulse on odd atoms -> U(pi/2)
%                 2nd pi pulse on even atoms
%                 3rd U^d(pi/2) (U dagger)
%
clear all;
format long

N  = 7; % total number of atoms


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

        if k_atom == N-1 % for G_{N-1}
            
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
Rho_red = Psi_final*Psi_final';
%-----------------------------------------------------------------------
% reduced density matrix rho_{N-1,N}
Rho_r_Nm1_N = Rho_red;

for cnt = N-1:-1:2
  
    
    psi_J_1 = one;
    psi_J_2 = one;
   
    
    for cnt2 = 1:cnt -1
     
        psi_J_1 = kron(one,psi_J_1);
        psi_J_2 = kron(one,psi_J_2);
                
    end
      psi_J_1 = kron(up,psi_J_1);
        psi_J_2 = kron(do,psi_J_2);
 
      Rho_r_Nm1_N = psi_J_1' * Rho_r_Nm1_N * psi_J_1 + psi_J_2' * Rho_r_Nm1_N * psi_J_2;

end




sigsig = kron(sigma_y,sigma_y);

rho_tilt = sigsig*Rho_r_Nm1_N*sigsig;

EV = sort(eig(Rho_r_Nm1_N*rho_tilt),'descend');

concurrence_Nm1_N = max(0,sqrt(EV(1))-sqrt(EV(2))-sqrt(EV(3))-sqrt(EV(4)))

R = sqrtm(sqrtm(Rho_r_Nm1_N)*rho_tilt*sqrtm(Rho_r_Nm1_N));
EV_R = sort(eig(R),'descend');
Con_Nm1_N  = max(0,EV_R(1)-EV_R(2)-EV_R(3)-EV_R(4));


% %---------------------------------------------------------------------
% % reduced density matrix rho_{1,N}
% Rho_r_1_N = Rho_red;
% 
% for cnt = N-1:-1:2
%   
%     
%     psi_J_1 = kron(up,one);
%     psi_J_2 = kron(do,one);
%     
%     for cnt2 = 1:cnt-1 
%      
%         psi_J_1 = kron(one,psi_J_1);
%         psi_J_2 = kron(one,psi_J_2);
%                 
%     end
%    
%       Rho_r_1_N = psi_J_1' * Rho_r_1_N * psi_J_1 + psi_J_2' * Rho_r_1_N * psi_J_2;
% 
% end
% 
% sigsig = kron(sigma_y,sigma_y);
% 
% rho_tilt = sigsig*Rho_r_1_N*sigsig;
% 
% EV = sort(eig(Rho_r_1_N*rho_tilt),'descend');
% 
% concurrence_1_N = max(0,sqrt(EV(1))-sqrt(EV(2))-sqrt(EV(3))-sqrt(EV(4)));
% 
% R = sqrtm(sqrtm(Rho_r_1_N)*rho_tilt*sqrtm(Rho_r_1_N));
% EV_R = sort(eig(R),'descend');
% Con_1_N  = max(0,EV_R(1)-EV_R(2)-EV_R(3)-EV_R(4));
% 
% 
% %---------------------------------------------------------------------
% % reduced density matrix rho_{N-2,N}
% Rho_r_Nm2_N = Rho_red;
% 
% for cnt = N-1:-1:3
%   
%     
%     psi_J_1 = one;
%     psi_J_2 = one;
%     
%     for cnt2 = 1:cnt-1 
%      
%         psi_J_1 = kron(one,psi_J_1);
%         psi_J_2 = kron(one,psi_J_2);
%                 
%     end
%       psi_J_1 = kron(up,psi_J_1);
%       psi_J_2 = kron(do,psi_J_2);
%       Rho_r_Nm2_N = psi_J_1' * Rho_r_Nm2_N * psi_J_1 + psi_J_2' * Rho_r_Nm2_N * psi_J_2;
% 
% end
% 
% psi_J_1 = kron(up,one);
% psi_J_1 = kron(one,psi_J_1);
% psi_J_2 = kron(do,one);
% psi_J_2 = kron(one,psi_J_2);
% Rho_r_Nm2_N = psi_J_1' * Rho_r_Nm2_N * psi_J_1 + psi_J_2' * Rho_r_Nm2_N * psi_J_2;
% 
% 
% sigsig = kron(sigma_y,sigma_y);
% 
% rho_tilt = sigsig*Rho_r_Nm2_N*sigsig;
% 
% EV = sort(eig(Rho_r_Nm2_N*rho_tilt),'descend');
% 
% concurrence_Nm2_N = max(0,sqrt(EV(1))-sqrt(EV(2))-sqrt(EV(3))-sqrt(EV(4)));
% 
% R = sqrtm(sqrtm(Rho_r_Nm2_N)*rho_tilt*sqrtm(Rho_r_Nm2_N));
% EV_R = sort(eig(R),'descend');
% Con_Nm2_N  = max(0,EV_R(1)-EV_R(2)-EV_R(3)-EV_R(4));
% 
% %---------------------------------------------------------------------
% % reduced density matrix rho_{2,N}
% Rho_r_2_N = Rho_red;
% 
% for cnt = N-1:-1:3
%   
%     
%     psi_J_1 = kron(up,one);
%     psi_J_2 = kron(do,one);
%     
%     for cnt2 = 1:cnt-1
%      
%         psi_J_1 = kron(one,psi_J_1);
%         psi_J_2 = kron(one,psi_J_2);
%                 
%     end
% 
%       Rho_r_2_N = psi_J_1' * Rho_r_2_N * psi_J_1 + psi_J_2' * Rho_r_2_N * psi_J_2;
% 
% end
% 
% psi_J_1 = kron(one,one);
% psi_J_1 = kron(up,psi_J_1);
% psi_J_2 = kron(one,one);
% psi_J_2 = kron(do,psi_J_2);
% 
% Rho_r_2_N = psi_J_1' * Rho_r_2_N * psi_J_1 + psi_J_2' * Rho_r_2_N * psi_J_2;
% 
% sigsig = kron(sigma_y,sigma_y);
% 
% rho_tilt = sigsig*Rho_r_2_N*sigsig;
% 
% EV = sort(eig(Rho_r_2_N*rho_tilt),'descend');
% 
% concurrence_2_N = max(0,sqrt(EV(1))-sqrt(EV(2))-sqrt(EV(3))-sqrt(EV(4)));
% 
% R = sqrtm(sqrtm(Rho_r_2_N)*rho_tilt*sqrtm(Rho_r_2_N));
% EV_R = sort(eig(R),'descend');
% Con_2_N  = max(0,EV_R(1)-EV_R(2)-EV_R(3)-EV_R(4));
% 
% 
% 
% 
%  %---------------------------------------------------------------------
% if N >5
% 
%     % reduced density matrix rho_{Nm3,N}
%     Rho_r_Nm3_N = Rho_red;
% 
%     for cnt = N-1:-1:4
% 
% 
%         psi_J_1 = one;
%         psi_J_2 = one;
% 
%         for cnt2 = 1:cnt-1
% 
%             psi_J_1 = kron(one,psi_J_1);
%             psi_J_2 = kron(one,psi_J_2);
% 
%         end
% 
%            psi_J_1 = kron(up,psi_J_1);
%            psi_J_2 = kron(do,psi_J_2);
%            Rho_r_Nm3_N = psi_J_1' * Rho_r_Nm3_N * psi_J_1 + psi_J_2' * Rho_r_Nm3_N * psi_J_2;
% 
%     end
%     
% psi_J_1 = kron(one,one);
% psi_J_1 = kron(up,psi_J_1);
% psi_J_1 = kron(one,psi_J_1);
% psi_J_2 = kron(one,one);
% psi_J_2 = kron(do,psi_J_2);
% psi_J_2 = kron(one,psi_J_2);
% Rho_r_Nm3_N = psi_J_1' * Rho_r_Nm3_N * psi_J_1 + psi_J_2' * Rho_r_Nm3_N * psi_J_2;
% 
% 
% psi_J_1 = kron(up,one);
% psi_J_1 = kron(one,psi_J_1);
% 
% psi_J_2 = kron(do,one);
% psi_J_2 = kron(one,psi_J_2);
% Rho_r_Nm3_N = psi_J_1' * Rho_r_Nm3_N * psi_J_1 + psi_J_2' * Rho_r_Nm3_N * psi_J_2;
% 
% 
% sigsig = kron(sigma_y,sigma_y);
% 
% rho_tilt = sigsig*Rho_r_Nm3_N*sigsig;
% 
% EV = sort(eig(Rho_r_Nm3_N*rho_tilt),'descend');
% 
% concurrence_Nm3_N = max(0,sqrt(EV(1))-sqrt(EV(2))-sqrt(EV(3))-sqrt(EV(4)));
% 
% R = sqrtm(sqrtm(Rho_r_Nm3_N)*rho_tilt*sqrtm(Rho_r_Nm3_N));
% EV_R = sort(eig(R),'descend');
% Con_Nm3_N  = max(0,EV_R(1)-EV_R(2)-EV_R(3)-EV_R(4));
% 
% end
% 
% 
% 
% if N == 7
%         % reduced density matrix rho_{2,4}
%     Rho_r_2_4 = Rho_red;
% 
%     for cnt = N-1:-1:4
% 
% 
%         psi_J_1 = kron(one,up);
%         psi_J_2 = kron(one,do);
% 
%         for cnt2 = 1:cnt-1
% 
%             psi_J_1 = kron(one,psi_J_1);
%             psi_J_2 = kron(one,psi_J_2);
% 
%         end
% 
%        
%            Rho_r_2_4 = psi_J_1' * Rho_r_2_4 * psi_J_1 + psi_J_2' * Rho_r_2_4 * psi_J_2;
% 
%     end
%     
%     
%    psi_J_1 = kron(up,one);
%    psi_J_2 = kron(do,one);
%    psi_J_1 = kron(one,psi_J_1);
%    psi_J_2 = kron(one,psi_J_2); 
%    psi_J_1 = kron(one,psi_J_1);
%    psi_J_2 = kron(one,psi_J_2); 
%     Rho_r_2_4 = psi_J_1' * Rho_r_2_4 * psi_J_1 + psi_J_2' * Rho_r_2_4 * psi_J_2;
%     
% 
%    psi_J_1 = kron(one,one);
%    psi_J_2 = kron(one,one); 
%    psi_J_1 = kron(up,psi_J_1);
%    psi_J_2 = kron(do,psi_J_2); 
%    
%    Rho_r_2_4 = psi_J_1' * Rho_r_2_4 * psi_J_1 + psi_J_2' * Rho_r_2_4 * psi_J_2;
%    
%    sigsig = kron(sigma_y,sigma_y);
% 
%    rho_tilt = sigsig*Rho_r_2_4*sigsig;
% 
%    EV = sort(eig(Rho_r_2_4*rho_tilt),'descend');
% 
%    concurrence_2_4 = max(0,sqrt(EV(1))-sqrt(EV(2))-sqrt(EV(3))-sqrt(EV(4)));
% 
%     R = sqrtm(sqrtm(Rho_r_2_4)*rho_tilt*sqrtm(Rho_r_2_4));
%     EV_R = sort(eig(R),'descend');
%     Con_2_4  = max(0,EV_R(1)-EV_R(2)-EV_R(3)-EV_R(4))
% 
% end
% 
% 
% 
% 
% % ANALYTICAL RESULTS FOR THE CASE N = 5
% %----------------------------------------------------------------------------
% % N = 5: 
% % Anaytical results for the reduced density matrix rho_15
% % 
% % Rho_r_1_N = kron(do,do)*kron(do,do)' -1/(2*sqrt(2))*kron(do,do)*kron(do,plus)'...
% %     -1/(2*sqrt(2))*kron(do,do)*kron(plus,do)'+1/4 * kron(do,do)*kron(plus,plus)'...
% %     -1/(2*sqrt(2))*kron(do,plus)*kron(do,do)'+1/2 * kron(do,plus) * kron(do,plus)'...
% %     +1/4* kron(do,plus)*kron(plus,do)' -1/(2*sqrt(2))*kron(do,plus)*kron(plus,plus)'...
% %     +1/4*kron(plus,do)*kron(do,plus)'-1/(2*sqrt(2))*kron(plus,do)*kron(plus,plus)'...
% %     +1/4* kron(plus,plus)*kron(do,do)' - 1/(2*sqrt(2))*kron(plus,plus)*kron(do,plus)'...
% %     - 1/(2*sqrt(2))*kron(plus,plus)*kron(plus,do)' +1/2 * kron(plus,plus)*kron(plus,plus)'...
% %     -1/(2*sqrt(2))*kron(plus,do)*kron(do,do)'+1/2 * kron(plus,do) *kron(plus,do)';
% %     
% % sigsig = kron(sigma_y,sigma_y);
% % 
% % rho_tilt = sigsig*Rho_r_1_N*sigsig;
% % 
% % EV = sort((eig(Rho_r_1_N*rho_tilt)),'descend');
% % 
% % concurrence_1_N = max(0,sqrt(EV(1))-sqrt(EV(2))-sqrt(EV(3))-sqrt(EV(4)))
% % 
% % R = sqrtm(sqrtm(Rho_r_1_N)*rho_tilt*sqrtm(Rho_r_1_N));
% % EV_R = sort(eig(R),'descend');
% % Con_1_N  = max(0,EV_R(1)-EV_R(2)-EV_R(3)-EV_R(4));
% 
% 
% %----------------------------------------------------------------------------
% % N = 5: 
% % Anaytical results for the reduced density matrix rho_35
% % 
% %  Rho_r_3_5 = kron(do,do)*kron(do,do)' -1/2 * kron(do,do)*kron(plus,plus)'...
% %  -1/(2*sqrt(2)) * kron(do,do) *kron(plus,do)' +1/4 * kron(do,do)*kron(plus,plus)'...
% %  -1/2 * kron(plus,plus)*kron(do,do)' +1/(4*sqrt(2)) * kron(plus,plus)*kron(plus,do)'...
% %  -1/(2*sqrt(2))* kron(plus,do)*kron(do,do)'+ 1/(4*sqrt(2))*kron(plus,do)*kron(plus,plus)'...
% %  +1/2 * kron(plus,do)*kron(plus,do)' - 1/(2*sqrt(2))* kron(plus,do)*kron(plus,plus)' ...
% %  +1/(4) * kron(plus,plus)*kron(do,do)' - 1/(2*sqrt(2))*kron(plus,plus)*kron(plus,do)'...
% %  + 1/2 * kron(plus,plus)* kron(plus,plus)';
% %  
% % sigsig = kron(sigma_y,sigma_y);
% % 
% % rho_tilt = sigsig*Rho_r_3_5*sigsig;
% % 
% % EV = sort(eig(Rho_r_3_5*rho_tilt),'descend');
% % 
% % concurrence_3_5 = max(0,sqrt(EV(1))-sqrt(EV(2))-sqrt(EV(3))-sqrt(EV(4)));
% % 
% % R = sqrtm(sqrtm(Rho_r_3_5)*rho_tilt*sqrtm(Rho_r_3_5));
% % EV_R = sort(eig(R),'descend');
% % Con_3_5  = max(0,EV_R(1)-EV_R(2)-EV_R(3)-EV_R(4));
% 
% 
% %----------------------------------------------------------------------------
% % N = 5: 
% % Anaytical results for the reduced density matrix rho_25
% 
% %   Rho_r_2_5  = kron(do,do)*kron(do,do)' - 1/(2*sqrt(2)) *kron(do,do)*kron(do,plus)'...
% %   -1/(2*sqrt(2))*kron(do,do)*kron(plus,do)' + 1/4 *kron(do,do)*kron(plus,plus)'...
% %   -1/(2*sqrt(2)) * kron(do,plus)*kron(do,do)' + 1/2 *kron(do,plus)*kron(do,plus)'...
% %   +1/4 * kron(do,plus)*kron(plus,do)' - 1/(2*sqrt(2))* kron(do,plus)*kron(plus,plus)'...
% %   -1/(2*sqrt(2))*kron(plus,do)*kron(do,do)' +1/4* kron(plus,do)*kron(do,plus)' ...
% %   + 1/2* kron(plus,do)*kron(plus,do)' - 1/(2*sqrt(2)) *kron(plus,do)*kron(plus,plus)'...
% %   +1/4 *kron(plus,plus)*kron(do,do)' -1/(2*sqrt(2)) * kron(plus,plus)*kron(do,plus)'...
% %   - 1/(2*sqrt(2)) *kron(plus,plus)*kron(plus,do)' + 1/2 * kron(plus,plus)*kron(plus,plus)';
% % 
% % sigsig = kron(sigma_y,sigma_y);
% % 
% % rho_tilt = sigsig*Rho_r_2_5*sigsig;
% % 
% % EV = sort(eig(Rho_r_2_5*rho_tilt),'descend');
% % 
% % concurrence_2_5 = max(0,sqrt(EV(1))-sqrt(EV(2))-sqrt(EV(3))-sqrt(EV(4)));
% % 
% % R = sqrtm(sqrtm(Rho_r_2_5)*rho_tilt*sqrtm(Rho_r_2_5));
% % EV_R = sort(eig(R),'descend');
% % Con_2_5  = max(0,EV_R(1)-EV_R(2)-EV_R(3)-EV_R(4));
% 
% 
% %----------------------------------------------------------------------------
% % N = 5: 
% % Anaytical results for the reduced density matrix rho_24
% 
% % Rho_r_24 = kron(do,do)*kron(do,do)' - 1/(2*sqrt(2))*kron(do,do)*kron(do,plus)'...
% %     -1/(2*sqrt(2)) * kron(do,do)*kron(plus,do)' +1/4 * kron(do,do)*kron(plus,plus)'...
% %     -1/(2*sqrt(2))* kron(do,plus)*kron(do,do)' + 1/2 * kron(do,plus)*kron(do,plus)'...
% %     +1/4 * kron(do,plus)*kron(plus,do)' -1/(2*sqrt(2))*kron(do,plus)*kron(plus,plus)'...
% %     -1/(2*sqrt(2)) * kron(plus,do)*kron(do,do)'+ 1/4 * kron(plus,do)*kron(do,plus)'...
% %     +1/2*kron(plus,do)*kron(plus,do)' - 1/(2*sqrt(2))*kron(plus,do)*kron(plus,plus)'...
% %     +1/4 *kron(plus,plus)*kron(do,do)' - 1/(2*sqrt(2))*kron(plus,plus)*kron(do,plus)'...
% %     -1/(2*sqrt(2)) *kron(plus,plus)*kron(plus,do)' +1/2*kron(plus,plus)*kron(plus,plus)';
% %     
% % sigsig = kron(sigma_y,sigma_y);
% % 
% % rho_tilt = sigsig*Rho_r_24*sigsig;
% % 
% % EV = sort(eig(Rho_r_24*rho_tilt),'descend');
% % 
% % concurrence_24 = max(0,sqrt(EV(1))-sqrt(EV(2))-sqrt(EV(3))-sqrt(EV(4)));
% % 
% % R = sqrtm(sqrtm(Rho_r_24)*rho_tilt*sqrtm(Rho_r_24));
% % EV_R = sort(eig(R),'descend');
% % Con_24  = max(0,EV_R(1)-EV_R(2)-EV_R(3)-EV_R(4));
