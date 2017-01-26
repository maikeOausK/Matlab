% Time Propagation of a chain of two spins

% Hamiltonian of the system
% H = sum_k Omega_k * sigma_x^(k) + sum_k Delta_k * n_k + sum_k,m V_{km} * n_k * n_m

% Omega_k : Rabi frequency of particle k
% Delta_k : Detuning of the rydberg state of particle k
% V_{km}  : interaction potential

% 1st approximation: use random numbers for the values

%
% Hamiltonian for 2 spins
% H = Omega * [ sigma_x^(1)+sigma_x^(2) ] + Delta * [ n_1 + n_2 ] + V_{12}
% * n_1 * n_2
% 
clear all;
format long;
Omega_1 = rand
Omega_2 = Omega_1;
Delta_1 = 0;%rand;
Delta_2 = Delta_1;
V_12  =10000*Omega_1
Period = 1/(sqrt(2)*Omega_1)

% Definition of the operators
sigma_x = [0 1; 1 0]; % Pauli matrix in x direction
one     = [1 0; 0 1]; % Identity
num     = [1 0; 0 0]; % Number operator

% Operators in 4 dimensional space (Product space)
% kron -> tensor product

sigma_x_1 = kron(sigma_x, one); % sigma_x for the 1st particle 
sigma_x_2 = kron(one, sigma_x); % sigma_x for the 2nd particle

num_1 = kron(num, one); % number operator for 1st particle
num_2 = kron(one, num); % number oerator for 2nd particle
num_12 = num_1+num_2;

part_1 = [0; 1]; % initial vector of 1st particle (spin down)
part_2 = [0; 1]; % initital vector of 2nd particle (spin down)

init_prod_state = kron(part_1, part_2); % initial product state for 2 spin down particle


Hamil = zeros(length(num)*length(num)); % initialisation of the hamiltonian matrix

Hamil = Omega_1 * sigma_x_1 + Omega_2 * sigma_x_2 + Delta_1 * num_1 +...
    Delta_2 * num_2 + V_12 * num_1*num_2 ;

%_______________________________________________________
% FINAL STATE VIA SOLVING THE TDSE (with ode solver)
%Psi = init_prod_state
%dt_Psi = HPsi(Psi, Hamil,num_1)
tspan = [0:0.01:10];
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
 [t,Psi] = ode45(@(t,Psi)HPsi(Psi,Hamil), tspan, init_prod_state,opts);

 Expect_num_1 = zeros(1,length(Psi));
  Expect_num_2 = zeros(1,length(Psi));
  Expect_num_12 = zeros(1,length(Psi));
  
 for cnt = 1: length(Psi)

     Psi(cnt,:) =Psi(cnt,:)/ sqrt(Psi(cnt,:)*Psi(cnt,:)');
   
     Expect_num_1(cnt) = Psi(cnt,:)* (num_1 * Psi(cnt,:)');
     Expect_num_2(cnt) = Psi(cnt,:)* (num_2 * Psi(cnt,:)');
 end

%  figure(1)
%  plot(t,Expect_num_1,'*',t,Expect_num_2,'o')
%  legend('Spin 1','Spin 2')
%  xlabel('Time t')
% ylabel('<n_{1/2}>')
%  
 
%_______________________________________________________
% FINAL STATE -> STATIONARY STATE EXPANSION
% |final(t)> = exp(-iHt) |initial>
%  

 

 for cnt = 1: length(tspan)
 final(cnt,:) =   expm(-Hamil*1i*tspan(cnt)) * init_prod_state;
 
     Expect_num_1(cnt) = final(cnt,:)* (num_1 * final(cnt,:)');
     Expect_num_2(cnt) = final(cnt,:)* (num_2 * final(cnt,:)');
     Expect_num_12(cnt) = final(cnt,:)* (num_12 * final(cnt,:)');
 end
 
 
 
 
%figure(2)
% plot(t,Expect_num_1,'-*',t,Expect_num_2,'-o')
% legend('Spin 1','Spin 2')
% xlabel('Time t')
% ylabel('<n_{1/2}>')
%  
figure(3)
plot(t,Expect_num_12,'-r*',t,sin(sqrt(2)*Omega_1*t).^2,'-b')
xlabel('Time t')
ylabel('<n_{1+2}>')
 
%____________________________________________________
% % Omega = 1
% % Delta = 0
% % V     = 0
% Hamil= sigma_x_1 + sigma_x_2 ;
%  tspan = [0:0.01:5];
% opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
%  [t,Psi] = ode45(@(t,Psi)HPsi(Psi,Hamil), tspan, init_prod_state,opts);
% 
%  Expect_num_1 = zeros(1,length(Psi));
%   Expect_num_2 = zeros(1,length(Psi));
%  
%  for cnt = 1: length(Psi)
% 
%      Psi(cnt,:) =Psi(cnt,:)/ sqrt(Psi(cnt,:)*Psi(cnt,:)');
%    
%      Expect_num_1(cnt) = Psi(cnt,:)* (num_1 * Psi(cnt,:)');
%      Expect_num_2(cnt) = Psi(cnt,:)* (num_2 * Psi(cnt,:)');
%  end
% 
% %  figure(3)
% %  plot(t,Expect_num_1,'-*',t,Expect_num_2,'-o')
% %  legend('Spin 1','Spin 2')
% %  xlabel('Time t')
% % ylabel('<n_{1/2}>')
% %  
% [V,D]=eig(-i*sigma_x_1);
% [V2,D2]=eig(-i*sigma_x_1);
% 
% Psi_x1 = zeros(4,1);
% Psi_x2 = zeros(4,1);
% 
% v_1 = [0;1;0;1];
% v_2 = [0;1;0;-1];
% 
% v_3 = [0;0;1;-1];
% v_4 = [0;0;1;1];
% 
% Psi_x1= 1/2*v_1*exp(-tspan*i) - 1/2*v_2*exp(tspan*i);
% Psi_x2= 1/2*v_4*exp(-tspan*i) - 1/2*v_3*exp(tspan*i);
% Psi_s = Psi_x1+Psi_x2;
% 
%  for cnt = 1: length(Psi_s)
% 
%    
%      Expect_num_1_s(cnt) = Psi_s(:,cnt)'* (num_1 * Psi_s(:,cnt));
%      Expect_num_2_s(cnt) = Psi_s(:,cnt)'* (num_2 * Psi_s(:,cnt));
%  end
% 
%  figure(3)
%   plot(t,Expect_num_1,'-m*',t,Expect_num_2,'xk','MarkerSize',12)
%   hold on;
%  plot(tspan',Expect_num_1_s,'rd',t,Expect_num_2_s,'b+')
% plot(t,(sin(t)).^2,'sg','MarkerSize',10)
%  legend('Spin 1 dgl','Spin 2 dgl','Spin 1 man','Spin 2 man','sin^2')
%  xlabel('Time t')
% ylabel('<n_{1/2}>')
