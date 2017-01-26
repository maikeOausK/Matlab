% 
clear all;
format long
up = [1;0];
down = [0;1];
plus = 1/sqrt(2)*[1;1];
minus = 1/sqrt(2)*[1;-1];
num = [1,0;0,0];
hadamard = 1/sqrt(2) * [-1,1;1,1];
one = eye(2);

dd = kron(down,down);
pp = kron(plus,plus);
dp = kron(down,plus);
pd = kron(plus,down);

rho_12 = dd*dd'-1/2*((dd*pp')+(pp*dd'))+ 1/2*(pp*pp');

sigma_y = [0 -i; i 0 ];
sigma_x = [0 1; 1 0 ];

sigsig = kron(sigma_y,sigma_y);

rho_tilt = sigsig*rho_12*sigsig;

EV = (eig(rho_12*rho_tilt));

concurrence = sqrt(EV(1))-sqrt(EV(2))-sqrt(EV(3))-sqrt(EV(4));

R = sqrtm(sqrtm(rho_12)*rho_tilt*sqrtm(rho_12));
EV_R = eig(R);
Con  = EV_R(1)-EV_R(2)-EV_R(3)-EV_R(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rho_12 = dd*dd'+ dp*dp'+pp*pp'...
%           -1/sqrt(2)*(dd*dp'+dp*dd'+dp*pp'+pp*dp');
%       
% 
% 
% rho_tilt = sigsig*rho_12*sigsig;
% 
% EV = (eig(rho_12*rho_tilt));
% 
% concurrence = sqrt(EV(1))-sqrt(EV(2))-sqrt(EV(3))-sqrt(EV(4));


% R = sqrtm(sqrtm(rho_12)*rho_tilt*sqrtm(rho_12));
% EV_R = eig(R);
% Con  = EV_R(1)-EV_R(2)-EV_R(3)-EV_R(4);

