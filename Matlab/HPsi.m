function dt_Psi = HPsi(Psi, Hamil)
dt_Psi = zeros(4,1);
%  
dt_Psi = -1i* Hamil * Psi;
end

