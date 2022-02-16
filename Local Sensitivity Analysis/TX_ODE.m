function dx= TX_Model(t,x,p)
global par 
dx = zeros(5,1);

par.P_x = 1e-9;
par.P_y = 1e-9;
par.P_z = 1e-9;

%%%% Molecule Species

X = x(1);  
Y = x(2);  
Z = x(3);
Pz_rep = x (4);
GFP = x(5);
Py_active = x(6);
Pz_active = x(7);

Px = par.P_x;
Py = par.P_y - Py_active;
Pz = par.P_z - Pz_active - Pz_rep;

% Py_active = p(10)*X*Py;
% Pz_active = p(11)*X*Pz;


dx(1) = p(1)*Px - p(11)*X*Pz - p(10)*X*Py - p(4)*X;

dx(2) = p(2)*Py_active - p(9)*Y*Pz - p(9)*Y*Pz_active - p(5)*Y; 

dx(3) = p(3)*Pz_active - p(6)*Z;

dx(4) = p(9)*Y*Pz + p(9)*Y*Pz_active;

dx(5) = p(7)*Z - p(8)*GFP; 

dx(6) = p(10)*X*Py - p(12)*Py_active;

dx(7) = p(11)*X*Pz - p(9)*Pz_active*Y - p(12)*Pz_active;





