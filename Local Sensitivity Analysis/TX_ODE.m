function dx= TX_Model(t,x,p)
global par 
dx = zeros(5,1);


%%%% Molecule Species

X = x(1);  
Y = x(2);  
Z = x(3);
Pz_rep = x (4);
GFP = x(5);

Px = par.P_x;
Py = par.P_y;
Pz = par.P_z - Pz_rep;

Py_active = p(10)*X*Py;
Pz_active = p(11)*X*Pz;


dx(1) = p(1)*Px - p(11)*X*Pz - p(10)*X*Py - p(4)*X;

dx(2) = p(2)*Py_active - p(9)*Y*Pz - p(9)*Y*Pz_active - p(5)*Y; 

dx(3) = p(3)*Pz_active - p(6)*Z;

dx(4) = p(9)*Y*Pz + p(9)*Y*Pz_active;

dx(5) = p(7)*Z - p(8)*GFP; 





