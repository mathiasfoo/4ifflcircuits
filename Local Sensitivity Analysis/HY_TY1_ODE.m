function dx= HY_TY1_ODE(t,x,p)
global par 

par.Px = 1e-9;
par.Py = 1e-9;
par.Pz = 1e-9;

dx = zeros(7,1);


%%%%molecule species

X = x(1);  
Y = x(2);  
Z = x (3);
GFP = x(4);
Py_active = x(5);
Pz_active = x(6);
Y_prime = x(7);


Px = par.P_x; 
Py = par.P_y - Py_active; 
Pz = par.P_z - Pz_active; 

% Py_active = p(9)*X*Py; 
% Pz_active = p(10)*X*Pz;

dx(1) = p(1)*Px - p(10)*X*Pz- p(9)*X*Py - p(4)*X;

dx(2) = p(2)*Py_active - p(5)*Y;

dx(3) = p(3)*Pz_active - p(11)*Y_prime*Z - p(6)*Z;

dx(4) = p(7)*Z - p(8)*GFP;

dx(5) = p(9)*X*Py - p(12)*Py_active;

dx(6) = p(10)*X*Pz - p(12)*Pz_active;

dx(7) = p(13)*Y - p(11)*Z*Y_prime - p(14)*Y_prime;



 
