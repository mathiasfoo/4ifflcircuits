function dx= HY_TY1_ODE(t,x,p)
global par 
dx = zeros(4,1);


%%%%molecule species

X = x(1);  
Y = x(2);  
Z = x (3);
GFP = x(4);


Px = par.P_x; 
Py = par.P_y; 
Pz = par.P_z; 

Py_active = p(9)*X*Py; 
Pz_active = p(10)*X*Pz;

dx(1) = p(1)*Px - p(10)*X*Pz- p(9)*X*Py - p(4)*X;

dx(2) = p(2)*Py_active - p(11)*Y*Z - p(5)*Y;

dx(3) = p(3)*Pz_active - p(11)*Y*Z - p(6)*Z;

dx(4) = p(7)*Z - p(8)*GFP;



 
