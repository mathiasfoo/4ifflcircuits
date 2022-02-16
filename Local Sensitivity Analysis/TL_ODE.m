function dx= TL_ODE(t,x,p)
global par 
dx = zeros(7,1);


%%%% Molecule Species

X = x(1);  
Y = x(2);  
Z = x(3);     
XY = x (4);
Y_prime = x(5); 
XZ = x (6);
GFP = x(7);

Px = par.P_x; 
Py = par.P_y; 
Pz = par.P_z; 

dx(1) = p(1)*Px - p(12)*X*Y - p(13)*X*Z - p(4)*X;

dx(2) = p(2)*Py - p(12)*X*Y - p(5)*Y; 

dx(3) = p(3)*Pz - p(11)*Z*Y_prime - p(13)*X*Z - p(6)*Z;

dx(4) = p(12)*X*Y - p(14)*XY;

dx(5) = p(9)*XY - p(11)*XZ*Y_prime - p(11)*Z*Y_prime - p(10)*Y_prime; 

dx(6) = p(13)*X*Z - p(11)*XZ*Y_prime - p(15)*XZ;

dx(7) = p(7)*XZ - p(8)*GFP; 
