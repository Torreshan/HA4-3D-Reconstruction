function [T1,R1,T2,R2] = TR_aus_E(E)
% In dieser Funktion sollen die moeglichen euklidischen Transformationen
% aus der Essentiellen Matrix extrahiert werden
[U,s,V] = svd(E);
%s = [1,0,0;0,1,0;0,0,0];
Rz_p = [0,-1,0;1,0,0;0,0,1];
Rz_n = [0,1,0;-1,0,0;0,0,1];
%*************************************************************************%
R1 = U*Rz_p'*V';
R2 = U*Rz_n'*V';
%*************************************************************************%
T1_hat = U*Rz_p*s*U'; 
T1 = [T1_hat(3,2);T1_hat(1,3);T1_hat(2,1)];
T2_hat = U*Rz_n*s*U';
T2 = [T2_hat(3,2);T2_hat(1,3);T2_hat(2,1)];
end