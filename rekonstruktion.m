function [T,R, lambdas, P1] = rekonstruktion(T1,T2,R1,R2, Korrespondenzen, K)
% Funktion zur Bestimmung der korrekten euklidischen Transformation, der
% Tiefeninformation und der 3D Punkte der Merkmalspunkte in Bild 1
[~,n] = size(Korrespondenzen);
x1 = [Korrespondenzen(1:2,:);ones(1,n)];
x1 = K\x1;
x2 = [Korrespondenzen(3:4,:);ones(1,n)];
x2 = K\x2;

%*************************************************************************%
%*********************    Initialisierung  *******************************%
M11 = zeros(3*n,n+1);
M12 = zeros(3*n,n+1);
M21 = zeros(3*n,n+1);
M22 = zeros(3*n,n+1);

m11 = zeros(3*n,n+1);
m12 = zeros(3*n,n+1);
m21 = zeros(3*n,n+1);
m22 = zeros(3*n,n+1);

R = zeros(3,3);
T = zeros(3,1);
lambda1 = {};
lambda2 = {};
V = zeros(2,2);

%*************************************************************************%
%********************** lambda 1 rechnen      ****************************%
for i = 1:n
M11((3*i-2):3*i,i) = skew(x2(:,i))*R1*x1(:,i);
M11((3*i-2):3*i,n+1) = skew(x2(:,i))*T1;
M12((3*i-2):3*i,i) = skew(x2(:,i))*R1*x1(:,i);
M12((3*i-2):3*i,n+1) = skew(x2(:,i))*T2;
M21((3*i-2):3*i,i) = skew(x2(:,i))*R2*x1(:,i);
M21((3*i-2):3*i,n+1)= skew(x2(:,i))*T1;
M22((3*i-2):3*i,i) = skew(x2(:,i))*R2*x1(:,i);
M22((3*i-2):3*i,n+1) = skew(x2(:,i))*T2;
end
lambda1{1,1} = MEIG(M11);
V(1,1) = length(find(lambda1{1,1}>0));
lambda1{1,2} = MEIG(M12);
V(1,2) = length(find(lambda1{1,2}>0));
lambda1{2,1} = MEIG(M21);
V(2,1) = length(find(lambda1{2,1}>0));
lambda1{2,2} = MEIG(M22);
V(2,2) = length(find(lambda1{2,2}>0));

%*************************************************************************%
%********************** lambda 2  rechnen ********************************%
for i = 1:n
m11((3*i-2):3*i,i) = skew(x1(:,i))*R1'*x2(:,i);
m11((3*i-2):3*i,n+1) = skew(x1(:,i))*R1'*T1;
m12((3*i-2):3*i,i) = skew(x1(:,i))*R1'*x2(:,i);
m12((3*i-2):3*i,n+1) = skew(x1(:,i))*R1'*T2;
m21((3*i-2):3*i,i) = skew(x1(:,i))*R2'*x2(:,i);
m21((3*i-2):3*i,n+1)= skew(x1(:,i))*R2'*T1;
m22((3*i-2):3*i,i) = skew(x1(:,i))*R2'*x2(:,i);
m22((3*i-2):3*i,n+1) = skew(x1(:,i))*R2'*T2;
end
lambda2{1,1} = MEIG(m11);
V(1,1) = V(1,1)+length(find(lambda2{1,1}>0));
lambda2{1,2} = MEIG(m12);
V(1,2) = V(1,2)+length(find(lambda2{1,2}>0));
lambda2{2,1} = MEIG(m21);
V(2,1) = V(2,1)+length(find(lambda2{2,1}>0));
lambda2{2,2} = MEIG(m22);
V(2,2) = V(2,2)+length(find(lambda2{2,2}>0));

%*************************************************************************%
%******************* die R,T,lambdas und P1 bestimmen   ******************%
[x,y]=find(V==max(max(V)),1);
switch x
    case 1
        R = R1;
    case 2
        R = R2;
end
switch y
    case 1
        T = T1;
    case 2
        T = T2;
end    
lambdas = lambda1{x,y};
gamma = lambdas(1,n+1); % gamma muss positiv sein 
T = gamma*T; % skalierung
P1 = x1.*repmat(lambdas(1,1:n),3,1);
end

function[mV] = MEIG(M)
    %die Rechnung der minimum sigular value und des entsprechenden Vector
   [V,D] = eig(M'*M);
   [~,index]=find(D ==min(min(diag(D))));
   nV = V(:,index);
   mV = nV';
   if( mV(1,end) < 0) %gamma
       mV = -mV;
   end
end
function [Vhat] = skew(V)
    % Umwandlung von V in eine schiefsymmetrische Matrix
    Vhat = [0 -V(3) V(2); V(3) 0 -V(1); -V(2) V(1) 0];
end