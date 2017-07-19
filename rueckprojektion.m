function repro_error = rueckprojektion(Korrespondenzen, P1, I2, T, R, K)
% Diese Funktion berechnet die projizierten Punkte in Kamera 2 und den
% mittleren Rueckprojektionsfehler
II = [1,0,0,0;0,1,0,0;0,0,1,0];
[~,n] = size(P1);
P2 = R*P1 + repmat(T,1,n);
P2 = [P2;ones(1,n)];
%x2_prediction = K*(P2./(repmat(P2(3,:),3,1)))
x2_prediction = K*II*P2./(repmat(P2(3,:),3,1));
% Zeige die robusten Korrespondenzpunktpaare
figure('name', 'Rueckprojektion');
imshow(uint8(I2))
hold on
plot(x2_prediction(1,:),x2_prediction(2,:),'r*')
hold on
plot(Korrespondenzen(3,:),Korrespondenzen(4,:),'g*')
for i=1:size(Korrespondenzen,2)
    hold on
    x_1 = [x2_prediction(1,i), Korrespondenzen(3,i)];
    x_2 = [x2_prediction(2,i), Korrespondenzen(4,i)];
    line(x_1,x_2);
end
hold off

%**********************    Error  rechnen    *****************************%
er = x2_prediction-[Korrespondenzen(3:4,:);ones(1,n)];
%er = x2_prediction-[Korrespondenzen(4,:);Korrespondenzen(3,:);ones(1,n)];
S = 0 ;
for i =1:n
S = S + norm(er(:,i));
end
repro_error = S/n;
end