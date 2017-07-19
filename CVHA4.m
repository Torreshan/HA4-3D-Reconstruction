%  Gruppennummer: M08
%  Gruppenmitglieder: Han,Fengze;Qu,Xingwei;Gao,Shangyin;Chen,Zhong;Wu,Yingyu.

%% Hausaufgabe 4
%  Bestimmung der euklidischen Bewegung und der 3D Rekonstruktion aus einem Stereobildpaar. 

%  Für die letztendliche Abgabe bitte die Kommentare in den folgenden Zeilen
%  enfernen und sicherstellen, dass alle optionalen Parameter über den
%  entsprechenden Funktionsaufruf fun('var',value) modifiziert werden können.


%% Bilder laden
Image1 = imread('szeneL.jpg');
IGray1 = rgb_to_gray(Image1);

Image2 = imread('szeneR.jpg');
IGray2 = rgb_to_gray(Image2);

%% Harris-Merkmale berechnen
Merkmale1 = harris_detektor(IGray1,'segment_length',9,'k',0.05,'min_dist',80,'N',50,'do_plot',false);
Merkmale2 = harris_detektor(IGray2,'segment_length',9,'k',0.05,'min_dist',80,'N',50,'do_plot',false);



%% Korrespondenzschätzung
tic;
Korrespondenzen = punkt_korrespondenzen(IGray1,IGray2,Merkmale1,Merkmale2,'window_length',25,'min_corr',0.9,'do_plot',false);
zeit_korrespondenzen = toc;
disp(['Es wurden ' num2str(size(Korrespondenzen,2)) ' Korrespondenzpunktpaare in ' num2str(zeit_korrespondenzen) 's gefunden.'])



%%  Finde robuste Korrespondenzpunktpaare mit Hilfe des RANSAC-Algorithmus
Korrespondenzen_robust = F_ransac(Korrespondenzen,'tolerance',0.015);
disp(['Es wurden ' num2str(size(Korrespondenzen_robust,2)) ' robuste Korrespondenzpunktpaare mittels RanSaC bestimmt.'])

% Zeige die robusten Korrespondenzpunktpaare
figure('name', 'Punkt-Korrespondenzen nach RANSAC');
imshow(uint8(IGray1))
hold on
plot(Korrespondenzen_robust(1,:),Korrespondenzen_robust(2,:),'r*')
imshow(uint8(IGray2))
alpha(0.5);
hold on
plot(Korrespondenzen_robust(3,:),Korrespondenzen_robust(4,:),'g*')
for i=1:size(Korrespondenzen_robust,2)
    hold on
    x_1 = [Korrespondenzen_robust(1,i), Korrespondenzen_robust(3,i)];
    x_2 = [Korrespondenzen_robust(2,i), Korrespondenzen_robust(4,i)];
    line(x_1,x_2);
end
hold off


%% Berechne die Essentielle Matrix
load('K.mat');
E = achtpunktalgorithmus(Korrespondenzen_robust,K);
disp(E);


%% Extraktion der möglichen euklidischen Bewegungen aus der Essentiellen Matrix und 3D-Rekonstruktion der Szene
 [T1,R1,T2,R2] = TR_aus_E(E);
 [T,R,lambdas,P1] = rekonstruktion(T1,T2,R1,R2,Korrespondenzen_robust,K);

%% Berechnung des mittleren Rückprojektionsfehlers auf der Bildebene von Kamera 2
repro_error = rueckprojektion(Korrespondenzen_robust, P1, IGray2, T, R, K);
disp('Error ist');
disp(repro_error)

