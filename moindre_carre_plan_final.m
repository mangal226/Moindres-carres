%% Initialisation

clc;clear;
pas=50; % Pas
x_init=-500:pas:500;% Tableau des coordonnées x
y_init=-500:pas:500; % Tableau des coordonnées y
z=-500:pas:500; % Tableau des coordonnées z
a=100; b=-500; c=1100; % Coefficients pour l'equation du plan

%Initialisation des sommes 
Sx=0;
Sy=0;
Sxx=0; 
Syy=0;
Sxy=0;
Sxz=0;
Syz=0;
Sz=0;

n_init=size(x_init,2); bruit=0.001; % Taille matrices coordonnées

n=1;% Initialisation pour calcul de la taille de la matrice de coordonnées globale

%% Allocation matrice coordonées 
for i=1:n_init
    for j=1:n_init
          n=n+1;
          M(1,n)=x_init(1,j);
          M(2,n)=y_init(1,i);
          M(3,n)=a*M(1,n)+b*M(2,n)+c;
    end
end

%% Tracé du nuage de points simples

% plot3(M(1,:),M(2,:),M(3,:),'x');
% xlabel('Valeur dentree x');
% ylabel('Valeur dentree y');
% zlabel('Valeur obtenue z');
% 
% hold on;

%% Insertion des bruits sur le nuage de points

% Insertion du bruit
rng(0,'twister');
ecart_type=0.25; % Parametrage de la valeur du bruit
moyenne=0;
valeur_bruit=ecart_type.*randn(n,1)+moyenne;

% Insertion des valeurs bruités
    M3_modif(3,:)=M(3,:);
% Bruitage pour des valeurs aleatoires
k=6;
while(k<n)
    M3_modif(3,k)=M3_modif(3,k)+valeur_bruit(k,1);
    M3_modif(3,k+3)=M3_modif(3,k-5)-valeur_bruit(k-5,1);
    k=k+10;
end
    

 %% Tracé du nuage de points bruités
 
%     plot3(M(1,:),M(2,:),M3_modif(3,:),'x');
%     xlabel('Valeur dentree x');
%     ylabel('Valeur dentree y');
%     zlabel('Valeur bruitée z');
%     hold on;
    
  
 %% Algorithme des moindres carrés
 
 % Calcul des sommes
    for k=1:n
        Sx=Sx+M(1,k);
        Sy=Sy+M(2,k);
        Sxx=Sxx+M(1,k)^2;
        Sxy=Sxy+M(1,k)*M(2,k);
        Syy=Syy+M(2,k)^2;
        Sxz=Sxz+M(1,k)*M3_modif(3,k);
        Syz=Syz+M(2,k)*M3_modif(3,k);
        Sz=Sz+M3_modif(3,k);
    end
    
 % Calcul des coefficients pour l'equation du plan
 B=[Sxz;Syz;Sz];
 A=[Sxx,Sxy,Sx;Sxy,Syy,Sy;Sx,Sy,n];
 Coeff=(inv(A))*B;% Inversion de matrice et obtention des coefficients
 
 
 %% Tracé du plan des moindres carrés
 
 [xPlan,yPlan]=meshgrid(x_init(1,:),y_init(1,:)); 
 zPlan=Coeff(1,1)*xPlan + Coeff(2,1)*yPlan+ Coeff(3,1);
 mesh(xPlan,yPlan,zPlan);
 
 xlabel('Valeur dentree x');
 ylabel('Valeur dentree y');
 zlabel('Plan des moindres carre z');
 hold on
 plot3(M(1,:),M(2,:),M3_modif(3,:),'x');
 hold off;
 %% Tracé des trois plans sur la même figure
%  
% tiledlayout(3,1);
% 
% nexttile;
% plot3(M(1,:),M(2,:),M(3,:));
% title('plan normal sans bruit');
% hold on
% 
% nexttile
% plot3(M(1,:),M(2,:),M3_modif(3,:),'x');
% title('nuage de points bruités');
% 
% nexttile;
% [xPlan,yPlan]=meshgrid(M(1,:),M(2,:));
% zPlan=Coeff(1,1)*xPlan + Coeff(2,1)*yPlan+Coeff(3,1);
% mesh(xPlan,yPlan,zPlan);
%  hold off;


%% Etude sur les bruits de mesure

% Calcul des erreurs de mésure
for t=1:n
    Erreur_mesure(1,t)= M3_modif(3,t)-Coeff(1,1)*M(1,t) -Coeff(2,1)*M(2,t)-Coeff(3,1); % Calcul erreur pour chaque mesure
end


% Calcul de l'erreur globale
    erreur=0;
    for q=1:n
        erreur=erreur+Erreur_mesure(1,q); 
    end

% Repartition des résidus autour de erreur=0
abscisse_bruit=1:1:n;
y_bruit=Erreur_mesure;

%     plot(abscisse_bruit,y_bruit,'o');
%     hold on;
%     yline(0);
% Determination de la densité
% for u=1:n
%     densite(1,u)=M3_modif(3,u)/Erreur_mesure(1,u);
% end
% histo_densite=histogram(densite);
%     plot(abscisse_bruit(1,:),densite(1,:),'x');
%     hold on;

% Nous avons fait plutot une etude sur le bruit de mesure lui même


% Trace du bruit de mesure

% plot(abscisse_bruit,y_bruit);
% hold on;


%Presentation de la repartition des erreurs de mesure sur un histogramme
% nombre_classe=50;
% histo_erreur=histogram(Erreur_mesure(1,:));

% Calcul du coefficient de determination 
% On essaye de voir si le modèle est fiable en calculant son coefficient de
% determination. Plus il tend vers 1 plus le modèle est fiable

% Calcul de la moyenne des Z mesurés
z_barre=Sz/n;

% Calcul de la variation expliquée par la régression 
variationeReg=0;
for p1=1:n
    variationeReg=variationeReg+ (Coeff(1,1)*M(1,p1)+ Coeff(2,1)*M(2,p1)+Coeff(3,1)-z_barre)^2;
    variance(1,p1)=variationeReg/(n-1);
end
% h2=histogram(variance);
% h3=histfit(variance);


% Calcul de la variation expliquée par les résidus
variationeRes=0;
for p2=1:n
    variationeRes=variationeRes+Erreur_mesure(1,p2)^2;
end
% 
 % Calcul de la variation totale
variationTotale=variationeReg+variationeRes;
% 
% Calcul du coefficient de determination/ Coefficient de correlation

R=variationeReg/variationTotale; % Coefficient de determination





    
            
    
        