%% Initialisation
clc;
clear;
x=-100:10:100; % Determination des coordonnées x

% Initialisation des sommes
Sx=0;
Sy=0;
Sxx=0;
Sxy=0;
a=20;
b=150;

n=size(x,2); % Determination de la taille de la matrice des coordonees x

%% Droite parfaite

% Creation des coordonnées y pour la droite
for i= 1:size(x,2)
    y(1,i)=a*x(1,i)+b;
end

% Tracé de la droite et du nuage de points la caracterisant
% plot(x,y,'x');
% hold on
% plot(x,y);
% xlabel('Valeurs dentrees x');
% ylabel('Valeur de sortie y');
% hold on;

%% Points bruités

% Alocation des valeurs avant bruitage
x_modif(1,:)=x(1,:);
y_modif(1,:)=y(1,:);


% Creation du nuage de points bruités
rng(0,'twister');

ecart_type=7; % Determination de la valeur du bruit
moyenne=0;

val_bruit=ecart_type.*randn(n,1)+moyenne; % Tableau contenant des valeurs de bruit
val_bruit=transpose(val_bruit); % Transposition du tableau obtenu pour une meilleure utilisation


% Bruitage des valeurs y de la droite
k=1;
while(k<n)
    y_modif(1,k)=y_modif(1,k)+100*val_bruit(1,k);
    y_modif(1,k+3)=y_modif(1,k+3)-100*val_bruit(1,k+3);
    k=k+7;
end
    
% Tracé du nuage de points bruité 
% plot (x_modif(1,:),y_modif(1,:),'x');
% xlabel('Valeurs dentrees x');
% ylabel('Valeur de sortie bruitees y');
% hold on;

%% Algorithme des moindres carrés

% Calcul des sommes
    for l=1:n
        Sx=Sx+x_modif(1,l);
        Sy=Sy+y_modif(1,l);
        Sxx=Sxx+(x_modif(1,l))^2;
        Sxy=Sxy+x(1,l)*y(1,l);
    end
    
    
% Calcul des coefficients de la droite 
A=[Sxy;Sy];
B=[Sxx,Sx;Sx,n];
Coeff=inv(B)*A; % Inversion de matrice et obtention des coefficients de la droite
            
% Droite des moindres carrés
for p=1:n
    y_moindre_carre(1,p)=Coeff(1,1)*x(1,p)+Coeff(2,1);
end

%Tracé de la droite des moindres carrés
plot(x_modif(1,:),y_moindre_carre(1,:),'r');
hold on;
plot(x_modif(1,:),y_modif(1,:),'x');
xlabel('Valeurs dentrees x');
ylabel('Droite des moindres carres y');
hold on;


% Tracé de toutes les droites simultanément
% tiledlayout(3,1)
% nexttile
% plot(x_modif(1,:),y(1,:));
% title('droite normale sans bruit');
% hold on
% nexttile
% plot(x_modif(1,:),y_modif(1,:),'x');
% title('nuage de points bruités');
% 
% nexttile
% plot(x_modif(1,:),y_moindre_carre(1,:));
% title('droite des moindres carres');



         
%% Calcul de l'erreur
    erreur=0;
    for q=1:n
        erreur=erreur+y_modif(1,q)-y_moindre_carre(1,q);
    end
        