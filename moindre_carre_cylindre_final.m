         %% Initialisation
clc;clear;
pas=50;%  Pas de deplacement vertical i.e. distance verticale entre 2 points
z_init=-500:pas:500; % Hauteur du cylindre
n_init=size(z_init,2); 
pas_angle=(2*pi/n_init); % Angle de deplacement sur le plan xy pour creer cercle
r=100; % Rayon du cylindre
teta=-pi:pas_angle:pi; 


%% Nuage de points reguliers pour un cercle
  for i=1:n_init
        x_init(1,i)=r*cos(teta(1,i));% Coordon�es x des points du cylindre
        y_init(1,i)=r*sin(teta(1,i));% Coordonn�es y des points du cylindre
  end
  
  
  
  %%  Creation du nuage de points constituant le cylindre
  j=0; 
  for k=1:n_init
      for l=1:n_init
          j=j+1;
          M(j,1)=x_init(1,l);
          M(j,2)=y_init(1,l);
          M(j,3)=z_init(1,k);
      end
  end
n=j;  % Receuil du nombre de points du cylindre; 

 
 %% Trac� du cylindre normal
%       plot3(M(:,1),M(:,2),M(:,3),'x'); 
%       xlabel('Valeur dentree X');
%       ylabel('Valeur dentree Y');
%       zlabel('Cylindre normal');
%       hold on;


%% Cr�ation du nuage de points bruit�s

    % Insertion du bruit
    rng(0,'twister');
    ecart_type=50;
    moyenne=0;
    bruit=ecart_type.*randn(n,1)+moyenne;
    M_modif=M;
    
    p=1;
    % Bruitage de certaines valeurs choisies aleatoirement
    
    while(p<n)
        M_modif(p,1)=M(p,1)+bruit(p,1); % Ajout de bruit sur la coordon�e x
        M_modif(p,2)=M(p,2)+bruit(p,1); % Ajout de bruit sur la coordonn�e y
        
        M_modif(p+1,1)=M(p+1,1)-bruit(p+1,1); % Ajout de bruit sur la coordon�e x
        M_modif(p+1,2)=M(p+1,2)-bruit(p+1,1); % Ajout de bruit sur la coordonn�e y
        p=p+50;
    end
          

%% Trac� du cylindre bruit�

%       plot3(M_modif(:,1),M_modif(:,2),M(:,3),'x');
%       xlabel('Valeur dentree X');
%       ylabel('Valeur dentree Y');
%       zlabel('Cylindre bruit�');
% hold on 
  
 %% Algorithme de Levenberg Marquard

% Initialisation des valeurs initiales 
Solution(1,1)=100;
Solution(2,1)=100;
Solution(3,1)=100;
Solution(4,1)=100;

% Affectation des valeurs initiales
a=Solution(1,1);
b=Solution(2,1);
xo=Solution(3,1);
yo=Solution(4,1);
ro(1,1)=250;
residu(1,1)=100;


s=2; % Initialisation de la condition sur la boucle while
temp=100; % Initialisation de la condition sur la boucle while

% Debut de la boucle while
while( s<n && abs(temp)>0.05)% Critere d'arret si la valeur absolue de l'erreur inf a 0.05
    
    % Intialisation des valeurs fixes de la matrice Jacobienne
    xi=M_modif(s,1);
    yi=M_modif(s,2);
    zi=M_modif(s,3);
    zo=0;
    
    % D�termination de la valeur de ro
    % ro est la moyenne des distances au point de coordonn�es xo, yo, zo apr�s chaque
    % it�ration
    somme_distance=0;
    for p=1:n
        xi=M_modif(p,1);
        yi=M_modif(p,2);
        zi=M_modif(p,3);
        somme_distance=somme_distance+((a*(yi-yo)-b*(xi-xo))^2+(a*(zi-zo)-(xi-xo))^2+(b*(zi-zo)-(yi-yo))^2)^0.5/((a^2+b^2+1)^0.5);
    end
    
    ro(1,s)=somme_distance/n; % Calcul de ro la valeur du rayon pour l'iteration s
    
    
    r=ro(1,s);
    % Ecriture de la fonction � deriver
    syms a b xo yo 
    val=(((a*(yi-yo)-b*(xi-xo))^2+(a*(zi-zo)-(xi-xo))^2+(b*(zi-zo)-(yi-yo))^2)^0.5/((a^2+b^2+1)^0.5)-r)^2; % Expression de la fonction � d�river

    % Calcul de la jacobienne
    J(1,1)=diff(val,a);
    J(1,2)=diff(val,b);
    J(1,3)=diff(val,xo);
    J(1,4)=diff(val,yo);
    

    % Calcul de la transpose de la jacobienne
    K=transpose(J);

    % Allocation des valeurs num�riques aux coefficients de la jacobienne
    a=Solution(1,s-1);
    b=Solution(2,s-1);
    xo=Solution(3,s-1);
    yo=Solution(4,s-1);
    
    % Calcul des valeurs num�riques de la Jacobienne obtenue et sa
    % transpos�e
    Jacob=double(subs(J));
    TransJacob=double(subs(K));

    % Valeur du residu
    residu(s,1)=(((a*(yi-yo)-b*(xi-xo))^2+(a*(zi-zo)-(xi-xo))^2+(b*(zi-zo)-(yi-yo))^2)^0.5/((a^2+b^2+1)^0.5)-r); % Calcul du residu pour l'iteration s
    temp=residu(s,1); % Allocation de valeur pour le crit�re d'arr�t
    
    
    % Allocation du coefficient d'ammortissement
    lambda_k=0.001;% Choix d'un coefficient d'ammortissement il est choisit entre 0 et 1

    % Calcul de la solution pour l'ordre suivant i.e. s+1
    Solution(:,s)=inv(double(TransJacob*Jacob+lambda_k*eye(4)))*TransJacob*residu(s,1);
    Solution_optimale=Solution(:,s);
    s=s+1;

end

% Presentation de la solution finale
Solution_finale =transpose(Solution(:,s-1));
r_optimal=ro(1,s-1);

 %% Trac� du cylindre des moindres carres

 %Creation du nuage de points du cylindre a partir du rayon optimal 
  j=0;
  for q1=1:n_init
      for q2=1:n_init
          j=j+1;
          M_cylindre_final(j,1)=r_optimal*cos(teta(1,q2))+xo;
          M_cylindre_final(j,2)=r_optimal*sin(teta(1,q2))+yo;
          M_cylindre_final(j,3)=z_init(1,q1);
      end
  end

  
% Trac� du cylindre des moindres carr�s
% plot3(M_cylindre_final(:,1),M_cylindre_final(:,2),M_cylindre_final(:,3));
% hold on;

% Trace de l'ensemble: cylindre moindres carres et points bruit�s
plot3(M_cylindre_final(:,1),M_cylindre_final(:,2),M_cylindre_final(:,3));
hold on;
      plot3(M_modif(:,1),M_modif(:,2),M(:,3),'x');

      xlabel('Valeur dentree X');
      ylabel('Valeur dentree Y');
      zlabel('Cylindre des moindres carres');
hold on;


% Trac� des trois cylindres sur la m�me figure
 
% tiledlayout(3,1)
% 
% nexttile
% plot3(M(:,1),M(:,2),M(:,3));
% title('cylindre normal sans bruit');
% hold on
% 
% nexttile
% plot3(M_modif(:,1),M_modif(:,2),M_modif(:,3),'x');
% title('nuage de points bruit�s');
% 
% nexttile
% plot3(M_cylindre_final(:,1),M_cylindre_final(:,2),M_cylindre_final(:,3),'x');
% title('cylindre des moindres carres');
% hold off;


