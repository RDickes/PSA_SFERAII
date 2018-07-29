% Fonction permettant de simuler un mod�le 2 en important auparavant les
% resultats d'un modele 1, de faire des modifs n�cessaire sur le modele , et enfin de modifier un parametre pour l'�tude de son
% influence

% ARGUMENTS:
% - mod1 : mod�le qui fournit les resultats import�s
% - mod2 : mod�le � simuler
% - tps : temps d'importation des r�sultats du mod�le 1
% - pack: package complet du mod�le 2 � simuler
% - modif : vecteur des grandeurs � modifier dans le modele 2 apres
% importation des resultats du modele 1
% - param: parametre � etudier avec le modele 2
% - startT : temps de d�marrage de la simul du mod�le 2
% - stopT : temps d'arr�t de la simul du mod�le 2
% - NOI : number of interval pour la simul du mod�le 2
% - File : nom du fichier resultat du mod�le 2
% - meth: solveur de la simul du mod�le 2
% - tol : tol�rance du solveur de la simul du mod�le 2


function Import_Modify_Simulate(mod1,mod2,tps,pack,modif,param,startT,stopT,NOI,File,meth,tol)

% Compilation du mod�le 2
commande=sprintf('translateModel("%s.%s");',pack,mod2);
dymolaM(commande);
clear('commande');

% Importation des r�sultats du modele 1
commande=sprintf('importInitialResult("%s.mat",%d);',mod1,tps);
dymolaM(commande);
clear('commande');

% Modification n�cessaire et changement du param�tre � �tudier
n=length(modif);
for i=1:n
   dymolaM(char(modif(i)));
end
dymolaM(char(param));

% Simulation du modele 2
commande=sprintf('simulateModel("%s.%s",startTime=%d,stopTime=%d,numberOfIntervals=%d,resultFile="%s",method="%s",tolerance=%f);',pack,mod2,startT,stopT,NOI,File,meth,tol);
dymolaM(commande);


end
