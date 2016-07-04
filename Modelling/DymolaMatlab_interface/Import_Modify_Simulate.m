% Fonction permettant de simuler un modèle 2 en important auparavant les
% resultats d'un modele 1, de faire des modifs nécessaire sur le modele , et enfin de modifier un parametre pour l'étude de son
% influence

% ARGUMENTS:
% - mod1 : modèle qui fournit les resultats importés
% - mod2 : modèle à simuler
% - tps : temps d'importation des résultats du modèle 1
% - pack: package complet du modèle 2 à simuler
% - modif : vecteur des grandeurs à modifier dans le modele 2 apres
% importation des resultats du modele 1
% - param: parametre à etudier avec le modele 2
% - startT : temps de démarrage de la simul du modèle 2
% - stopT : temps d'arrêt de la simul du modèle 2
% - NOI : number of interval pour la simul du modèle 2
% - File : nom du fichier resultat du modèle 2
% - meth: solveur de la simul du modèle 2
% - tol : tolérance du solveur de la simul du modèle 2


function Import_Modify_Simulate(mod1,mod2,tps,pack,modif,param,startT,stopT,NOI,File,meth,tol)

% Compilation du modèle 2
commande=sprintf('translateModel("%s.%s");',pack,mod2);
dymolaM(commande);
clear('commande');

% Importation des résultats du modele 1
commande=sprintf('importInitialResult("%s.mat",%d);',mod1,tps);
dymolaM(commande);
clear('commande');

% Modification nécessaire et changement du paramètre à étudier
n=length(modif);
for i=1:n
   dymolaM(char(modif(i)));
end
dymolaM(char(param));

% Simulation du modele 2
commande=sprintf('simulateModel("%s.%s",startTime=%d,stopTime=%d,numberOfIntervals=%d,resultFile="%s",method="%s",tolerance=%f);',pack,mod2,startT,stopT,NOI,File,meth,tol);
dymolaM(commande);


end
