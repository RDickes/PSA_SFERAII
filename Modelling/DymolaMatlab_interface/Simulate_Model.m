% Fonction permettant d'appeler la fonction SimulateModel de Dymola
%- Utilise le modèle du répertoire courant et génère les fichiers résultats
% dans le même répertoire
%- Génere un fichier "nom_fic".mat
%- Nécessite d'avoir le modele en question ouvert dans Dymola

% ARGUMENTS:

% mod: nom du modele sans le .mo (char)
% startT: temps initial (reel)
% stopT: temps final (reel)
% meth: solveur dymola (char)
% pack : nom(total) du package auquel appartient le modele (char)
% chem_res : repertoire du fichier .mat (char): dossier/dossier/
% nom_fic : nom du fichier resultats généré

function Simulate_Model(mod,startT,stopT,NOI,meth,pack,chem_res,nom_fic,tol)
disp(sprintf('EXECUTION DU MODELE DYMOLA %s',mod));

disp('nettoyage environnement Dymola');
%Nettoyage de l'environnement Dymola
% dymolaM('clearlog();');
% dymolaM('removeResults();');
% dymolaM('removePlots();');

disp('preparation des arguments pour la fonction dymola');
% Préparation des arguments pour la fonction Dymola
modele=sprintf('%s.%s',pack,mod);
File=sprintf('%s%s',chem_res,mod);

disp('appel de la fonction Dymola');
% Appel de la fonction Dymola
commande=sprintf('simulateModel("%s",startTime=%d,stopTime=%d,numberOfIntervals=%d,resultFile="%s",method="%s",tolerance=%f);',modele,startT,stopT,NOI,File,meth,tol);
dymolaM(commande);

end
