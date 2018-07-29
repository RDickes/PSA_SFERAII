% Fonction permettant d'appeler la fonction SimulateModel de Dymola
%- Utilise le mod�le du r�pertoire courant et g�n�re les fichiers r�sultats
% dans le m�me r�pertoire
%- G�nere un fichier "nom_fic".mat
%- N�cessite d'avoir le modele en question ouvert dans Dymola

% ARGUMENTS:

% mod: nom du modele sans le .mo (char)
% startT: temps initial (reel)
% stopT: temps final (reel)
% meth: solveur dymola (char)
% pack : nom(total) du package auquel appartient le modele (char)
% chem_res : repertoire du fichier .mat (char): dossier/dossier/
% nom_fic : nom du fichier resultats g�n�r�

function Simulate_Model(mod,startT,stopT,NOI,meth,pack,chem_res,nom_fic,tol)
disp(sprintf('EXECUTION DU MODELE DYMOLA %s',mod));

disp('nettoyage environnement Dymola');
%Nettoyage de l'environnement Dymola
% dymolaM('clearlog();');
% dymolaM('removeResults();');
% dymolaM('removePlots();');

disp('preparation des arguments pour la fonction dymola');
% Pr�paration des arguments pour la fonction Dymola
modele=sprintf('%s.%s',pack,mod);
File=sprintf('%s%s',chem_res,mod);

disp('appel de la fonction Dymola');
% Appel de la fonction Dymola
commande=sprintf('simulateModel("%s",startTime=%d,stopTime=%d,numberOfIntervals=%d,resultFile="%s",method="%s",tolerance=%f);',modele,startT,stopT,NOI,File,meth,tol);
dymolaM(commande);

end
