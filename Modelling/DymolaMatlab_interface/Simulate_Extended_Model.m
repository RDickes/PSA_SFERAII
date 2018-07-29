% Fonction permettant d'appeler la fonction SimulateExtendedModel de Dymola
%- Génere des fichiers de la forme "nom_fic"_"valeur_param".mat
%- Nécessite que le modèle en question soit ouvert dans Dymola

% ARGUMENTS:
% mod: nom du modele Dymola sans le .mo (char)
% startT: temps initial (entier)
% stopT: temps final (entier)
% NOI: number of intervals (entier)
% meth: solveur dymola (char)
% param: parametre a etudier (char)
% valeur_param : valeur du parametre à étudier (reel)
% nom_fic: nom générique du fichier généré (char)
% var : variable de sortie (char)
% pack : nom(total) du package auquel appartient le modele (char)
% chem_res: repertoire d'ecriture du fichier .mat(char) dossier/dossier/

function Simulate_Extended_Model(mod,startT,stopT,NOI,meth,param,valeur_param,nom_fic,var,pack,chem_res)
disp(sprintf('EXECUTION DU MODELE DYMOLA %s',mod));

disp('-Preparation des arguments pour la fonction dymola');
% Préparation des arguments pour la fonction Dymola
modele=sprintf('%s.%s',pack,mod);
File=sprintf('%s%s_%f',chem_res,nom_fic,valeur_param);
initial_names=sprintf('{"%s"}',param);
initial_values=sprintf('{%f}',valeur_param);

disp('-Appel de la fonction Dymola');
% Appel de la fonction Dymola
commande=sprintf('simulateExtendedModel(problem="%s",startTime=%d,stopTime=%d,numberOfIntervals=%d,method="%s",initialNames=%s,initialValues=%s,finalNames={"%s"},resultFile="%s");',modele,startT,stopT,NOI,meth,initial_names,initial_values,var,File);
dymolaM(commande);

end
