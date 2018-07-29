% Fonction permettant d'appeler la fonction importInitialResult de Dymola
%- Importe un fichier de resulat .mat dans l'environnement actuel Dymola


% ARGUMENTS:
% fic: chemin et nom du fichier resultat .mat (char)
% tps : temps des données à importer (integer)


function importInitialResult(fic,tps)
disp(sprintf('IMPORTATION DES RESULTATS DU FICHIER %s',fic));

disp('Appel de la fonction Dymola importation');
% Appel de la fonction Dymola
commande=sprintf('importInitialResult("%s",%d);',fic,tps);
dymolaM(commande);

end