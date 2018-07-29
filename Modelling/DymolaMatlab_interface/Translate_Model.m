% Fonction permettant d'appeler la fonction TranslateModel de Dymola
%- Compile le mod�le en argument dans le r�pertoire courant
%- N�cessite que le mod�le en question soit ouvert dans Dymola

% ARGUMENTS:
% mod: nom du modele Dymola sans le .mo (char)
% pack : nom(total) du package auquel appartient le modele (char)


function Translate_Model(mod,pack)
disp(sprintf('COMPILATION DU MODELE DYMOLA %s',mod));

%Nettoyage de l'environnement Dymola
dymolaM('clearlog();');
dymolaM('removeResults();');
dymolaM('removePlots();');

% Pr�paration des arguments pour la fonction Dymola
modele=sprintf('%s.%s',pack,mod);

% Appel de la fonction Dymola
commande=sprintf('translateModel("%s");',modele);
dymolaM(commande);

end