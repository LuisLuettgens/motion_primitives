clear
clc

mkdir('output');
pointFile = fopen('output/legendreLobattoPoints.h','w');
weightFile = fopen('output/legendreLobattoWeights.h','w');

fprintf(pointFile,'/* MATLAB-Skript: LegendreLobatto.m */\n\n');
fprintf(pointFile,'#include <vector>\n\n');
fprintf(pointFile,'static std::vector<double> getLegendreLobattoPoints(const int i) {\n');
fprintf(pointFile,'\tswitch(i) {\n');

fprintf(weightFile,'/* MATLAB-Skript: LegendreLobatto.m */\n\n');
fprintf(weightFile,'#include <vector>\n\n');
fprintf(weightFile,'static std::vector<double> getLegendreLobattoWeights(const int i) {\n');
fprintf(weightFile,'\tswitch(i) {\n');

% erster Punkt {-1,1} von Hand eintragen
fprintf(pointFile,'\t\tcase %i:\n', 2);
fprintf(pointFile,'\t\t\treturn {\n');
fprintf(pointFile,'\t\t\t\t%.16f,\n',-1);
fprintf(pointFile,'\t\t\t\t%.16f\n',+1);
fprintf(pointFile,'\t\t\t};\n');

% Gewichte fuer ersten Punkt
fprintf(weightFile,'\t\tcase %i:\n', 2);
fprintf(weightFile,'\t\t\treturn {\n');
fprintf(weightFile,'\t\t\t\t%.16f,\n',1);
fprintf(weightFile,'\t\t\t\t%.16f\n',1);
fprintf(weightFile,'\t\t\t};\n');

% Genauigkeit erhoehen (noetig fuer Bestimmung der Gewichte)
digits(128);

N = 200;
for i = 3:N
    fprintf('step %i of %i\n', i-2,N-2)
    
    syms x
    
    % loesen
    roots = vpasolve(diff(legendreP(i-1,x),x) == 0);
    
    % Gewichte bestimmen
    x = roots;
    weights = vpa(2)./(i*(i-1)*(legendreP(i-1,x).^2));
    
    % Gewichte schreiben ANFANG
    fprintf(weightFile,'\t\tcase %i:\n', i);
    fprintf(weightFile,'\t\t\treturn {\n');
    fprintf(weightFile,'\t\t\t\t%.16f,\n',2/(i*(i-1)));
    fprintf(weightFile,'\t\t\t\t%.16f,\n',weights);
    fprintf(weightFile,'\t\t\t\t%.16f\n',2/(i*(i-1)));
    fprintf(weightFile,'\t\t\t};\n');
    % Gewichte schreiben ENDE
    
    % Nullstellen schreiben ANFANG
    fprintf(pointFile,'\t\tcase %i:\n', i);
    fprintf(pointFile,'\t\t\treturn {\n');
    fprintf(pointFile,'\t\t\t\t%.16f,\n',-1);
    fprintf(pointFile,'\t\t\t\t%.16f,\n',roots);
    fprintf(pointFile,'\t\t\t\t%.16f\n',+1);
    fprintf(pointFile,'\t\t\t};\n');
    % Nullstellen schreiben ENDE
end

fprintf(pointFile,'\t\tdefault:\n');
fprintf(pointFile,'\t\t\tstd::cout << "There are no roots for " << i << "!" << std::endl;\n');
fprintf(pointFile,'\t\t\texit(1);\n');
fprintf(pointFile,'\t\t\treturn {\n');
fprintf(pointFile,'\t\t\t\t%.16f\n',0);
fprintf(pointFile,'\t\t\t};\n');
fprintf(pointFile,'\t}\n');
fprintf(pointFile,'}\n');

fprintf(weightFile,'\t\tdefault:\n');
fprintf(weightFile,'\t\t\tstd::cout << "There are no weights for " << i << "!" << std::endl;\n');
fprintf(weightFile,'\t\t\texit(1);\n');
fprintf(weightFile,'\t\t\treturn {\n');
fprintf(weightFile,'\t\t\t\t%.16f\n',0);
fprintf(weightFile,'\t\t\t};\n');
fprintf(weightFile,'\t}\n');
fprintf(weightFile,'}\n');

fclose(pointFile);
fclose(weightFile);

% vpasolve(diff(legendreP(i,x),x) == 0)