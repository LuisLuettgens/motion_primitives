clear
clc

mkdir('output');
pointFile = fopen('output/legendreGaussPoints.h','w');
weightFile = fopen('output/legendreGaussWeights.h','w');

fprintf(pointFile,'/* MATLAB-Skript: LegendreGauss.m */\n\n');
fprintf(pointFile,'#include <vector>\n\n');
fprintf(pointFile,'static std::vector<double> getLegendreGaussPoints(const int i) {\n');
fprintf(pointFile,'\tswitch(i) {\n');

fprintf(weightFile,'/* MATLAB-Skript: LegendreGauss.m */\n\n');
fprintf(weightFile,'#include <vector>\n\n');
fprintf(weightFile,'static std::vector<double> getLegendreGaussWeights(const int i) {\n');
fprintf(weightFile,'\tswitch(i) {\n');

% Genauigkeit erhoehen (noetig fuer Bestimmung der Gewichte)
digits(128);

N = 200;
for i = 1:N
    fprintf('step %i of %i\n', i,N);
    
    syms x;
    
    % loesen
    roots = vpasolve(legendreP(i,x) == 0);
    
    % Gewichte bestimmen
    deriv = diff(legendreP(i,x));
    x = roots;
    weights = vpa(2)./((vpa(1)-(roots.^2)).*subs(deriv).^2);
    
    % Gewichte schreiben ANFANG
    fprintf(weightFile,'\t\tcase %i:\n', i);
    fprintf(weightFile,'\t\t\treturn {\n');
    if (length(weights) > 1)
    fprintf(weightFile,'\t\t\t\t%.16f,\n',weights(1:end-1));
    end
    fprintf(weightFile,'\t\t\t\t%.16f\n',weights(end));
    fprintf(weightFile,'\t\t\t};\n');
    % Gewichte schreiben ENDE
    
    % Nullstellen schreiben ANFANG
    fprintf(pointFile,'\t\tcase %i:\n', i);
    fprintf(pointFile,'\t\t\treturn {\n');
    if (length(roots) > 1)
    fprintf(pointFile,'\t\t\t\t%.16f,\n',roots(1:end-1));
    end
    fprintf(pointFile,'\t\t\t\t%.16f\n',roots(end));
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
