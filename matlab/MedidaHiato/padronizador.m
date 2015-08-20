function padronizado = padronizador(dados)
% Padroniza "dados"

for iii = 1:size(dados,2)
    padronizado(:,iii) = (dados(:,iii) - mean(dados(:,iii)))/std(dados(:,iii));
end