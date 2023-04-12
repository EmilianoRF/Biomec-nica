 function [AMostrar] = InterpolaA100Muestras(AMostrar); 
%  recibe una lista de ordenadas en la variable AMostrar y la interpola a 100 muestras de ordenadas;
% Genera la absisa original con la variable Tiempo
% Genera la nueva absisa  con la variable Porcentaje
% Devuelve sobreescriviendo la variable AMostrar con 10 muestras
 Tiempo = 1:1:length(AMostrar);
 
 Porcentaje = 1:length(AMostrar)/100:length(AMostrar);
 
AMostrar  = interp1(Tiempo,AMostrar,Porcentaje,'spline');