%% Practica 1 (08/03/23)

clear all % borra las variables
close all % cierra las ventanas
clc % limpia la consola


ruta_archivo = 'C:\Users\alumno\Desktop\Clase_1\0037_Davis_MarchaDavis_Walking02b2021.c3d';

%
h = btkReadAcquisition(ruta_archivo);

[marcadores, informacionCine] = btkGetMarkers(h);

[fuerzas, informacionFuerzas] = btkGetForcePlatforms(h);

antropometria = btkFindMetaData(h,'Antropometria');

Eventos=btkGetEvents(h);


% A1 es la variable del modelo que se usa para guardar el peso
A1 = antropometria.children.PESO.info.values';

%
A2 = antropometria.children.LONGITUD_ASIS.info.values';



% Se lee la frecuencia de muestreo de la camara

fm = informacionCine.frequency(1);

% Se calculan las muestras correspondientes a los eventos del pie
% derecho e izquierdo, con esto despues se ingresa al vector de marcadores y se
% sacan las posiciones.
RHS1 = round(Eventos.Derecho_RHS(1)*fm);
RHS2 = round(Eventos.Derecho_RHS(2)*fm);

LHS1 = round(Eventos.Izquierdo_LHS(1)*fm);
LHS2 = round(Eventos.Izquierdo_LHS(1*fm);


inicio = RHS1;
fin = LHS2;

% marcadores.r_met(fila_inicial:fila_final,columna inicial:columna final)
% si no se pone nada se toma todo
p1=marcadores.r_met(inicio:fin,:);









