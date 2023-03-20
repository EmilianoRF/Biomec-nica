clear all % borra las variables
close all % cierra las ventanas
clc % limpia la consola

% ===========================================================================================================================================================%
%  VARIABLES EXTRAS

rojo  = [1 0 0];
verde  = [0.051 0.529 0.051 ];
azul = [0 0 1];
negro = [0 0 0];
factor_vectores = 1/20;
% ===========================================================================================================================================================%
%  SE IMPORTAN LOS DATOS DE LA MARCHA

ruta_archivo = 'C:\Users\Emiliano\Desktop\Biomec�nica\Pr�ctica\Biomecanica\0037_Davis_MarchaDavis_Walking02b2021.c3d';
h = btkReadAcquisition(ruta_archivo);
[marcadores, informacionCine] = btkGetMarkers(h);
[fuerzas, informacionFuerzas] = btkGetForcePlatforms(h);
antropometria = btkFindMetaData(h,'Antropometria');
Eventos=btkGetEvents(h);

% ===========================================================================================================================================================%
%      DATOS ANTROPOMETRICOS

% Las medidas estan en cm entonces se usa un factor para pasarlas a metros
factor = 1/100;
%  A1 es el peso de la persona.
A1 = antropometria.children.PESO.info.values';

%  A2 es la distacia entre crestas il�acas.
A2 = antropometria.children.LONGITUD_ASIS.info.values*factor;

% A7 es la longitud de la pierna derecha.
A7 =antropometria.children.LONGITUD_PIERNA_DERECHA.info.values*factor;

% A8 es la longitud de la pierna derecha.
A8 =antropometria.children.LONGITUD_PIERNA_IZQUIERDA.info.values*factor;

% A11 es el di�metro de la rodilla derecha.
A11 =antropometria.children.DIAMETRO_RODILLA_DERECHA.info.values*factor;

% A12 es el di�metro de la rodilla izquierda.
A12 =antropometria.children.DIAMETRO_RODILLA_IZQUIERDA.info.values*factor;

% A13 es la longitud del pie derecho.
A13 =antropometria.children.LONGITUD_PIE_DERECHO.info.values*factor;

% A14 es la longitud del pie izquierdo.
A14 =antropometria.children.LONGITUD_PIE_IZQUIERDO.info.values*factor;

% A15 es la altura del maleolo derecho.
A15 =antropometria.children.ALTURA_MALEOLOS_DERECHO.info.values*factor;

% A16 es la altura del maleolo derecho.
A16 =antropometria.children.ALTURA_MALEOLOS_IZQUIERDO.info.values*factor;

% A17 es el ancho del maleolo derecho.
A17 =antropometria.children.ANCHO_MALEOLOS_DERECHO.info.values*factor;

% A18 es el ancho del maleolo izquierdo.
A18 =antropometria.children.ANCHO_MALEOLOS_IZQUIERDO.info.values*factor;

% A19 es el ancho del pie derecho.
A19 =antropometria.children.ANCHO_PIE_DERECHO.info.values*factor;

% A20 es el ancho del pie izquierdo.
A20 =antropometria.children.ANCHO_PIE_IZQUIERDO.info.values*factor;

% ===========================================================================================================================================================%
%  RECORTE DE LOS MARCADORES

% Se lee la frecuencia de muestreo de la camara
fm = informacionCine.frequency(1);

% Se encuentran los valores de muestras en los que comienza a registrarse la marcha. Con estos valores se ingresan a los vectores de marcadores
% y se recortan para quedarnos con la informacion util.

% Pie derecho (RHS = right heel strike)
RHS1 = round(Eventos.Derecho_RHS(1)*fm);
RHS2 = round(Eventos.Derecho_RHS(2)*fm);
% Pie izquierdo (LHS = left heel strike)
LHS1 = round(Eventos.Izquierdo_LHS(1)*fm);
LHS2 = round(Eventos.Izquierdo_LHS(2)*fm);
% En este caso el comienzo de la marcha se da con el apoyo del pie derecho
inicio = RHS1;
fin= LHS2;

% Se hace el recorte. La sintaxis es (fila_inicial:fila_final,columna inicial:columna final) si no se pone nada se toman todas las columnas.

% Cabeza del 2do metatarsiano derecho.
p1=marcadores.r_met(inicio:fin,:);

% Talon derecho.
p2=marcadores.r_heel(inicio:fin,:);

% Maleolo medial derecho.
p3=marcadores.r_mall(inicio:fin,:);

% Banda tibial derecha.
p4=marcadores.r_bar_2(inicio:fin,:);

% Epicondilo femoral derecho.
p5=marcadores.r_knee_1(inicio:fin,:);

% Banda femoral derecha.
p6=marcadores.r_bar_1(inicio:fin,:);

% Espina iliaca anterior superior derecha.
p7=marcadores.r_asis(inicio:fin,:);

% Cabeza del 2do metatarsino izquierdo.
p8=marcadores.l_met(inicio:fin,:);

% Talon izquierdo.
p9=marcadores.l_heel(inicio:fin,:);

% Maleolo lateral izquierdo.
p10=marcadores.l_mall(inicio:fin,:);

% Banda tibial izquierda.
p11=marcadores.l_bar_2(inicio:fin,:);

% Epicondilo femoral izquierdo.
p12=marcadores.l_knee_2(inicio:fin,:);

% Banda femoral izquierda.
p13=marcadores.l_bar_1(inicio:fin,:);

% Espina iliaca anterior superior izquierda.
p14=marcadores.l_asis(inicio:fin,:);

% Sacro.
p15=marcadores.sacrum(inicio:fin,:);

% ===========================================================================================================================================================%
%  FILTRADO DE LOS MARCADORES

% Frecuencia de Nyquist.
fn= fm/2;
% Frecuencia de corte del filtro.
fc = 10;
% b y a son el numerador y denomindor de la funcion de transferencia del filtro.
[b,a] = butter(2,fc/fn);
% Ahora filtramos cada uno de los marcadores por columna.
for i=1:3
    p1(:,i)    =  filtfilt(b,a, p1(:,i));
    p2(:,i)    =  filtfilt(b,a, p2(:,i));
    p3(:,i)    =  filtfilt(b,a, p3(:,i));
    p4(:,i)    =  filtfilt(b,a, p4(:,i));
    p5(:,i)    =  filtfilt(b,a, p5(:,i));
    p6(:,i)    =  filtfilt(b,a, p6(:,i));
    p7(:,i)    =  filtfilt(b,a, p7(:,i));
    p8(:,i)    =  filtfilt(b,a, p8(:,i));
    p9(:,i)    =  filtfilt(b,a, p9(:,i));
    p10(:,i)  =  filtfilt(b,a, p10(:,i));
    p11(:,i)  =  filtfilt(b,a, p11(:,i));
    p12(:,i)  =  filtfilt(b,a, p12(:,i));
    p13(:,i)  =  filtfilt(b,a, p13(:,i));
    p14(:,i)  =  filtfilt(b,a, p14(:,i));
    p15(:,i)  =  filtfilt(b,a, p15(:,i));
end

% Plot de testeo de los marcadores de los tobillos
%figure
%plot(p2(:,3),'LineWidth', 1.5);
%hold on
%plot(p9(:,3),'LineWidth', 1.5);
%legend({'Talon derecho','Talon izquierdo'})

% ===========================================================================================================================================================%
%  SE CALCULAN U,V,W 

%     Calculo para la PELVIS
% ------------   vPelvis = (p14 - p7)/|p14 - p7| 
numerador      = p14-p7;
denominador  = sqrt(sum(numerador.^2,2));
vPelvis               = numerador./denominador;

% ------------   wPelvis = (p7 - p15) x( p14-p15)/|(p7 - p15) x( p14-p15)| 
numerador      = cross((p7-p15),(p14-p15));
denominador  = sqrt(sum(numerador.^2,2));
wPelvis             = numerador./denominador;

% ------------   uPelvis = vPelvis x wPelvis
uPelvis             = cross(vPelvis,wPelvis);

% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%     Calculo para la PIERNA DERECHA
% ------------   vPierna_derecha =  (p3 - p5)/|p3 - p5| 
numerador           = p3-p5;
denominador       = sqrt(sum(numerador.^2,2));
vPierna_derecha  = numerador./denominador;

% ------------   uPierna_derecha = (p4 - p5) x( p3-p5)/|(p4 - p5) x( p3-p5)| 
numerador            = cross((p4-p5),(p3-p5));
denominador        = sqrt(sum(numerador.^2,2));
uPierna_derecha  = numerador./denominador;

% ------------   wPierna_derecha = uPierna_derecha x vPierna_derecha
wPierna_derecha   = cross(uPierna_derecha,vPierna_derecha);
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%

%     Calculo para la PIERNA IZQUIERDA
% ------------   vPierna_izquierda =  (p10 - p12)/|p10 - p12| 
numerador             = p10-p12;
denominador         = sqrt(sum(numerador.^2,2));
vPierna_izquierda  = numerador./denominador;

% ------------   uPierna_izquierda = (p11 - p12) x( p10-p12)/|(p11 - p12) x( p10-p12)| 
numerador              = cross((p11-p12),(p10-p12));
denominador          = sqrt(sum(numerador.^2,2));
uPierna_izquierda   = numerador./denominador;

% ------------   wPierna_izquierda = uPierna_izquierda x vPierna_izquierda
wPierna_izquierda   = cross(uPierna_izquierda,vPierna_izquierda);
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%

%     Calculo para el PIE DERECHO
% ------------   uPie_derecho =  (p1 - p2)/|p1 - p2| 
numerador           = p1-p3;
denominador       = sqrt(sum(numerador.^2,2));
uPie_derecho       = numerador./denominador;

% ------------   wPie_derecho = (p1 - p3) x( p2-p3)/|(p1 - p3) x( p2-p3)| 
numerador            = cross((p1-p3),(p2-p3));
denominador        = sqrt(sum(numerador.^2,2));
wPie_derecho       = numerador./denominador;

% ------------   vPie_derecho = wPie_derecho x uPie_derecho
vPie_derecho       = cross(wPie_derecho,uPie_derecho);
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%     Calculo para el PIE IZQUIERDO
% ------------   uPie_izquierdo =  (p8 - p9)/|p8 - p9| 
numerador           = p8-p9;
denominador       = sqrt(sum(numerador.^2,2));
uPie_izquierdo   = numerador./denominador;

% ------------   wPie_izquierdo = (p8 - p10) x( p9-p10)/|(p8- p10) x( p9-p10)| 
numerador            = cross((p8-p10),(p9-p10));
denominador        = sqrt(sum(numerador.^2,2));
wPie_izquierdo      = numerador./denominador;

% ------------   vPie_izquierdo = wPie_izquierdo x uPie_izquierdo
vPie_izquierdo       = cross(wPie_izquierdo,uPie_izquierdo);

% ===========================================================================================================================================================%
%  SE CALCULAN LAS POSICIONES ARTICULARES

%     Calculo para la CADERA DERECHA
  p_cadera_derecha  = p15+0.598*A2*uPelvis - 0.344*A2*vPelvis - 0.29*A2*wPelvis;
  
  %     Calculo para la CADERA IZQUIERDA
  p_cadera_izquierda  = p15+0.598*A2*uPelvis + 0.344*A2*vPelvis - 0.29*A2*wPelvis;
  
  %     Calculo para la RODILLA DERECHA
  p_rodilla_derecha  = p5+0.0*A11*uPierna_derecha+0.0*A11*vPierna_derecha + 0.5*A11*wPierna_derecha;
  
  %     Calculo para la RODILLA IZQUIERDA
  p_rodilla_izquierda  = p12+0.0*A12*uPierna_izquierda+0.0*A12*vPierna_izquierda - 0.5*A12*wPierna_izquierda;

  %     Calculo para el TOBILLO DERECHO
  p_tobillo_derecho  = p3+0.016*A13*uPie_derecho + 0.392*A15*vPie_derecho+0.478*A17*wPie_derecho;
  
  %     Calculo para el TOBILLO IZQUIERDO
  p_tobillo_izquierdo  = p10+0.016*A14*uPie_izquierdo+0.392*A16*vPie_izquierdo-0.478*A18*wPie_izquierdo;

  % Calculo para el DEDO DERECHO
  p_dedo_derecho = p3+0.742*A13*uPie_derecho + 1.074*A15*vPie_derecho - 0.187*A19*wPie_derecho;
  
    % Calculo para el DEDO IZQUIERDO
  p_dedo_izquierdo = p10+0.742*A14*uPie_izquierdo + 1.074*A16*vPie_izquierdo + 0.187*A20*wPie_izquierdo;
  
 % ===========================================================================================================================================================%
  
 figure
 plot3(p15(:,1),p15(:,2),p15(:,3),'LineWidth',2,'color',negro); 
 hold on
 plot3(p7(:,1),p7(:,2),p7(:,3),'LineWidth',2); 
 hold on
 plot3(p14(:,1),p14(:,2),p14(:,3),'LineWidth',2); 
 xlabel('x [m]','Interpreter','latex');
 ylabel('y [m]','Interpreter','latex');
 zlabel('z [m]','Interpreter','latex');
 hold on
 grid on
 for i=1:10:length(p15)
     quiver3(p15(i,1),p15(i,2),p15(i,3),uPelvis(i,1)*factor_vectores,uPelvis(i,2)*factor_vectores,uPelvis(i,3)*factor_vectores,'color',rojo,'LineWidth',1.5)
     quiver3(p15(i,1),p15(i,2),p15(i,3),vPelvis(i,1)*factor_vectores,vPelvis(i,2)*factor_vectores,vPelvis(i,3)*factor_vectores,'color',verde,'LineWidth',1.5)
     quiver3(p15(i,1),p15(i,2),p15(i,3),wPelvis(i,1)*factor_vectores,wPelvis(i,2).*factor_vectores,wPelvis(i,3)*factor_vectores,'color',azul,'LineWidth',1.5)
 end
 legend({'Marcador sacro','Marcador espina iliaca derecha','Marcador espina iliaca izquierda','$\vec{u}$','$\vec{v}$','$\vec{w}$'},'Interpreter','latex')