%% Limpiar variables
clear all;

%%  Variables extras.
rojo                          = [1 0 0];
verde                       = [0.051 0.529 0.051 ];
azul                          = [0 0 1];
negro                       = [0 0 0];
purpura                   = [1, 0.11, 0.612, 0.2];
factor_vectores     = 1/30;
paso_vectores       = 10;
line_width               = 1.75;
factor_a_grados    = 180/pi;

Biomecanica         = FuncionesBiomecanica;
% ===========================================================================================================================================================%
%%  Se importan los datos de la marcha.
ruta_archivo = '0029_Davis_Marcha01_Walking10b2021.c3d';
h = btkReadAcquisition(ruta_archivo);
[marcadores, informacionCine] = btkGetMarkers(h);
[fuerzas, informacionFuerzas] = btkGetForcePlatforms(h);
antropometria = btkFindMetaData(h,'Antropometria');
Eventos=btkGetEvents(h);
% ===========================================================================================================================================================%
%%  Datos antropometricos.
% Las medidas estan en cm entonces se usa un factor para pasarlas a metros
factor = 1/100;
%  A1 es el peso de la persona.
A1 = antropometria.children.PESO.info.values';

%  A2 es la distacia entre crestas ilíacas.
A2 = antropometria.children.LONGITUD_ASIS.info.values*factor;

% A7 es la longitud de la pierna derecha.
A7 =antropometria.children.LONGITUD_PIERNA_DERECHA.info.values*factor;

% A8 es la longitud de la pierna derecha.
A8 =antropometria.children.LONGITUD_PIERNA_IZQUIERDA.info.values*factor;

% A11 es el diámetro de la rodilla derecha.
A11 =antropometria.children.DIAMETRO_RODILLA_DERECHA.info.values*factor;

% A12 es el diámetro de la rodilla izquierda.
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
%%  Recorte de los marcadores.

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
%%  Filtrado de los marcadores.

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


% ===========================================================================================================================================================%
%%  Se calculan u,v,w.

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
%     Calculo para la PIERNA DERECHA (rodilla)
% ------------   vPierna_derecha =  (p3 - p5)/|p3 - p5|
numerador           = p3-p5;
denominador       = sqrt(sum(numerador.^2,2));
vRodilla_derecha  = numerador./denominador;

% ------------   uRodilla_derecha = (p4 - p5) x( p3-p5)/|(p4 - p5) x( p3-p5)|
numerador            = cross((p4-p5),(p3-p5));
denominador        = sqrt(sum(numerador.^2,2));
uRodilla_derecha  = numerador./denominador;

% ------------   wRodilla_derecha = uRodilla_derecha x vRodilla_derecha
wRodilla_derecha   = cross(uRodilla_derecha,vRodilla_derecha);
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%

%     Calculo para la RODILLA IZQUIERDA
% ------------   vPierna_izquierda =  (p10 - p12)/|p10 - p12|
numerador             = p10-p12;
denominador         = sqrt(sum(numerador.^2,2));
vRodilla_izquierda  = numerador./denominador;

% ------------   uRodilla_izquierda = (p10 - p12) x( p11-p12)/|(p10 - p12) x( p11-p12)|
numerador              = cross((p10-p12),(p11-p12));
denominador          = sqrt(sum(numerador.^2,2));
uRodilla_izquierda   = numerador./denominador;

% ------------   wRodilla_izquierda = uRodilla_izquierda x vRodilla_izquierda
wRodilla_izquierda   = cross(uRodilla_izquierda,vRodilla_izquierda);
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
%%  Calculo de las posiciones articulares.

%     Calculo para la CADERA DERECHA
p_cadera_derecha  = p15+0.598*A2*uPelvis - 0.344*A2*vPelvis - 0.29*A2*wPelvis;

%     Calculo para la CADERA IZQUIERDA
p_cadera_izquierda  = p15+0.598*A2*uPelvis + 0.344*A2*vPelvis - 0.29*A2*wPelvis;

%     Calculo para la RODILLA DERECHA
p_rodilla_derecha  = p5+0.0*A11*uRodilla_derecha+0.0*A11*vRodilla_derecha + 0.5*A11*wRodilla_derecha;

%     Calculo para la RODILLA IZQUIERDA
p_rodilla_izquierda  = p12+0.0*A12*uRodilla_izquierda+0.0*A12*vRodilla_izquierda - 0.5*A12*wRodilla_izquierda;

%     Calculo para el TOBILLO DERECHO
p_tobillo_derecho  = p3+0.016*A13*uPie_derecho + 0.392*A15*vPie_derecho+0.478*A17*wPie_derecho;

%     Calculo para el TOBILLO IZQUIERDO
p_tobillo_izquierdo  = p10+0.016*A14*uPie_izquierdo+0.392*A16*vPie_izquierdo-0.478*A18*wPie_izquierdo;

% Calculo para el DEDO DERECHO
p_dedo_derecho = p3+0.742*A13*uPie_derecho + 1.074*A15*vPie_derecho - 0.187*A19*wPie_derecho;

% Calculo para el DEDO IZQUIERDO
p_dedo_izquierdo = p10+0.742*A14*uPie_izquierdo + 1.074*A16*vPie_izquierdo + 0.187*A20*wPie_izquierdo;

% ===========================================================================================================================================================%

%%  Plot de u,v,w y las posiciones articulares.
close all % cierra las ventanas
clc % limpia la consola

    figure
    % Posicion del marcador del sacro como linea
    plot3(p15(:,1),p15(:,2),p15(:,3),'LineWidth',0.75,'color',negro,'HandleVisibility','off');
    hold on
    % Posicion de la cadera derecha como linea
    plot3(p_cadera_derecha(:,1),p_cadera_derecha(:,2),p_cadera_derecha(:,3),'LineWidth',0.75,'color',negro,'HandleVisibility','off');
    hold on
    % Posicion de la cadera izquierda como linea
    plot3(p_cadera_izquierda(:,1),p_cadera_izquierda(:,2),p_cadera_izquierda(:,3),'LineWidth',0.75,'color',negro,'HandleVisibility','off');
    hold on
    % Posicion de la rodilla derecha como linea
    plot3(p_rodilla_derecha(:,1),p_rodilla_derecha(:,2),p_rodilla_derecha(:,3),'LineWidth',0.75,'color',negro,'HandleVisibility','off');
    hold on
    % Posicion de la rodilla izquierda como linea
    plot3(p_rodilla_izquierda(:,1),p_rodilla_izquierda(:,2),p_rodilla_izquierda(:,3),'LineWidth',0.75,'color',negro,'HandleVisibility','off');
    hold on
    % Posicion del tobillo derecho como linea
    plot3(p_tobillo_derecho(:,1),p_tobillo_derecho(:,2),p_tobillo_derecho(:,3),'LineWidth',0.75,'color',negro,'HandleVisibility','off');
    hold on
    % Posicion del tobillo izquierdo como linea
    plot3(p_tobillo_izquierdo(:,1),p_tobillo_izquierdo(:,2),p_tobillo_izquierdo(:,3),'LineWidth',0.75,'color',negro,'HandleVisibility','off');
    hold on
    for i=1:paso_vectores:length(p_cadera_derecha)
        
        % Marcadores del sacro
        plot3(p15(i,1),p15(i,2),p15(i,3),'diamond','LineWidth',line_width,'color',negro);
        hold on
        % Marcadores de la cadera derecha
        plot3(p_cadera_derecha(i,1),p_cadera_derecha(i,2),p_cadera_derecha(i,3),'o','LineWidth',line_width,'color',negro);
        hold on
        % Marcadores de la cadera izquierda
        plot3(p_cadera_izquierda(i,1),p_cadera_izquierda(i,2),p_cadera_izquierda(i,3),'square','LineWidth',line_width,'color',negro);
        hold on
        % Marcadores de la rodilla derecha
        plot3(p_rodilla_derecha(i,1),p_rodilla_derecha(i,2),p_rodilla_derecha(i,3),'*','LineWidth',line_width,'color',negro);
        hold on
        % Marcadores de la rodilla izquierda
        plot3(p_rodilla_izquierda(i,1),p_rodilla_izquierda(i,2),p_rodilla_izquierda(i,3),'x','LineWidth',line_width,'color',negro);
        hold on
        % Marcadores del tobillo derecho
        plot3(p_tobillo_derecho(i,1),p_tobillo_derecho(i,2),p_tobillo_derecho(i,3),'pentagram','LineWidth',line_width,'color',negro);
        hold on
        % Marcadores del tobillo izquierdo
        plot3(p_tobillo_izquierdo(i,1),p_tobillo_izquierdo(i,2),p_tobillo_izquierdo(i,3),'>','LineWidth',line_width,'color',negro);
        hold on
        
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
        % CADERA DERECHA
        % Vector u
        quiver3(p_cadera_derecha(i,1),...
            p_cadera_derecha(i,2),...
            p_cadera_derecha(i,3),...
            uPelvis(i,1)*factor_vectores,uPelvis(i,2)*factor_vectores,uPelvis(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
        % Vector v
        quiver3(p_cadera_derecha(i,1),...
            p_cadera_derecha(i,2),...
            p_cadera_derecha(i,3),...
            vPelvis(i,1)*factor_vectores,vPelvis(i,2)*factor_vectores,vPelvis(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
        % Vector w
        quiver3(p_cadera_derecha(i,1),...
            p_cadera_derecha(i,2),...
            p_cadera_derecha(i,3),...
            wPelvis(i,1)*factor_vectores,wPelvis(i,2).*factor_vectores,wPelvis(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
        % CADERA IZQUIERDA
        % Vector u
        quiver3(p_cadera_izquierda(i,1),...
            p_cadera_izquierda(i,2),...
            p_cadera_izquierda(i,3),...
            uPelvis(i,1)*factor_vectores,uPelvis(i,2)*factor_vectores,uPelvis(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
        % Vector v
        quiver3(p_cadera_izquierda(i,1),...
            p_cadera_izquierda(i,2),...
            p_cadera_izquierda(i,3),vPelvis(i,1)*factor_vectores,vPelvis(i,2)*factor_vectores,vPelvis(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
        % Vector w
        quiver3(p_cadera_izquierda(i,1),...
            p_cadera_izquierda(i,2),...
            p_cadera_izquierda(i,3),...
            wPelvis(i,1)*factor_vectores,wPelvis(i,2).*factor_vectores,wPelvis(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
        % RODILLA DERECHA
        % Vector u
        quiver3(p_rodilla_derecha(i,1),...
            p_rodilla_derecha(i,2),...
            p_rodilla_derecha(i,3),...
            uRodilla_derecha(i,1)*factor_vectores,uRodilla_derecha(i,2)*factor_vectores,uRodilla_derecha(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
        % Vector v
        quiver3(p_rodilla_derecha(i,1),...
            p_rodilla_derecha(i,2),...
            p_rodilla_derecha(i,3),...
            vRodilla_derecha(i,1)*factor_vectores,vRodilla_derecha(i,2)*factor_vectores,vRodilla_derecha(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
        % Vector w
        quiver3(p_rodilla_derecha(i,1),...
            p_rodilla_derecha(i,2),...
            p_rodilla_derecha(i,3),...
            wRodilla_derecha(i,1)*factor_vectores,wRodilla_derecha(i,2).*factor_vectores,wRodilla_derecha(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
        % RODILLA IZQUIERDA
        % Vector u
        quiver3(p_rodilla_izquierda(i,1),...
            p_rodilla_izquierda(i,2),...
            p_rodilla_izquierda(i,3),...
            uRodilla_izquierda(i,1)*factor_vectores,uRodilla_izquierda(i,2)*factor_vectores,uRodilla_izquierda(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
        % Vector v
        quiver3(p_rodilla_izquierda(i,1),...
            p_rodilla_izquierda(i,2),...
            p_rodilla_izquierda(i,3),...
            vRodilla_izquierda(i,1)*factor_vectores,uRodilla_izquierda(i,2)*factor_vectores,uRodilla_izquierda(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
        % Vector w
        quiver3(p_rodilla_izquierda(i,1),...
            p_rodilla_izquierda(i,2),...
            p_rodilla_izquierda(i,3),...
            wRodilla_izquierda(i,1)*factor_vectores,wRodilla_izquierda(i,2).*factor_vectores,wRodilla_izquierda(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
        % TOBILLO DERECHO
        % Vector u        
        quiver3(p_tobillo_derecho(i,1),...
            p_tobillo_derecho(i,2),...
            p_tobillo_derecho(i,3),...
            uPie_derecho(i,1)*factor_vectores,uPie_derecho(i,2)*factor_vectores,uPie_derecho(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
        % Vector v
        quiver3(p_tobillo_derecho(i,1),...
            p_tobillo_derecho(i,2),...
            p_tobillo_derecho(i,3),...
            vPie_derecho(i,1)*factor_vectores,vPie_derecho(i,2)*factor_vectores,vPie_derecho(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
        % Vector w
        quiver3(p_tobillo_derecho(i,1),...
            p_tobillo_derecho(i,2),...
            p_tobillo_derecho(i,3),wPie_derecho(i,1)*factor_vectores,wPie_derecho(i,2).*factor_vectores,wPie_derecho(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
        % TOBILLO IZQUIERDO
        % Vector u
        quiver3(p_tobillo_izquierdo(i,1),...
            p_tobillo_izquierdo(i,2),...
            p_tobillo_izquierdo(i,3),...
            uPie_izquierdo(i,1)*factor_vectores,uPie_izquierdo(i,2)*factor_vectores,uPie_izquierdo(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
        % Vector v
        quiver3(p_tobillo_izquierdo(i,1),...
            p_tobillo_izquierdo(i,2),...
            p_tobillo_izquierdo(i,3),...
            vPie_izquierdo(i,1)*factor_vectores,vPie_izquierdo(i,2)*factor_vectores,vPie_izquierdo(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
        % Vector w
        quiver3(p_tobillo_izquierdo(i,1),...
            p_tobillo_izquierdo(i,2),...
            p_tobillo_izquierdo(i,3),...
            wPie_izquierdo(i,1)*factor_vectores,wPie_izquierdo(i,2).*factor_vectores,wPie_izquierdo(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
        
    end
    grid on
    title('Trayectoria de los segmentos articulares junto con los vectores $\vec{u}$, $\vec{v}$ y $\vec{w}$','Interpreter','latex')
    xlabel('x [m]','Interpreter','latex');
    ylabel('y [m]','Interpreter','latex');
    zlabel('z [m]','Interpreter','latex');
    legend({'Marcador sacro (p15)',...,
        'Posicion cadera derecha ',...,
        'Posicion cadera izquierda ',...,
        'Posicion rodilla derecha ',...,
        'Posicion rodilla izquierda ',...,
        'Posicion tobillo derecho ',...,
        'Posicion tobillo izquierdo ',...,
        '$\vec{u}$',...,
        '$\vec{v}$',...,
        '$\vec{w}$'},'Interpreter','latex')
%%  Plot u,v,w pelvis
close all % cierra las ventanas
clc % limpia la consola

    figure
    % Posicion de la cadera derecha como linea
    plot3(p_cadera_derecha(:,1),p_cadera_derecha(:,2),p_cadera_derecha(:,3),'LineWidth',line_width);
    hold on
    % Posicion de la cadera izquierda como linea
    plot3(p_cadera_izquierda(:,1),p_cadera_izquierda(:,2),p_cadera_izquierda(:,3),'LineWidth',line_width);
    hold on
    for i=1:paso_vectores:length(p_cadera_derecha)       
       %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
        % CADERA DERECHA
        % Vector u
        quiver3(p_cadera_derecha(i,1),...
            p_cadera_derecha(i,2),...
            p_cadera_derecha(i,3),...
            uPelvis(i,1)*factor_vectores,uPelvis(i,2)*factor_vectores,uPelvis(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
        hold on
        % Vector v
        quiver3(p_cadera_derecha(i,1),...
            p_cadera_derecha(i,2),...
            p_cadera_derecha(i,3),...
            vPelvis(i,1)*factor_vectores,vPelvis(i,2)*factor_vectores,vPelvis(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
        hold on
        % Vector w
        quiver3(p_cadera_derecha(i,1),...
            p_cadera_derecha(i,2),...
            p_cadera_derecha(i,3),...
            wPelvis(i,1)*factor_vectores,wPelvis(i,2).*factor_vectores,wPelvis(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
        hold on
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
        % CADERA IZQUIERDA
        % Vector u
        quiver3(p_cadera_izquierda(i,1),...
            p_cadera_izquierda(i,2),...
            p_cadera_izquierda(i,3),...
            uPelvis(i,1)*factor_vectores,uPelvis(i,2)*factor_vectores,uPelvis(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
        hold on
        % Vector v
        quiver3(p_cadera_izquierda(i,1),...
            p_cadera_izquierda(i,2),...
            p_cadera_izquierda(i,3),vPelvis(i,1)*factor_vectores,vPelvis(i,2)*factor_vectores,vPelvis(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
        hold on
        % Vector w
        quiver3(p_cadera_izquierda(i,1),...
            p_cadera_izquierda(i,2),...
            p_cadera_izquierda(i,3),...
            wPelvis(i,1)*factor_vectores,wPelvis(i,2).*factor_vectores,wPelvis(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    end
    grid on
    title('Posiciones articulares (PA) de las caderas y vectores $\vec{u}$,$\vec{v}$,$\vec{w}$','Interpreter','latex')
    xlabel('x [m]','Interpreter','latex');
    ylabel('y [m]','Interpreter','latex');
    zlabel('z [m]','Interpreter','latex');
    legend({'PA cadera derecha',...,
        'PA cadera izquierda',...,
        '$\vec{u}$',...,
        '$\vec{v}$',...,
        '$\vec{w}$'},'Interpreter','latex')
%%  Plot u,v,w rodillas
close all % cierra las ventanas
clc % limpia la consola

    figure
    % Posicion de la rodilla derecha como linea
    plot3(p_rodilla_derecha(:,1),p_rodilla_derecha(:,2),p_rodilla_derecha(:,3),'LineWidth',line_width);
    hold on
    % Posicion de la rodilla izquierda como linea
    plot3(p_rodilla_izquierda(:,1),p_rodilla_izquierda(:,2),p_rodilla_izquierda(:,3),'LineWidth',line_width);
    hold on
    for i=1:paso_vectores:length(p_cadera_derecha)       
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
        % RODILLA DERECHA
        % Vector u
        quiver3(p_rodilla_derecha(i,1),...
            p_rodilla_derecha(i,2),...
            p_rodilla_derecha(i,3),...
            uRodilla_derecha(i,1)*factor_vectores,uRodilla_derecha(i,2)*factor_vectores,uRodilla_derecha(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
        % Vector v
        quiver3(p_rodilla_derecha(i,1),...
            p_rodilla_derecha(i,2),...
            p_rodilla_derecha(i,3),...
            vRodilla_derecha(i,1)*factor_vectores,vRodilla_derecha(i,2)*factor_vectores,vRodilla_derecha(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
        % Vector w
        quiver3(p_rodilla_derecha(i,1),...
            p_rodilla_derecha(i,2),...
            p_rodilla_derecha(i,3),...
            wRodilla_derecha(i,1)*factor_vectores,wRodilla_derecha(i,2).*factor_vectores,wRodilla_derecha(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
        % RODILLA IZQUIERDA
        % Vector u
        quiver3(p_rodilla_izquierda(i,1),...
            p_rodilla_izquierda(i,2),...
            p_rodilla_izquierda(i,3),...
            uRodilla_izquierda(i,1)*factor_vectores,uRodilla_izquierda(i,2)*factor_vectores,uRodilla_izquierda(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
        % Vector v
        quiver3(p_rodilla_izquierda(i,1),...
            p_rodilla_izquierda(i,2),...
            p_rodilla_izquierda(i,3),...
            vRodilla_izquierda(i,1)*factor_vectores,uRodilla_izquierda(i,2)*factor_vectores,uRodilla_izquierda(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
        % Vector w
        quiver3(p_rodilla_izquierda(i,1),...
            p_rodilla_izquierda(i,2),...
            p_rodilla_izquierda(i,3),...
            wRodilla_izquierda(i,1)*factor_vectores,wRodilla_izquierda(i,2).*factor_vectores,wRodilla_izquierda(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    end
    grid on
    title('Posiciones articulares (PA) de las rodillas y $\vec{u}$, $\vec{v}$ y $\vec{w}$','Interpreter','latex')
    xlabel('x [m]','Interpreter','latex');
    ylabel('y [m]','Interpreter','latex');
    zlabel('z [m]','Interpreter','latex');
    legend({'PA rodilla derecha ',...,
        'PA rodilla izquierda ',...,
        '$\vec{u}$',...,
        '$\vec{v}$',...,
        '$\vec{w}$'},'Interpreter','latex') 
%%  Plot u,v,w tobillos
close all % cierra las ventanas
clc % limpia la consola

    figure
    % Posicion del tobillo derecho como linea
    plot3(p_tobillo_derecho(:,1),p_tobillo_derecho(:,2),p_tobillo_derecho(:,3),'LineWidth',line_width);
    hold on
    % Posicion del tobillo izquierdo como linea
    plot3(p_tobillo_izquierdo(:,1),p_tobillo_izquierdo(:,2),p_tobillo_izquierdo(:,3),'LineWidth',line_width);
    hold on
    for i=1:paso_vectores:length(p_cadera_derecha)       
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
        % TOBILLO DERECHO
        % Vector u        
        quiver3(p_tobillo_derecho(i,1),...
            p_tobillo_derecho(i,2),...
            p_tobillo_derecho(i,3),...
            uPie_derecho(i,1)*factor_vectores,uPie_derecho(i,2)*factor_vectores,uPie_derecho(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
        % Vector v
        quiver3(p_tobillo_derecho(i,1),...
            p_tobillo_derecho(i,2),...
            p_tobillo_derecho(i,3),...
            vPie_derecho(i,1)*factor_vectores,vPie_derecho(i,2)*factor_vectores,vPie_derecho(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
        % Vector w
        quiver3(p_tobillo_derecho(i,1),...
            p_tobillo_derecho(i,2),...
            p_tobillo_derecho(i,3),wPie_derecho(i,1)*factor_vectores,wPie_derecho(i,2).*factor_vectores,wPie_derecho(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
        % TOBILLO IZQUIERDO
        % Vector u
        quiver3(p_tobillo_izquierdo(i,1),...
            p_tobillo_izquierdo(i,2),...
            p_tobillo_izquierdo(i,3),...
            uPie_izquierdo(i,1)*factor_vectores,uPie_izquierdo(i,2)*factor_vectores,uPie_izquierdo(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
        % Vector v
        quiver3(p_tobillo_izquierdo(i,1),...
            p_tobillo_izquierdo(i,2),...
            p_tobillo_izquierdo(i,3),...
            vPie_izquierdo(i,1)*factor_vectores,vPie_izquierdo(i,2)*factor_vectores,vPie_izquierdo(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
        % Vector w
        quiver3(p_tobillo_izquierdo(i,1),...
            p_tobillo_izquierdo(i,2),...
            p_tobillo_izquierdo(i,3),...
            wPie_izquierdo(i,1)*factor_vectores,wPie_izquierdo(i,2).*factor_vectores,wPie_izquierdo(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
           end
    grid on
    title('Posiciones articulares (PA) de los tobillos y $\vec{u}$, $\vec{v}$ y $\vec{w}$','Interpreter','latex')
    xlabel('x [m]','Interpreter','latex');
    ylabel('y [m]','Interpreter','latex');
    zlabel('z [m]','Interpreter','latex');
    legend({'PA tobillo derecho ',...,
        'PA tobillo izquierda ',...,
        '$\vec{u}$',...,
        '$\vec{v}$',...,
        '$\vec{w}$'},'Interpreter','latex')
    
%%  Calculos de los centros de masa.
% Muslo derecho
centro_masa_muslo_derecho = p_cadera_derecha + 0.39*(p_rodilla_derecha - p_cadera_derecha);
% Muslo derecho
centro_masa_muslo_izquierdo = p_cadera_izquierda + 0.39*(p_rodilla_izquierda - p_cadera_izquierda);

% Pierna derecha
centro_masa_pierna_derecha = p_rodilla_derecha + 0.42*(p_tobillo_derecho - p_rodilla_derecha);
% Pierna izquierda
centro_masa_pierna_izquierda = p_rodilla_izquierda + 0.42*(p_tobillo_izquierdo - p_rodilla_izquierda);

% Pie derecho
centro_masa_pie_derecho = p_tobillo_derecho + 0.44*(p_dedo_derecho - p_tobillo_derecho);
% Pie izquierdo
centro_masa_pie_izquierdo = p_tobillo_izquierdo + 0.44*(p_dedo_izquierdo - p_tobillo_izquierdo);

% ===========================================================================================================================================================%
%%  Calculos de i,j,k.

%     Calculo para la PELVIS
% ------------
iPelvis  = wPelvis;
jPelvis   = uPelvis;
kPelvis = vPelvis;
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%     Calculo para el MUSLO DERECHO

% ------------   iMuslo_derecho  = ( pR.Hip - pR.Knee) / |pR.Hip - pR.Knee|
numerador           = p_cadera_derecha - p_rodilla_derecha;
denominador       = sqrt(sum(numerador.^2,2));
iMuslo_derecho   = numerador./denominador;

% ------------   jMuslo_derecho  = ( (p6 - pR.Hip) x (pR.Knee - pR.Hip) ) / |( (p6 - pR.Hip) x (pR.Knee - pR.Hip) )|
numerador            = cross(p6-p_cadera_derecha, p_rodilla_derecha - p_cadera_derecha);
denominador        = sqrt(sum(numerador.^2,2));
jMuslo_derecho    = numerador./denominador;

% ------------   kMuslo_derecho  = iMuslo_derecho x jMuslo_derecho
kMuslo_derecho    = cross(iMuslo_derecho, jMuslo_derecho);

%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%     Calculo para el MUSLO IZQUIERDO

% ------------   iMuslo_izquierdo  = (pL.Hip - pL.Knee)/ |pL.Hip - pL.Knee|
numerador           = p_cadera_izquierda - p_rodilla_izquierda;
denominador       = sqrt(sum(numerador.^2,2));
iMuslo_izquierdo   = numerador./denominador;

% ------------   jMuslo_izquierdo  = ( (pL.Knee - pL.Hip) x (p13 - pL.Hip) ) / | (pL.Knee - pL.Hip) x (p13 - pL.Hip)|
numerador            = cross(p_rodilla_izquierda-p_cadera_izquierda, p13 - p_cadera_izquierda);
denominador        = sqrt(sum(numerador.^2,2));
jMuslo_izquierdo    = numerador./denominador;

% ------------   kMuslo_izquierdo  = iMuslo_izquierdo x jMuslo_izquierdo
kMuslo_izquierdo    = cross(iMuslo_izquierdo, jMuslo_izquierdo);

%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%     Calculo para la PIERNA DERECHA

% ------------   iPierna_derecha  = (pR.Knee - pR.Ankle) / |pR.Knee - pR.Ankle|
numerador           = p_rodilla_derecha - p_tobillo_derecho;
denominador       = sqrt(sum(numerador.^2,2));
iPierna_derecha   = numerador./denominador;

% ------------   jPierna_derecha  = ( (p5 - pR.Knee) x (pR.Ankle - pR.Knee) ) / |(p5 - pR.Knee) x (pR.Ankle - pR.Knee)|
numerador            = cross(p5-p_rodilla_derecha, p_tobillo_derecho - p_rodilla_derecha);
denominador        = sqrt(sum(numerador.^2,2));
jPierna_derecha    = numerador./denominador;

% ------------   kPierna_derecha  = iPierna_derecha x jPierna_derecha
kPierna_derecha    = cross(iPierna_derecha, jPierna_derecha);

%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%     Calculo para la PIERNA IZQUIERDA

% ------------   iPierna_izquierda  = (pL.Knee - pL.Ankle) / |pL.Knee - pL.Ankle|
numerador           = p_rodilla_izquierda - p_tobillo_izquierdo;
denominador       = sqrt(sum(numerador.^2,2));
iPierna_izquierda   = numerador./denominador;

% ------------   jPierna_derecha  = ( (pL.Ankle - pL.Knee) x (p12 - pL.Knee)) / |(pL.Ankle - pL.Knee) x (p12 - pL.Knee)|
numerador            = cross(p_tobillo_izquierdo-p_rodilla_izquierda, p12 - p_rodilla_izquierda);
denominador        = sqrt(sum(numerador.^2,2));
jPierna_izquierda    = numerador./denominador;

% ------------   kPierna_derecha  = iPierna_derecha x jPierna_derecha
kPierna_izquierda    = cross(iPierna_izquierda, jPierna_izquierda);

%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%     Calculo para el PIE DERECHO

% ------------   iPie_derecho  = (p2 - pR.Toe)/ |p2 - pR.Toe|
numerador           = p2 - p_dedo_derecho;
denominador       = sqrt(sum(numerador.^2,2));
iPie_derecho   = numerador./denominador;

% ------------   kPie_derecho  = ( (pR.Ankle - p2) x (pR.Toe - p2) ) / |(pR.Ankle - p2) x (pR.Toe - p2) |
numerador            = cross(p_tobillo_derecho-p2, p_dedo_derecho - p2);
denominador        = sqrt(sum(numerador.^2,2));
kPie_derecho    = numerador./denominador;

% ------------   jPie_derecho  = kPie_derecho x iPie_derecho
jPie_derecho    = cross(kPie_derecho, iPie_derecho);

%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%     Calculo para el PIE IZQUIERDO

% ------------   iPie_izquierdo  = (p9 - pL.Toe) / |p9 - pL.Toe|
numerador           = p9 - p_dedo_izquierdo;
denominador       = sqrt(sum(numerador.^2,2));
iPie_izquierdo   = numerador./denominador;

% ------------   kPie_izquierdo  = ( (pL.Ankle - p9) x (pL.Toe - p9) ) / |(pL.Ankle - p9) x (pL.Toe - p9) |
numerador            = cross(p_tobillo_izquierdo-p9, p_dedo_izquierdo - p9);
denominador        = sqrt(sum(numerador.^2,2));
kPie_izquierdo    = numerador./denominador;

% ------------   jPie_izquierdo  = kPie_izquierdo x iPie_izquierdo
jPie_izquierdo    = cross(kPie_izquierdo, iPie_izquierdo);
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%

%%  Plot de i,j,k, posiciones articulares y los centros de masa.

close all % cierra las ventanas
clc % limpia la consola
figure
% Posicion articular de la cadera derecha
plot3(p_cadera_derecha(:,1),p_cadera_derecha(:,2),p_cadera_derecha(:,3),'LineWidth',line_width);
hold on
% Posicion articular de la cadera izquierda
plot3(p_cadera_izquierda(:,1),p_cadera_izquierda(:,2),p_cadera_izquierda(:,3),'LineWidth',line_width);
hold on
% Posicion articular de la rodilla derecha
plot3(p_rodilla_derecha(:,1),p_rodilla_derecha(:,2),p_rodilla_derecha(:,3),'LineWidth',line_width);
hold on
% Posicion articular de la rodilla izquierda
plot3(p_rodilla_izquierda(:,1),p_rodilla_izquierda(:,2),p_rodilla_izquierda(:,3),'LineWidth',line_width);
hold on
% Posicion articular del tobillo derecho
plot3(p_tobillo_derecho(:,1),p_tobillo_derecho(:,2),p_tobillo_derecho(:,3),'LineWidth',line_width);
hold on
% Posicion articular del tobillo izquierdo
plot3(p_tobillo_izquierdo(:,1),p_tobillo_izquierdo(:,2),p_tobillo_izquierdo(:,3),'LineWidth',line_width);
hold on
% Se hace un plot de los centros de masa
% Marcadores de la pelvis
plot3(p15(:,1),p15(:,2),p15(:,3),'LineWidth',0.75,'color',negro,'HandleVisibility','off');
hold on
% Marcadores del centro de masa del muslo derecho
plot3(centro_masa_muslo_derecho(:,1),centro_masa_muslo_derecho(:,2),centro_masa_muslo_derecho(:,3),'LineWidth',0.75,'color',negro,'HandleVisibility','off');
hold on
% Marcadores del centro de masa del muslo izquierdo
plot3(centro_masa_muslo_izquierdo(:,1),centro_masa_muslo_izquierdo(:,2),centro_masa_muslo_izquierdo(:,3),'LineWidth',0.75,'color',negro,'HandleVisibility','off');
hold on
% Marcadores del centro de masa de la pierna derecha
plot3(centro_masa_pierna_derecha(:,1),centro_masa_pierna_derecha(:,2),centro_masa_pierna_derecha(:,3),'LineWidth',0.75,'color',negro,'HandleVisibility','off');
hold on
% Marcadores del centro de masa de la pierna izquierda
plot3(centro_masa_pierna_izquierda(:,1),centro_masa_pierna_izquierda(:,2),centro_masa_pierna_izquierda(:,3),'LineWidth',0.75,'color',negro,'HandleVisibility','off');
hold on
% Marcadores del centro de masa del pie derecho
plot3(centro_masa_pie_derecho(:,1),centro_masa_pie_derecho(:,2),centro_masa_pie_derecho(:,3),'LineWidth',0.75,'color',negro,'HandleVisibility','off');
hold on
% Marcadores del centro de masa del pie izquierdo
plot3(centro_masa_pie_izquierdo(:,1),centro_masa_pie_izquierdo(:,2),centro_masa_pie_izquierdo(:,3),'LineWidth',0.75,'color',negro,'HandleVisibility','off');
hold on

% Se plotean los vectores
for i=1:paso_vectores:length(centro_masa_muslo_derecho)
    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    % Se plotean los marcadores donde luego se van a plotear los vectores
    %Marcadores de la pelvis
    plot3(p15(i,1),p15(i,2),p15(i,3),'diamond','LineWidth',line_width,'color',negro);
    hold on
    
    % Marcadores del centro de masa del muslo derecho
    plot3(centro_masa_muslo_derecho(i,1),centro_masa_muslo_derecho(i,2),centro_masa_muslo_derecho(i,3),'o','LineWidth',line_width,'color',negro);
    hold on
    
    % Marcadores del centro de masa del muslo izquierdo
    plot3(centro_masa_muslo_izquierdo(i,1),centro_masa_muslo_izquierdo(i,2),centro_masa_muslo_izquierdo(i,3),'square','LineWidth',line_width,'color',negro);
    hold on
    
    % Marcadores del centro de masa de la pierna derecha
    plot3(centro_masa_pierna_derecha(i,1),centro_masa_pierna_derecha(i,2),centro_masa_pierna_derecha(i,3),'*','LineWidth',line_width,'color',negro);
    hold on
    
    % Marcadores del centro de masa de la pierna izquierda
    plot3(centro_masa_pierna_izquierda(i,1),centro_masa_pierna_izquierda(i,2),centro_masa_pierna_izquierda(i,3),'x','LineWidth',line_width,'color',negro);
    hold on
    
    % Marcadores del centro de masa del pie derecho
    plot3(centro_masa_pie_derecho(i,1),centro_masa_pie_derecho(i,2),centro_masa_pie_derecho(i,3),'pentagram','LineWidth',line_width,'color',negro);
    hold on
    
    % Marcadores del centro de masa del pie izquierdo
    plot3(centro_masa_pie_izquierdo(i,1),centro_masa_pie_izquierdo(i,2),centro_masa_pie_izquierdo(i,3),'>','LineWidth',line_width,'color',negro);
    hold on
    
    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    % SE PLOTEAN LOS VECTORES i j k
    % i,j,k de la pelvis
    quiver3(p15(i,1),p15(i,2),p15(i,3),iPelvis(i,1)*factor_vectores,iPelvis(i,2)*factor_vectores,iPelvis(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
    quiver3(p15(i,1),p15(i,2),p15(i,3),jPelvis(i,1)*factor_vectores,jPelvis(i,2)*factor_vectores,jPelvis(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
    quiver3(p15(i,1),p15(i,2),p15(i,3),kPelvis(i,1)*factor_vectores,kPelvis(i,2).*factor_vectores,kPelvis(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
    
    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    % MUSLO DERECHO
    % i del muslo derecho
    quiver3(centro_masa_muslo_derecho(i,1),...
        centro_masa_muslo_derecho(i,2),...
        centro_masa_muslo_derecho(i,3),...
        iMuslo_derecho(i,1)*factor_vectores,iMuslo_derecho(i,2)*factor_vectores,iMuslo_derecho(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
    % j del muslo derecho
    quiver3(centro_masa_muslo_derecho(i,1),...
        centro_masa_muslo_derecho(i,2),...
        centro_masa_muslo_derecho(i,3),...
        jMuslo_derecho(i,1)*factor_vectores,jMuslo_derecho(i,2)*factor_vectores,jMuslo_derecho(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
    % k del muslo derecho
    quiver3(centro_masa_muslo_derecho(i,1),...
        centro_masa_muslo_derecho(i,2),...
        centro_masa_muslo_derecho(i,3),...
        kMuslo_derecho(i,1)*factor_vectores,kMuslo_derecho(i,2).*factor_vectores,kMuslo_derecho(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
    
    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    % MUSLO IZQUIERDO
    % i del muslo izquierdo
    quiver3(centro_masa_muslo_izquierdo(i,1),...
        centro_masa_muslo_izquierdo(i,2),...
        centro_masa_muslo_izquierdo(i,3),...
        iMuslo_izquierdo(i,1)*factor_vectores,iMuslo_izquierdo(i,2)*factor_vectores,iMuslo_izquierdo(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
    % j del muslo izquierdo
    quiver3(centro_masa_muslo_izquierdo(i,1),...
        centro_masa_muslo_izquierdo(i,2),...
        centro_masa_muslo_izquierdo(i,3),...
        jMuslo_izquierdo(i,1)*factor_vectores,jMuslo_izquierdo(i,2)*factor_vectores,jMuslo_izquierdo(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
    % k del muslo izquierdo
    quiver3(centro_masa_muslo_izquierdo(i,1),...
        centro_masa_muslo_izquierdo(i,2),...
        centro_masa_muslo_izquierdo(i,3),...
        kMuslo_izquierdo(i,1)*factor_vectores,kMuslo_izquierdo(i,2).*factor_vectores,kMuslo_izquierdo(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    % PIERNA DERECHA
    % i de la pierna derecha
    quiver3(centro_masa_pierna_derecha(i,1),...
        centro_masa_pierna_derecha(i,2),...
        centro_masa_pierna_derecha(i,3),...
        iPierna_derecha(i,1)*factor_vectores,iPierna_derecha(i,2)*factor_vectores,iPierna_derecha(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
    % j de la pierna derecha
    quiver3(centro_masa_pierna_derecha(i,1),...
        centro_masa_pierna_derecha(i,2),...
        centro_masa_pierna_derecha(i,3),...
        jPierna_derecha(i,1)*factor_vectores,jPierna_derecha(i,2)*factor_vectores,jPierna_derecha(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
    % k de la pierna derecha
    quiver3(centro_masa_pierna_derecha(i,1),...
        centro_masa_pierna_derecha(i,2),...
        centro_masa_pierna_derecha(i,3),...
        kPierna_derecha(i,1)*factor_vectores,kPierna_derecha(i,2).*factor_vectores,kPierna_derecha(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    % PIERNA IZQUIERDA
    % i de la pierna izquierda
    quiver3(centro_masa_pierna_izquierda(i,1),...
        centro_masa_pierna_izquierda(i,2),...
        centro_masa_pierna_izquierda(i,3),...
        iPierna_izquierda(i,1)*factor_vectores,iPierna_izquierda(i,2)*factor_vectores,iPierna_izquierda(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
    % j de la pierna izquierda
    quiver3(centro_masa_pierna_izquierda(i,1),...
        centro_masa_pierna_izquierda(i,2),...
        centro_masa_pierna_izquierda(i,3),...
        jPierna_izquierda(i,1)*factor_vectores,jPierna_izquierda(i,2)*factor_vectores,jPierna_izquierda(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
    % k de la pierna izquierda
    quiver3(centro_masa_pierna_izquierda(i,1),...
        centro_masa_pierna_izquierda(i,2),...
        centro_masa_pierna_izquierda(i,3),...
        kPierna_izquierda(i,1)*factor_vectores,kPierna_izquierda(i,2).*factor_vectores,kPierna_izquierda(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    %  PIE DERECHO
    % i del pie derecho
    quiver3(centro_masa_pie_derecho(i,1),...
        centro_masa_pie_derecho(i,2),...
        centro_masa_pie_derecho(i,3),...
        iPie_derecho(i,1)*factor_vectores,iPie_derecho(i,2)*factor_vectores,iPie_derecho(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
    % j del pie derecho
    quiver3(centro_masa_pie_derecho(i,1),...
        centro_masa_pie_derecho(i,2),...
        centro_masa_pie_derecho(i,3),...
        jPie_derecho(i,1)*factor_vectores,jPie_derecho(i,2)*factor_vectores,jPie_derecho(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
    % k del pie derecho
    quiver3(centro_masa_pie_derecho(i,1),...
        centro_masa_pie_derecho(i,2),...
        centro_masa_pie_derecho(i,3),...
        kPie_derecho(i,1)*factor_vectores,kPie_derecho(i,2).*factor_vectores,kPie_derecho(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    % PIE IZQUIERDO
    % i del pie izquierdo
    quiver3(centro_masa_pie_izquierdo(i,1),...
        centro_masa_pie_izquierdo(i,2),...
        centro_masa_pie_izquierdo(i,3),...
        iPie_izquierdo(i,1)*factor_vectores,iPie_izquierdo(i,2)*factor_vectores,iPie_izquierdo(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
    % j del pie izquierdo
    quiver3(centro_masa_pie_izquierdo(i,1),...
        centro_masa_pie_izquierdo(i,2),...
        centro_masa_pie_izquierdo(i,3),...
        jPie_izquierdo(i,1)*factor_vectores,jPie_izquierdo(i,2)*factor_vectores,jPie_izquierdo(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
    % k del pie izquierdo
    quiver3(centro_masa_pie_izquierdo(i,1),...
        centro_masa_pie_izquierdo(i,2),...
        centro_masa_pie_izquierdo(i,3),...
        kPie_izquierdo(i,1)*factor_vectores,kPie_izquierdo(i,2).*factor_vectores,kPie_izquierdo(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
end

grid on
title('Centros de masa junto con los vectores $\vec{i}$, $\vec{j}$ y $\vec{k}$ (PA = posicion articular).','Interpreter','latex')
xlabel('x [m]','Interpreter','latex');
ylabel('y [m]','Interpreter','latex');
zlabel('z [m]','Interpreter','latex');

title('Centros de masa junto con los vectores $\vec{i}$, $\vec{j}$ y $\vec{k}$ (PA = posicion articular).','Interpreter','latex')
legend({'PA cadera derecha',...
    'PA cadera izquierda',...
    'PA rodilla derecha',...
    'PA rodilla izquierda',...
    'PA tobillo derecho',...
    'PA tobillo izquierdo',...
    'Marcador sacro (p15)',...,
    'Muslo derecho ',...,
    'Muslo izquierdo ',...,
    'Pierna derecha',...,
    'Pierna izquierda ',...,
    'Pie derecho ',...,
    'Pie izquierdo ',...,
    '$\vec{i}$',...,
    '$\vec{j}$',...,
    '$\vec{k}$'},'Interpreter','latex')
%%  Plot i,j,k  y centros de masa.
close all % cierra las ventanas
clc % limpia la consola
figure
% Marcadores de la pelvis
plot3(p15(:,1),p15(:,2),p15(:,3),'LineWidth',0.75,'color',negro,'HandleVisibility','off');
hold on
% Marcadores del centro de masa del muslo derecho
plot3(centro_masa_muslo_derecho(:,1),centro_masa_muslo_derecho(:,2),centro_masa_muslo_derecho(:,3),'LineWidth',0.75,'color',negro,'HandleVisibility','off');
hold on
% Marcadores del centro de masa del muslo izquierdo
plot3(centro_masa_muslo_izquierdo(:,1),centro_masa_muslo_izquierdo(:,2),centro_masa_muslo_izquierdo(:,3),'LineWidth',0.75,'color',negro,'HandleVisibility','off');
hold on
% Marcadores del centro de masa de la pierna derecha
plot3(centro_masa_pierna_derecha(:,1),centro_masa_pierna_derecha(:,2),centro_masa_pierna_derecha(:,3),'LineWidth',0.75,'color',negro,'HandleVisibility','off');
hold on
% Marcadores del centro de masa de la pierna izquierda
plot3(centro_masa_pierna_izquierda(:,1),centro_masa_pierna_izquierda(:,2),centro_masa_pierna_izquierda(:,3),'LineWidth',0.75,'color',negro,'HandleVisibility','off');
hold on
% Marcadores del centro de masa del pie derecho
plot3(centro_masa_pie_derecho(:,1),centro_masa_pie_derecho(:,2),centro_masa_pie_derecho(:,3),'LineWidth',0.75,'color',negro,'HandleVisibility','off');
hold on
% Marcadores del centro de masa del pie izquierdo
plot3(centro_masa_pie_izquierdo(:,1),centro_masa_pie_izquierdo(:,2),centro_masa_pie_izquierdo(:,3),'LineWidth',0.75,'color',negro,'HandleVisibility','off');
hold on

% Se plotean los vectores
for i=1:paso_vectores:length(centro_masa_muslo_derecho)
    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    % Se plotean los marcadores donde luego se van a plotear los vectores
    %Marcadores de la pelvis
    plot3(p15(i,1),p15(i,2),p15(i,3),'diamond','LineWidth',line_width,'color',negro);
    hold on
    
    % Marcadores del centro de masa del muslo derecho
    plot3(centro_masa_muslo_derecho(i,1),centro_masa_muslo_derecho(i,2),centro_masa_muslo_derecho(i,3),'o','LineWidth',line_width,'color',negro);
    hold on
    
    % Marcadores del centro de masa del muslo izquierdo
    plot3(centro_masa_muslo_izquierdo(i,1),centro_masa_muslo_izquierdo(i,2),centro_masa_muslo_izquierdo(i,3),'square','LineWidth',line_width,'color',negro);
    hold on
    
    % Marcadores del centro de masa de la pierna derecha
    plot3(centro_masa_pierna_derecha(i,1),centro_masa_pierna_derecha(i,2),centro_masa_pierna_derecha(i,3),'*','LineWidth',line_width,'color',negro);
    hold on
    
    % Marcadores del centro de masa de la pierna izquierda
    plot3(centro_masa_pierna_izquierda(i,1),centro_masa_pierna_izquierda(i,2),centro_masa_pierna_izquierda(i,3),'x','LineWidth',line_width,'color',negro);
    hold on
    
    % Marcadores del centro de masa del pie derecho
    plot3(centro_masa_pie_derecho(i,1),centro_masa_pie_derecho(i,2),centro_masa_pie_derecho(i,3),'pentagram','LineWidth',line_width,'color',negro);
    hold on
    
    % Marcadores del centro de masa del pie izquierdo
    plot3(centro_masa_pie_izquierdo(i,1),centro_masa_pie_izquierdo(i,2),centro_masa_pie_izquierdo(i,3),'>','LineWidth',line_width,'color',negro);
    hold on
    
    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    % SE PLOTEAN LOS VECTORES i j k
    % i,j,k de la pelvis
    quiver3(p15(i,1),p15(i,2),p15(i,3),iPelvis(i,1)*factor_vectores,iPelvis(i,2)*factor_vectores,iPelvis(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
    quiver3(p15(i,1),p15(i,2),p15(i,3),jPelvis(i,1)*factor_vectores,jPelvis(i,2)*factor_vectores,jPelvis(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
    quiver3(p15(i,1),p15(i,2),p15(i,3),kPelvis(i,1)*factor_vectores,kPelvis(i,2).*factor_vectores,kPelvis(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
    
    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    % MUSLO DERECHO
    % i del muslo derecho
    quiver3(centro_masa_muslo_derecho(i,1),...
        centro_masa_muslo_derecho(i,2),...
        centro_masa_muslo_derecho(i,3),...
        iMuslo_derecho(i,1)*factor_vectores,iMuslo_derecho(i,2)*factor_vectores,iMuslo_derecho(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
    % j del muslo derecho
    quiver3(centro_masa_muslo_derecho(i,1),...
        centro_masa_muslo_derecho(i,2),...
        centro_masa_muslo_derecho(i,3),...
        jMuslo_derecho(i,1)*factor_vectores,jMuslo_derecho(i,2)*factor_vectores,jMuslo_derecho(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
    % k del muslo derecho
    quiver3(centro_masa_muslo_derecho(i,1),...
        centro_masa_muslo_derecho(i,2),...
        centro_masa_muslo_derecho(i,3),...
        kMuslo_derecho(i,1)*factor_vectores,kMuslo_derecho(i,2).*factor_vectores,kMuslo_derecho(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
    
    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    % MUSLO IZQUIERDO
    % i del muslo izquierdo
    quiver3(centro_masa_muslo_izquierdo(i,1),...
        centro_masa_muslo_izquierdo(i,2),...
        centro_masa_muslo_izquierdo(i,3),...
        iMuslo_izquierdo(i,1)*factor_vectores,iMuslo_izquierdo(i,2)*factor_vectores,iMuslo_izquierdo(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
    % j del muslo izquierdo
    quiver3(centro_masa_muslo_izquierdo(i,1),...
        centro_masa_muslo_izquierdo(i,2),...
        centro_masa_muslo_izquierdo(i,3),...
        jMuslo_izquierdo(i,1)*factor_vectores,jMuslo_izquierdo(i,2)*factor_vectores,jMuslo_izquierdo(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
    % k del muslo izquierdo
    quiver3(centro_masa_muslo_izquierdo(i,1),...
        centro_masa_muslo_izquierdo(i,2),...
        centro_masa_muslo_izquierdo(i,3),...
        kMuslo_izquierdo(i,1)*factor_vectores,kMuslo_izquierdo(i,2).*factor_vectores,kMuslo_izquierdo(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    % PIERNA DERECHA
    % i de la pierna derecha
    quiver3(centro_masa_pierna_derecha(i,1),...
        centro_masa_pierna_derecha(i,2),...
        centro_masa_pierna_derecha(i,3),...
        iPierna_derecha(i,1)*factor_vectores,iPierna_derecha(i,2)*factor_vectores,iPierna_derecha(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
    % j de la pierna derecha
    quiver3(centro_masa_pierna_derecha(i,1),...
        centro_masa_pierna_derecha(i,2),...
        centro_masa_pierna_derecha(i,3),...
        jPierna_derecha(i,1)*factor_vectores,jPierna_derecha(i,2)*factor_vectores,jPierna_derecha(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
    % k de la pierna derecha
    quiver3(centro_masa_pierna_derecha(i,1),...
        centro_masa_pierna_derecha(i,2),...
        centro_masa_pierna_derecha(i,3),...
        kPierna_derecha(i,1)*factor_vectores,kPierna_derecha(i,2).*factor_vectores,kPierna_derecha(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    % PIERNA IZQUIERDA
    % i de la pierna izquierda
    quiver3(centro_masa_pierna_izquierda(i,1),...
        centro_masa_pierna_izquierda(i,2),...
        centro_masa_pierna_izquierda(i,3),...
        iPierna_izquierda(i,1)*factor_vectores,iPierna_izquierda(i,2)*factor_vectores,iPierna_izquierda(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
    % j de la pierna izquierda
    quiver3(centro_masa_pierna_izquierda(i,1),...
        centro_masa_pierna_izquierda(i,2),...
        centro_masa_pierna_izquierda(i,3),...
        jPierna_izquierda(i,1)*factor_vectores,jPierna_izquierda(i,2)*factor_vectores,jPierna_izquierda(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
    % k de la pierna izquierda
    quiver3(centro_masa_pierna_izquierda(i,1),...
        centro_masa_pierna_izquierda(i,2),...
        centro_masa_pierna_izquierda(i,3),...
        kPierna_izquierda(i,1)*factor_vectores,kPierna_izquierda(i,2).*factor_vectores,kPierna_izquierda(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    %  PIE DERECHO
    % i del pie derecho
    quiver3(centro_masa_pie_derecho(i,1),...
        centro_masa_pie_derecho(i,2),...
        centro_masa_pie_derecho(i,3),...
        iPie_derecho(i,1)*factor_vectores,iPie_derecho(i,2)*factor_vectores,iPie_derecho(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
    % j del pie derecho
    quiver3(centro_masa_pie_derecho(i,1),...
        centro_masa_pie_derecho(i,2),...
        centro_masa_pie_derecho(i,3),...
        jPie_derecho(i,1)*factor_vectores,jPie_derecho(i,2)*factor_vectores,jPie_derecho(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
    % k del pie derecho
    quiver3(centro_masa_pie_derecho(i,1),...
        centro_masa_pie_derecho(i,2),...
        centro_masa_pie_derecho(i,3),...
        kPie_derecho(i,1)*factor_vectores,kPie_derecho(i,2).*factor_vectores,kPie_derecho(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
    %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
    % PIE IZQUIERDO
    % i del pie izquierdo
    quiver3(centro_masa_pie_izquierdo(i,1),...
        centro_masa_pie_izquierdo(i,2),...
        centro_masa_pie_izquierdo(i,3),...
        iPie_izquierdo(i,1)*factor_vectores,iPie_izquierdo(i,2)*factor_vectores,iPie_izquierdo(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
    % j del pie izquierdo
    quiver3(centro_masa_pie_izquierdo(i,1),...
        centro_masa_pie_izquierdo(i,2),...
        centro_masa_pie_izquierdo(i,3),...
        jPie_izquierdo(i,1)*factor_vectores,jPie_izquierdo(i,2)*factor_vectores,jPie_izquierdo(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
    % k del pie izquierdo
    quiver3(centro_masa_pie_izquierdo(i,1),...
        centro_masa_pie_izquierdo(i,2),...
        centro_masa_pie_izquierdo(i,3),...
        kPie_izquierdo(i,1)*factor_vectores,kPie_izquierdo(i,2).*factor_vectores,kPie_izquierdo(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
end

grid on
title('Centros de masa junto con los vectores $\vec{i}$, $\vec{j}$ y $\vec{k}$ (PA = posicion articular).','Interpreter','latex')
xlabel('x [m]','Interpreter','latex');
ylabel('y [m]','Interpreter','latex');
zlabel('z [m]','Interpreter','latex');
title('Centros de masa (CM) junto con los vectores $\vec{i}$, $\vec{j}$ y $\vec{k}$.','Interpreter','latex')
legend({'CM pelvis (p15)',...,
    'CM muslo derecho ',...,
    'CM muslo izquierdo ',...,
    'CM pierna derecha',...,
    'CM pierna izquierda ',...,
    'CM pie derecho ',...,
    'CM pie izquierdo ',...,
    '$\vec{i}$',...,
    '$\vec{j}$',...,
    '$\vec{k}$'},'Interpreter','latex')
%=======================================================================================================================================================%
%%  Plot de los centros de masa y posiciones articulares.
close all % cierra las ventanas
clc % limpia la consola
% Plot centros de masa
figure
plot3(centro_masa_muslo_derecho(:,1),...
    centro_masa_muslo_derecho(:,2),...
    centro_masa_muslo_derecho(:,3),'LineWidth',line_width);
hold on
plot3(centro_masa_muslo_izquierdo(:,1),...
    centro_masa_muslo_izquierdo(:,2),...
    centro_masa_muslo_izquierdo(:,3),'LineWidth',line_width);
hold on
plot3(centro_masa_pierna_derecha(:,1),...
    centro_masa_pierna_derecha(:,2),...
    centro_masa_pierna_derecha(:,3),'LineWidth',line_width);
hold on
plot3(centro_masa_pierna_izquierda(:,1),...
    centro_masa_pierna_izquierda(:,2),...
    centro_masa_pierna_izquierda(:,3),'LineWidth',line_width);
hold on
plot3(centro_masa_pie_derecho(:,1),...
    centro_masa_pie_derecho(:,2),...
    centro_masa_pie_derecho(:,3),'LineWidth',line_width);
hold on
plot3(centro_masa_pie_izquierdo(:,1),...
    centro_masa_pie_izquierdo(:,2),...
    centro_masa_pie_izquierdo(:,3),'LineWidth',line_width);
hold on
%Plot posiciones articulares
plot3(p_cadera_derecha(:,1),p_cadera_derecha(:,2),p_cadera_derecha(:,3),'--','LineWidth',line_width);
hold on
plot3(p_cadera_izquierda(:,1),p_cadera_izquierda(:,2),p_cadera_izquierda(:,3),'--','LineWidth',line_width);
hold on
plot3(p_rodilla_derecha(:,1),p_rodilla_derecha(:,2),p_rodilla_derecha(:,3),'--','LineWidth',line_width);
hold on
plot3(p_rodilla_izquierda(:,1),p_rodilla_izquierda(:,2),p_rodilla_izquierda(:,3),'--','LineWidth',line_width);
hold on
plot3(p_tobillo_derecho(:,1),p_tobillo_derecho(:,2),p_tobillo_derecho(:,3),'--','LineWidth',line_width);
hold on
plot3(p_tobillo_izquierdo(:,1),p_tobillo_izquierdo(:,2),p_tobillo_izquierdo(:,3),'--','LineWidth',line_width);
title('Centros de masa (CM) y posiciones articulares (PA)','Interpreter','latex')
grid on
xlabel('x [m]','Interpreter','latex');
ylabel('y [m]','Interpreter','latex');
zlabel('z [m]','Interpreter','latex');
legend({'CM muslo derecho',...
    'CM muslo izquierdo',...
    'CM pierna derecha',...
    'CM pierna izquierda',...
    'CM pie derecho',...
    'CM pie izquierdo',...
    'PA cadera derecha',...
    'PA cadera izquierda',...
    'PA rodilla derecha',...
    'PA rodilla izquierda',...
    'PA tobillo derecho',...
    'PA tobillo izquierdo'
    },'Interpreter','latex')
%===========================================================================================================================================================%
%%  Plot i,j,k pelvis.
figure
plot3(p15(:,1),p15(:,2),p15(:,3),'LineWidth',line_width,'color',negro)
hold on
for i=1:paso_vectores:length(centro_masa_muslo_derecho)
    % i,j,k de la pelvis
    quiver3(p15(i,1),p15(i,2),p15(i,3),iPelvis(i,1)*factor_vectores,iPelvis(i,2)*factor_vectores,iPelvis(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
    hold on
    quiver3(p15(i,1),p15(i,2),p15(i,3),jPelvis(i,1)*factor_vectores,jPelvis(i,2)*factor_vectores,jPelvis(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
    hold on
    quiver3(p15(i,1),p15(i,2),p15(i,3),kPelvis(i,1)*factor_vectores,kPelvis(i,2).*factor_vectores,kPelvis(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
end
title('Centro de masa de la pelvis (p15)  y vectores $\vec{i}$,$\vec{j}$,$\vec{k}$','Interpreter','latex')
grid on
xlabel('x [m]','Interpreter','latex');
ylabel('y [m]','Interpreter','latex');
zlabel('z [m]','Interpreter','latex');
legend({'Centro de masa',...
    '$\vec{i}$',...
    '$\vec{j}$',...
    '$\vec{k}$'...
    },'Interpreter','latex')
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%%  Plot i,j,k muslo derecho.
close all % cierra las ventanas
clc % limpia la consola
figure
% Plot de la posicion del centro de masa del muslo
plot3(centro_masa_muslo_derecho(:,1),centro_masa_muslo_derecho(:,2),centro_masa_muslo_derecho(:,3),'LineWidth',line_width,'color',negro)
hold on
for i=1:paso_vectores:length(centro_masa_muslo_derecho)
    % i del muslo derecho
    quiver3(centro_masa_muslo_derecho(i,1),...
        centro_masa_muslo_derecho(i,2),...
        centro_masa_muslo_derecho(i,3),...
        iMuslo_derecho(i,1)*factor_vectores,iMuslo_derecho(i,2)*factor_vectores,iMuslo_derecho(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
    % j del muslo derecho
    hold  on
    quiver3(centro_masa_muslo_derecho(i,1),...
        centro_masa_muslo_derecho(i,2),...
        centro_masa_muslo_derecho(i,3),...
        jMuslo_derecho(i,1)*factor_vectores,jMuslo_derecho(i,2)*factor_vectores,jMuslo_derecho(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
    % k del muslo derecho
    hold on
    quiver3(centro_masa_muslo_derecho(i,1),...
        centro_masa_muslo_derecho(i,2),...
        centro_masa_muslo_derecho(i,3),...
        kMuslo_derecho(i,1)*factor_vectores,kMuslo_derecho(i,2).*factor_vectores,kMuslo_derecho(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
end
title('Centro de masa del muslo derecho y vectores $\vec{i}$,$\vec{j}$,$\vec{k}$','Interpreter','latex')
grid on
xlabel('x [m]','Interpreter','latex');
ylabel('y [m]','Interpreter','latex');
zlabel('z [m]','Interpreter','latex');
legend({'Centro de masa',...
    '$\vec{i}$',...
    '$\vec{j}$',...
    '$\vec{k}$'...
    },'Interpreter','latex')
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%%  Plot i,j,k muslo izquierdo
close all % cierra las ventanas
clc % limpia la consola
figure
% Plot de centro de masa
plot3(centro_masa_muslo_izquierdo(:,1),centro_masa_muslo_izquierdo(:,2),centro_masa_muslo_izquierdo(:,3),'LineWidth',line_width,'color',negro)
hold on
for i=1:paso_vectores:length(centro_masa_muslo_izquierdo)
    % i del muslo izquierdo
    quiver3(centro_masa_muslo_izquierdo(i,1),...
        centro_masa_muslo_izquierdo(i,2),...
        centro_masa_muslo_izquierdo(i,3),...
        iMuslo_izquierdo(i,1)*factor_vectores,iMuslo_izquierdo(i,2)*factor_vectores,iMuslo_izquierdo(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
    % j del muslo izquierdo
    hold  on
    quiver3(centro_masa_muslo_izquierdo(i,1),...
        centro_masa_muslo_izquierdo(i,2),...
        centro_masa_muslo_izquierdo(i,3),...
        jMuslo_izquierdo(i,1)*factor_vectores,jMuslo_izquierdo(i,2)*factor_vectores,jMuslo_izquierdo(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
    % k del muslo izquierdo
    hold on
    quiver3(centro_masa_muslo_izquierdo(i,1),...
        centro_masa_muslo_izquierdo(i,2),...
        centro_masa_muslo_izquierdo(i,3),...
        kMuslo_izquierdo(i,1)*factor_vectores,kMuslo_izquierdo(i,2).*factor_vectores,kMuslo_izquierdo(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
end
title('Centro de masa del muslo izquierdo y vectores $\vec{i}$,$\vec{j}$,$\vec{k}$','Interpreter','latex')
grid on
xlabel('x [m]','Interpreter','latex');
ylabel('y [m]','Interpreter','latex');
zlabel('z [m]','Interpreter','latex');
legend({'Centro masa',...
    '$\vec{i}$',...
    '$\vec{j}$',...
    '$\vec{k}$'...
    },'Interpreter','latex')
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%%  Plot i,j,k pierna derecha
close all % cierra las ventanas
clc % limpia la consola
figure
% Plot de centro de masa
plot3(centro_masa_pierna_derecha(:,1),...
    centro_masa_pierna_derecha(:,2),...
    centro_masa_pierna_derecha(:,3),'LineWidth',line_width,'color',negro)
hold on
for i=1:paso_vectores:length(centro_masa_pierna_derecha)
    % i de la pierna derecha
    quiver3(centro_masa_pierna_derecha(i,1),...
        centro_masa_pierna_derecha(i,2),...
        centro_masa_pierna_derecha(i,3),...
        iPierna_derecha(i,1)*factor_vectores,iPierna_derecha(i,2)*factor_vectores,iPierna_derecha(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
    % j de la pierna derecha
    hold  on
    quiver3(centro_masa_pierna_derecha(i,1),...
        centro_masa_pierna_derecha(i,2),...
        centro_masa_pierna_derecha(i,3),...
        jPierna_derecha(i,1)*factor_vectores,jPierna_derecha(i,2)*factor_vectores,jPierna_derecha(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
    % k de la pierna derecha
    hold on
    quiver3(centro_masa_pierna_derecha(i,1),...
        centro_masa_pierna_derecha(i,2),...
        centro_masa_pierna_derecha(i,3),...
        kPierna_derecha(i,1)*factor_vectores,kPierna_derecha(i,2).*factor_vectores,kPierna_derecha(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
end
title('Centro de masa de la pierna derecha y vectores $\vec{i}$,$\vec{j}$,$\vec{k}$','Interpreter','latex')
grid on
xlabel('x [m]','Interpreter','latex');
ylabel('y [m]','Interpreter','latex');
zlabel('z [m]','Interpreter','latex');
legend({'Centro masa',...
    '$\vec{i}$',...
    '$\vec{j}$',...
    '$\vec{k}$'...
    },'Interpreter','latex')
%-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%%  Plot i,j,k pierna izquierda.
close all % cierra las ventanas
clc % limpia la consola
figure
% Plot de centro de masa
plot3(centro_masa_pierna_izquierda(:,1),...
    centro_masa_pierna_izquierda(:,2),...
    centro_masa_pierna_izquierda(:,3),'LineWidth',line_width,'color',negro)
hold on
for i=1:paso_vectores:length(centro_masa_pierna_izquierda)
    % i de la pierna izquierda
    quiver3(centro_masa_pierna_izquierda(i,1),...
        centro_masa_pierna_izquierda(i,2),...
        centro_masa_pierna_izquierda(i,3),...
        iPierna_izquierda(i,1)*factor_vectores,iPierna_izquierda(i,2)*factor_vectores,iPierna_izquierda(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
    % j de la pierna izquierda
    hold  on
    quiver3(centro_masa_pierna_izquierda(i,1),...
        centro_masa_pierna_izquierda(i,2),...
        centro_masa_pierna_izquierda(i,3),...
        jPierna_izquierda(i,1)*factor_vectores,jPierna_izquierda(i,2)*factor_vectores,jPierna_izquierda(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
    % k de la pierna izquierda
    hold on
    quiver3(centro_masa_pierna_izquierda(i,1),...
        centro_masa_pierna_izquierda(i,2),...
        centro_masa_pierna_izquierda(i,3),...
        kPierna_izquierda(i,1)*factor_vectores,kPierna_izquierda(i,2).*factor_vectores,kPierna_izquierda(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
end
title('Centro de masa de la pierna izquierda y vectores $\vec{i}$,$\vec{j}$,$\vec{k}$','Interpreter','latex')
grid on
xlabel('x [m]','Interpreter','latex');
ylabel('y [m]','Interpreter','latex');
zlabel('z [m]','Interpreter','latex');
legend({'Centro masa',...
    '$\vec{i}$',...
    '$\vec{j}$',...
    '$\vec{k}$'...
    },'Interpreter','latex')
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%-
%%  Plot i,j,k pie derecho.
close all % cierra las ventanas
clc % limpia la consola
figure
% Plot de centro de masa
plot3(centro_masa_pie_derecho(:,1),...
    centro_masa_pie_derecho(:,2),...
    centro_masa_pie_derecho(:,3),'LineWidth',line_width,'color',negro)
hold on
for i=1:paso_vectores:length(centro_masa_pie_derecho)
    % i del pie derecho
    quiver3(centro_masa_pie_derecho(i,1),...
        centro_masa_pie_derecho(i,2),...
        centro_masa_pie_derecho(i,3),...
        iPie_derecho(i,1)*factor_vectores,iPie_derecho(i,2)*factor_vectores,iPie_derecho(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
    % j del pie derecho
    hold  on
    quiver3(centro_masa_pie_derecho(i,1),...
        centro_masa_pie_derecho(i,2),...
        centro_masa_pie_derecho(i,3),...
        jPie_derecho(i,1)*factor_vectores,jPie_derecho(i,2)*factor_vectores,jPie_derecho(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
    % k del pie derecho
    hold on
    quiver3(centro_masa_pie_derecho(i,1),...
        centro_masa_pie_derecho(i,2),...
        centro_masa_pie_derecho(i,3),...
        kPie_derecho(i,1)*factor_vectores,kPie_derecho(i,2).*factor_vectores,kPie_derecho(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
end
title('Centro de masa del pie derecho y vectores $\vec{i}$,$\vec{j}$,$\vec{k}$','Interpreter','latex')
grid on
xlabel('x [m]','Interpreter','latex');
ylabel('y [m]','Interpreter','latex');
zlabel('z [m]','Interpreter','latex');
legend({'Centro masa',...
    '$\vec{i}$',...
    '$\vec{j}$',...
    '$\vec{k}$'...
    },'Interpreter','latex')
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%-
%%  Plot i,j,k pie izquierdo.
close all % cierra las ventanas
clc % limpia la consola
figure
% Plot de centro de masa
plot3(centro_masa_pie_izquierdo(:,1),...
    centro_masa_pie_izquierdo(:,2),...
    centro_masa_pie_izquierdo(:,3),'LineWidth',line_width,'color',negro)
hold on
for i=1:paso_vectores:length(centro_masa_pie_izquierdo)
    % i del pie izquierdo
    quiver3(centro_masa_pie_izquierdo(i,1),...
        centro_masa_pie_izquierdo(i,2),...
        centro_masa_pie_izquierdo(i,3),...
        iPie_izquierdo(i,1)*factor_vectores,iPie_izquierdo(i,2)*factor_vectores,iPie_izquierdo(i,3)*factor_vectores,'color',rojo,'LineWidth',line_width)
    % j del pie izquierdo
    hold  on
    quiver3(centro_masa_pie_izquierdo(i,1),...
        centro_masa_pie_izquierdo(i,2),...
        centro_masa_pie_izquierdo(i,3),...
        jPie_izquierdo(i,1)*factor_vectores,jPie_izquierdo(i,2)*factor_vectores,jPie_izquierdo(i,3)*factor_vectores,'color',verde,'LineWidth',line_width)
    % k del pie izquierdo
    hold on
    quiver3(centro_masa_pie_izquierdo(i,1),...
        centro_masa_pie_izquierdo(i,2),...
        centro_masa_pie_izquierdo(i,3),...
        kPie_izquierdo(i,1)*factor_vectores,kPie_izquierdo(i,2).*factor_vectores,kPie_izquierdo(i,3)*factor_vectores,'color',azul,'LineWidth',line_width)
end
title('Centro de masa del pie izquierdo y vectores $\vec{i}$,$\vec{j}$,$\vec{k}$','Interpreter','latex')
grid on
xlabel('x [m]','Interpreter','latex');
ylabel('y [m]','Interpreter','latex');
zlabel('z [m]','Interpreter','latex');
legend({'Centro masa',...
    '$\vec{i}$',...
    '$\vec{j}$',...
    '$\vec{k}$'...
    },'Interpreter','latex')
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%-

%% Renombre de los i,j,k para calculos de los angulos

% Muslos
i1  = iMuslo_derecho;
j1  = jMuslo_derecho;
k1 = kMuslo_derecho;

i2  = iMuslo_izquierdo;
j2  = jMuslo_izquierdo;
k2 = kMuslo_izquierdo;

% Piernas
i3  = iPierna_derecha;
j3  = jPierna_derecha;
k3 = kPierna_derecha;

i4  = iPierna_izquierda;
j4  = jPierna_izquierda;
k4 = kPierna_izquierda;

% Pies
i5  = iPie_derecho;
j5  = jPie_derecho;
k5 = kPie_derecho;

i6  = iPie_izquierdo;
j6  = jPie_izquierdo;
k6 = kPie_izquierdo;
%% Calculos angulos articulares
clc
%-----------------------------------------------  Cadera derecha 
numerador = cross(kPelvis,i1);
denominador = sqrt(sum(numerador.^2,2));
I_cadera_derecha = numerador./denominador;

% alfaR.Hip = sin-1[lR.Hip * iPelvis] 
% alfa: flexion/extension
alfa_cadera_derecha = asin(dot(I_cadera_derecha,iPelvis,2));

% betaR.Hip = sin-1[kPelvis * i1] 
% beta: abduccion/aduccion
beta_cadera_derecha = asin(dot(kPelvis,i1,2));

% gamaR.Hip = -sin-1[lR.Hip * k1] 
% gama : rotacion interna/externa
gama_cadera_derecha = -asin(dot(I_cadera_derecha,k1,2));

%-----------------------------------------------  Cadera izquierda 
numerador = cross(kPelvis,i2);
denominador = sqrt(sum(numerador.^2,2));
I_cadera_izquierda = numerador./denominador;

% alfaL.Hip = sin-1[lL.Hip * iPelvis] 
% alfa: flexion/extension
alfa_cadera_izquierda = asin(dot(I_cadera_izquierda,iPelvis,2));
%alfa_cadera_izquierda=alfa_cadera_izquierda(LHS1-inicio:LHS2-inicio);
% betaL.Hip = -sin-1[kPelvis * i2] 
% beta: abduccion/aduccion
beta_cadera_izquierda = -asin(dot(kPelvis,i2,2));
% gamaL.Hip = sin-1[lL.Hip * k2] 
% gama : rotacion interna/externa
gama_cadera_izquierda = asin(dot(I_cadera_izquierda,k2,2));


%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%-----------------------------------------------  Rodilla derecha 
numerador = cross(k1,i3);
denominador = sqrt(sum(numerador.^2,2));
I_rodilla_derecha = numerador./denominador;

% alfaR.Knee = -sin-1[lR.Knee * i1] 
% alfa: flexion/extension
alfa_rodilla_derecha = -asin(dot(I_rodilla_derecha,i1,2));

% betaR.Knee = sin-1[k1 * i3] 
% beta: abduccion/aduccion
beta_rodilla_derecha = asin(dot(k1,i3,2));

% gamaR.Knee = -sin-1[lR.Knee * k3] 
% gama: rotacion interna/externa
gama_rodilla_derecha = -asin(dot(I_rodilla_derecha,k3,2));

%-----------------------------------------------  Rodilla izquierda 
numerador = cross(k2,i4);
denominador = sqrt(sum(numerador.^2,2));
I_rodilla_izquierda = numerador./denominador;

% alfaL.Knee = -sin-1[lL.Knee * i2] 
% alfa: flexion/extension
alfa_rodilla_izquierda = -asin(dot(I_rodilla_izquierda,i2,2));

% betaL.Knee = -sin-1[k2 * i4]  
% beta: abduccion/aduccion
beta_rodilla_izquierda = -asin(dot(k2,i4,2));

% gamaL.Knee = sin-1[lL.Knee * k4] 
% gama : rotacion interna/externa
gama_rodilla_izquierda  = asin(dot(I_rodilla_izquierda,k4,2));


%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%-----------------------------------------------  Tobillo derecho
numerador = cross(k3,i5);
denominador = sqrt(sum(numerador.^2,2));
I_tobillo_derecho = numerador./denominador;

% alfaR.Ankle = -sin-1[lR.Ankle * j3]  
% alfa: flexion/extension
alfa_tobillo_derecho = -asin(dot(I_tobillo_derecho,j3,2));

% betaR.Ankle = sin-1[k3 * i5] 
% beta: abduccion/aduccion
beta_tobillo_derecho = asin(dot(k3,i5,2));

% gamaR.Ankle = sin-1[lR.Ankle * k5] 
% gama: rotacion interna/externa
gama_tobillo_derecho = asin(dot(I_tobillo_derecho,k5,2));

%-----------------------------------------------  Tobillo izquierdo
numerador = cross(k4,i6);
denominador = sqrt(sum(numerador.^2,2));
I_tobillo_izquierdo = numerador./denominador;

% alfaL.Ankle = -sin-1[lL.Ankle * j4] 
% alfa: flexion/extension
alfa_tobillo_izquierdo= -asin(dot(I_tobillo_izquierdo,j4,2));

% betaL.Ankle = -sin-1[k4 * i6] 
% beta: abduccion/aduccion
beta_tobillo_izquierdo= -asin(dot(k4,i6,2));

% gamaL.Ankle = -sin-1[lL.Ankle * k6] 
% gama: rotacion interna/externa
gama_tobillo_izquierdo= -asin(dot(I_tobillo_izquierdo,k6,2));


%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%

%Recorte angulos articulares

% En este caso la marcha comienza con el pie derecho

%---------- CADERAS
alfa_cadera_derecha      = alfa_cadera_derecha(1:RHS2-inicio);
alfa_cadera_izquierda    = alfa_cadera_izquierda(LHS1-inicio:LHS2-inicio);

beta_cadera_derecha     = beta_cadera_derecha(1:RHS2-inicio);
beta_cadera_izquierda   = beta_cadera_izquierda(LHS1-inicio:LHS2-inicio);

gama_cadera_derecha   = gama_cadera_derecha(1:RHS2-inicio);
gama_cadera_izquierda = gama_cadera_izquierda(LHS1-inicio:LHS2-inicio);

%---------- Rodillas
alfa_rodilla_derecha      = alfa_rodilla_derecha(1:RHS2-inicio);
alfa_rodilla_izquierda    = alfa_rodilla_izquierda(LHS1-inicio:LHS2-inicio);

beta_rodilla_derecha      = beta_rodilla_derecha(1:RHS2-inicio);
beta_rodilla_izquierda    = beta_rodilla_izquierda(LHS1-inicio:LHS2-inicio);

gama_rodilla_derecha   = gama_rodilla_derecha(1:RHS2-inicio);
gama_rodilla_izquierda  = gama_rodilla_izquierda(LHS1-inicio:LHS2-inicio);

%---------- Tobillos
alfa_tobillo_derecho     = alfa_tobillo_derecho(1:RHS2-inicio);
alfa_tobillo_izquierdo    = alfa_tobillo_izquierdo(LHS1-inicio:LHS2-inicio);

beta_tobillo_derecho      = beta_tobillo_derecho(1:RHS2-inicio);
beta_tobillo_izquierdo     = beta_tobillo_izquierdo(LHS1-inicio:LHS2-inicio);

gama_tobillo_derecho    = gama_tobillo_derecho(1:RHS2-inicio);
gama_tobillo_izquierdo   = gama_tobillo_izquierdo(LHS1-inicio:LHS2-inicio);

%% Plot angulos articulares
close all % cierra las ventanas
clc % limpia la consola
figure
subplot(3,3,1);

% Plot de alfa
plot(InterpolaA100Muestras(alfa_cadera_derecha.*factor_a_grados),'LineWidth',line_width)
hold on
plot(InterpolaA100Muestras(alfa_cadera_izquierda.*factor_a_grados),'LineWidth',line_width)
hold on 
plot([60 60 ], [-20 40],'LineWidth',line_width*0.5,'LineStyle','--','color',negro,'HandleVisibility','off');
title('Cadera: flexion (+) $\|$ extesion (-) ','Interpreter','latex','FontSize',12)
grid on
xticks(0:10:100);
xlabel('CM [$\%$]','Interpreter','latex');
ylabel('$\alpha$ [ grados ]','Interpreter','latex','FontSize',13);
legend({'Derecha', 'Izquierda'},'Interpreter','latex','Location','best')
% Plot de beta
subplot(3,3,2);
plot(InterpolaA100Muestras(beta_cadera_derecha.*factor_a_grados),'LineWidth',line_width)
hold on
plot(InterpolaA100Muestras(beta_cadera_izquierda.*factor_a_grados),'LineWidth',line_width)
hold on 
plot([60 60], [-10 10],'LineWidth',line_width*0.5,'LineStyle','--','color',negro,'HandleVisibility','off');
hold on 
plot([60 60 ], [-20 20],'LineWidth',line_width*0.5,'LineStyle','--','color',negro,'HandleVisibility','off');
grid on
xticks(0:10:100);
title('Cadera: abduccion (+) $\|$ aduccion (-) ','Interpreter','latex','FontSize',12)
xlabel('CM [$\%$]','Interpreter','latex');
ylabel('$\beta$ [ grados ]','Interpreter','latex','FontSize',13);
legend({'Derecha', 'Izquierda'},'Interpreter','latex','Location','best')
% Plot de gama
subplot(3,3,3);
plot(InterpolaA100Muestras(gama_cadera_derecha.*factor_a_grados),'LineWidth',line_width)
hold on
plot(InterpolaA100Muestras(gama_cadera_izquierda.*factor_a_grados),'LineWidth',line_width)
hold on 
plot([60 60 ], [-30 20],'LineWidth',line_width*0.5,'LineStyle','--','color',negro,'HandleVisibility','off');
grid on
xticks(0:10:100);
title('Cadera: rot. interna (+) $\|$ rot. externa (-) ','Interpreter','latex','FontSize',12)
xlabel('CM [$\%$]','Interpreter','latex');
ylabel('$\gamma$ [ grados ]','Interpreter','latex','FontSize',13);
legend({'Derecha', 'Izquierda'},'Interpreter','latex','Location','best')

subplot(3,3,4);
% Plot de alfa
plot(InterpolaA100Muestras(alfa_rodilla_derecha.*factor_a_grados),'LineWidth',line_width)
hold on
plot(InterpolaA100Muestras(alfa_rodilla_izquierda.*factor_a_grados),'LineWidth',line_width)
hold on 
plot([60 60 ], [0 60],'LineWidth',line_width*0.5,'LineStyle','--','color',negro,'HandleVisibility','off');
title('Rodilla: flexion (+) $\|$ extesion (-) ','Interpreter','latex','FontSize',12)
grid on
xticks(0:10:100);
xlabel('CM [$\%$]','Interpreter','latex');
ylabel('$\alpha$ [ grados ]','Interpreter','latex','FontSize',13);
legend({'Derecha', 'Izquierda'},'Interpreter','latex','Location','best')
% Plot de beta
subplot(3,3,5);
plot(InterpolaA100Muestras(beta_rodilla_derecha.*factor_a_grados),'LineWidth',line_width)
hold on
plot(InterpolaA100Muestras(beta_rodilla_izquierda.*factor_a_grados),'LineWidth',line_width)
hold on 
plot([60 60 ], [-20 5],'LineWidth',line_width*0.5,'LineStyle','--','color',negro,'HandleVisibility','off');
grid on
xticks(0:10:100);
title('Rodilla: abduccion (+) $\|$ aduccion (-) ','Interpreter','latex','FontSize',12)
xlabel('CM [$\%$]','Interpreter','latex');
ylabel('$\beta$ [ grados ]','Interpreter','latex','FontSize',13);
legend({'Derecha', 'Izquierda'},'Interpreter','latex','Location','best')
% Plot de gama
subplot(3,3,6);
plot(InterpolaA100Muestras(gama_rodilla_derecha.*factor_a_grados),'LineWidth',line_width)
hold on
plot(InterpolaA100Muestras(gama_rodilla_izquierda.*factor_a_grados),'LineWidth',line_width)
hold on 
plot([60 60 ], [-30 20],'LineWidth',line_width*0.5,'LineStyle','--','color',negro,'HandleVisibility','off');
grid on
xticks(0:10:100);
title('Rodilla: rot. interna (+) $\|$ rot. externa (-) ','Interpreter','latex','FontSize',12)
xlabel('CM [$\%$]','Interpreter','latex');
ylabel('$\gamma$ [ grados ]','Interpreter','latex','FontSize',13);
legend({'Derecha', 'Izquierda'},'Interpreter','latex','Location','best')


subplot(3,3,7);
% Plot de alfa
plot(InterpolaA100Muestras(alfa_tobillo_derecho.*factor_a_grados),'LineWidth',line_width)
hold on
plot(InterpolaA100Muestras(alfa_tobillo_izquierdo.*factor_a_grados),'LineWidth',line_width)
hold on 
plot([60 60 ], [-30 20],'LineWidth',line_width*0.5,'LineStyle','--','color',negro,'HandleVisibility','off');
title('Tobillo: flexion (+) $\|$ extesion (-) ','Interpreter','latex','FontSize',12)
grid on
xticks(0:10:100);
xlabel('CM [$\%$]','Interpreter','latex');
ylabel('$\alpha$ [ grados ]','Interpreter','latex','FontSize',13);
legend({'Derecha', 'Izquierda'},'Interpreter','latex','Location','best')
% Plot de beta
subplot(3,3,9);
plot(InterpolaA100Muestras(beta_tobillo_derecho.*factor_a_grados),'LineWidth',line_width)
hold on
plot(InterpolaA100Muestras(beta_tobillo_izquierdo.*factor_a_grados),'LineWidth',line_width)
hold on 
plot([60 60 ], [-30 20],'LineWidth',line_width*0.5,'LineStyle','--','color',negro,'HandleVisibility','off');
grid on
xticks(0:10:100);
title('Tobillo: abduccion (+) $\|$ aduccion (-) ','Interpreter','latex','FontSize',12)
xlabel('CM [$\%$]','Interpreter','latex');
ylabel('$\beta$ [ grados ]','Interpreter','latex','FontSize',13);
legend({'Derecha', 'Izquierda'},'Interpreter','latex','Location','best')
% Plot de gama
subplot(3,3,8);
plot(InterpolaA100Muestras(gama_tobillo_derecho.*factor_a_grados),'LineWidth',line_width)
hold on
plot(InterpolaA100Muestras(gama_tobillo_izquierdo.*factor_a_grados),'LineWidth',line_width)
hold on 
plot([60 60 ], [15 40],'LineWidth',line_width*0.5,'LineStyle','--','color',negro,'HandleVisibility','off');
grid on
xticks(0:10:100);
title('Tobillo: inversion (+) $\|$ eversion (-) ','Interpreter','latex','FontSize',12)
xlabel('CM [$\%$]','Interpreter','latex');
ylabel('$\gamma$ [ grados ]','Interpreter','latex','FontSize',13);
legend({'Derecha', 'Izquierda'},'Interpreter','latex','Location','best')

%% Calculos de los angulos de Euler

plotear_angulos_euler = true;
   % Muslo derecho
 [  alfa_s_md, alfa_c_md,...
     beta_s_md, beta_c_md,...
     gamma_s_md, gamma_c_md] = Biomecanica.AngulosEuler(i1,j1,k1,plotear_angulos_euler);
 
    % Muslo izquierdo
 [  alfa_s_mi, alfa_c_mi,...
     beta_s_mi, beta_c_mi,...
     gamma_s_mi, gamma_c_mi] = Biomecanica.AngulosEuler(i2,j2,k2,plotear_angulos_euler);
 
    % Pierna derecha
 [  alfa_s_pd, alfa_c_pd,...
     beta_s_pd, beta_c_pd,...
     gamma_s_pd, gamma_c_pd] = Biomecanica.AngulosEuler(i3,j3,k3,plotear_angulos_euler);
 
    % Pierna izquierda
 [  alfa_s_pi, alfa_c_pi,...
     beta_s_pi, beta_c_pi,...
     gamma_s_pi, gamma_c_pi] = Biomecanica.AngulosEuler(i4,j4,k4,plotear_angulos_euler);
  
    % Pie derecho
 [  alfa_s_pid, alfa_c_pid,...
     beta_s_pid, beta_c_pid,...
     gamma_s_pid, gamma_c_pid] = Biomecanica.AngulosEuler(i5,j5,k5,plotear_angulos_euler);
 
     % Pie izquierda
 [  alfa_s_pii, alfa_c_pii,...
     beta_s_pii, beta_c_pii,...
     gamma_s_pii, gamma_c_pii] = Biomecanica.AngulosEuler(i6,j6,k6,plotear_angulos_euler);
 %% Calculos de las velocidades angulares
    % Muslo derecho
 [  derivada_alfa_md,...
     derivada_beta_md,...
     derivada_gamma_md] = Biomecanica.Derivadas(alfa_s_md,alfa_c_md,alfa_s_md,fm);
 
 % Muslo izquierdo
  [  derivada_alfa_mi,...
     derivada_beta_mi,...
     derivada_gamma_mi] = Biomecanica.Derivadas(alfa_s_mi,alfa_c_mi,alfa_s_mi,fm);
 
  % Pierna derecho
  [  derivada_alfa_pd,...
     derivada_beta_pd,...
     derivada_gamma_pd] = Biomecanica.Derivadas(alfa_s_pd,alfa_c_pd,alfa_s_pd,fm);
  
  % Pierna izquierda
  [  derivada_alfa_pi,...
     derivada_beta_pi,...
     derivada_gamma_pi] = Biomecanica.Derivadas(alfa_s_pi,alfa_c_pi,alfa_c_pi,fm);
  
  % Pie derecho
  [  derivada_alfa_pid,...
     derivada_beta_pid,...
     derivada_gamma_pid] = Biomecanica.Derivadas(alfa_s_pid,alfa_c_pid,alfa_s_pid,fm);
 
   % Pie izquierdo
  [  derivada_alfa_pii,...
     derivada_beta_pii,...
     derivada_gamma_pii] = Biomecanica.Derivadas(alfa_s_pii,alfa_c_pii,alfa_s_pii,fm);
 
 
    