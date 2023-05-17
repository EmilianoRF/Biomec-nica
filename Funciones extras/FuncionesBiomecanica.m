classdef FuncionesBiomecanica
   methods (Static)
      function [alfa_s, alfa_c, beta_s, beta_c, gamma_s, gamma_c] = AngulosEuler(i,j,k,graficar)
         % Se definen I, J, K
         I           = zeros(size(i));
         I(:,1)    = 1;  
         J           = zeros(size(j));
         J(:,2)    = 1;
         K         = zeros(size(k));
         K(:,3)  = 1;
         %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
         % Calculo de la linea de nodos
         numerador      = cross(K,k);
         denominador  = sqrt(sum(numerador.^2,2));
         nodo                  = numerador./denominador;
        %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
         % Calculo de los angulos con seno
         alfa_s         = asind(dot(cross(I,nodo),K,2));
         beta_s        = asind(dot(cross(K,k),nodo,2));
         gamma_s  = asind(dot(cross(nodo,i),k,2));
        %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%         
         % Calculo de los angulos con coseno
         numerador      = dot(J,nodo,2);
         denominador  = sqrt(sum(numerador.^2,2));
         factor                = numerador./denominador;
         alfa_c                = -dot(factor,acosd(dot(I,nodo,2)),2);
        %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%         
         beta_c              = acosd(dot(K,k,2));
        %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%        
         numerador      = dot(j,nodo,2);
         denominador  = sqrt(sum(numerador.^2,2));
         factor                = numerador./denominador;
         gamma_c        = -dot(factor,acosd(dot(i,nodo,2)),2);
         
         if graficar
             line_width = 1.75;
             factor_a_grados    = 1;
             figure
             subplot(2,3,1)
             plot(factor_a_grados*alfa_c,'LineWidth',line_width,'color','b')
             title('Coseno - $\alpha$','Interpreter','latex')
             xlabel('Muetras','Interpreter','latex')
             ylabel('Grados','Interpreter','latex')
             grid on
             
             subplot(2,3,2)
             plot(factor_a_grados*beta_c,'LineWidth',line_width,'color','b')
             title('Coseno - $\beta$','Interpreter','latex')
             xlabel('Muetras','Interpreter','latex')
             ylabel('Grados','Interpreter','latex')
             grid on
             
             subplot(2,3,3)
             plot(factor_a_grados*gamma_c,'LineWidth',line_width,'color','b')
             title('Coseno - $\gamma$','Interpreter','latex')
             xlabel('Muetras','Interpreter','latex')
             ylabel('Grados','Interpreter','latex')
             grid on
            
             subplot(2,3,4)
             plot(factor_a_grados*alfa_s,'LineWidth',line_width,'color','r')
             title('Seno - $\alpha$','Interpreter','latex')
             xlabel('Muetras','Interpreter','latex')
             ylabel('Grados','Interpreter','latex')
             grid on
             
            subplot(2,3,5)
            plot(factor_a_grados*beta_s,'LineWidth',line_width,'color','r')
            title('Seno - $\beta$','Interpreter','latex')
             xlabel('Muetras','Interpreter','latex')
             ylabel('Grados','Interpreter','latex')
             grid on
             
            subplot(2,3,6)
            plot(factor_a_grados*gamma_s,'LineWidth',line_width,'color','r')
            title('Seno - $\gamma$','Interpreter','latex')
            xlabel('Muetras','Interpreter','latex')
             ylabel('Grados','Interpreter','latex')
             grid on
         end
      end
      function [deruno,derdos,dertres] = Derivadas(uno,dos,tres,fm)
    % function [deruno,derdos,dertres] = derivadas(uno,dos,tres,fm)
    % Esta funci�n devuelve las derivadas de las se�ales uno, dos y tres, las
    % cuales fueron muestreadas con una frecuencia de muestreo fm. Uno, dos y
    % tres ser�n vectores de una columna y n filas. Las se�ales derivadas
    % deruno, derdos y dertres tambi�n ser�n vectores de nx1.
    % El primer y el segundo valor del vector derivada ser�n iguales.
    % Entradas:
    % uno, dos y tres = cada uno es un vector [n x1] que se quiera derivar.
    % fm = es la frecuencia de muestreo utilizada.

% Salidas
% deruno = vector de tama�o [n x 1], y que representa la derivada del
% vector uno. deruno[1]=deruno[n]=0.
% derdos = vector de tama�o [n x 1], y que representa la derivada del
% vector dos. derdos[1]=derdos[n]=0.
% dertres = vector de tama�o [n x 1], y que representa la derivada del
% vector tres. dertres[1]=dertres[n]=0.
% Paola Catalfamo 15-09-2020


deruno = zeros(size (uno));        % Armo las derivadas con la misma longitud, para no tener inconvenientes con la longitud de las se�ales
derdos = zeros(size (dos));
dertres = zeros(size (tres));

deltat=1/fm;

for n=2: (length (uno)) -1                 % n va desde 2 hasta fin-1
    
    deruno(n) = (uno(n+1) - uno(n-1)) ./ (2*deltat);
    
    derdos(n) = (dos(n+1) - dos(n-1)) ./ (2*deltat);
    
    dertres (n) = (tres(n+1) - tres(n-1)) ./ (2*deltat);
    
end;

deruno(1)= deruno(2);
derdos (1) = derdos (2);
dertres (1) = dertres (2);
      end
      function [wx,wy,wz] = VelocidadAngularSegmento (alfa, deralfa, beta, derbeta, gama,dergama)
% function [wx,wy,wz] = velocidadangularsegmento (alfa, deralfa, beta, derbeta, gama,dergama);
% Esta funci�n calcula las velocidades angulares de los segmentos en
% coordenadas locales, a partir de los �ngulos de Euler alfa, beta, gama y
% sus derivadas expresadas en coordenadas locales.

% Entradas:
% Todas las entradas son vectores de tama�o n x 1, donde n es el n�mero de
% muestras correspondientes a dos ciclos completos de la marcha.
% alfa= angulo alfa del segmento de inter�s
% deralfa = derivada del �ngulo alfa del segmento de inter�s
% beta= angulo beta del segmento de inter�s
% derbeta = derivada del �ngulo beta del segmento de inter�s
% gama= angulo gamma del segmento de inter�s
% dergama = derivada del �ngulo gamma del segmento de inter�s

% Salidas de la funci�n
% Todas las salidas son vectores de tama�o nx1, donde n es el n�mero de
% muestras correspondientes a dos ciclos completos de la marcha.

% wx = velocidad angular del segmento de inter�s en la direcci�n x (es
% decir en la direcci�n i local), wx est� expresado en coordenadas locales

% wy = velocidad angular del segmento de inter�s en la direcci�n y (es
% decir en la direcci�n j local), wy est� expresado en coordenadas locales

% wz = velocidad angular del segmento de inter�s en la direcci�n z (es
% decir en la direcci�n k local), wx est� expresado en coordenadas locales

% Ejemplo de c�mo llamar a la funci�n para el segmento MUSLO DERECHO:
% [wxmusder,wymusder,wzmusder] = velocidadangularsegmento (alfamusder, deralfamusder, betamusder, derbetamusder, gamamusder,dergamamusder)

wx= deralfa.*sind(beta).*sind(gama) + derbeta.*cosd(gama);
wy = deralfa.* sind(beta).*cosd(gama) - derbeta.*sind(gama);
wz = deralfa.*cosd(beta) + dergama;

      end
     function [senialnormalizada] = InterpolaA100Muestras(senial)
% function [senialnormalizada] = InterpolaA100Muestras(senial);
% Esta funci�n recibe una se�al en la variable "senial"
% y la interpola a 100 muestras en la variable senialporciento;
% La se�al que ingresa es un vector nx1
% Y la senialnormalizada es un vector que contiene 1x100 datos
% y que nosotros entendemos como normalizada al ciclo de la marcha.
% Ejemplo de uso:
% [angulocaderanorm] = InterpolaA100Muestras(angulocadera);

 muestras = 1:1:length(senial);
 
 Porcentaje = 1:length(senial)/100:length(senial);
 
senialnormalizada  = interp1(muestras,senial,Porcentaje,'spline');
     end
     function [aceleracion] = Aceleracion(senial,fm)
         aceleracion = zeros(size(senial));
         for n=2:(length(senial)-1)
             aceleracion(n,:)   = (senial(n+1,:)-2*senial(n,:) + senial(n-1,:))./(1/fm)^2;
         end
         aceleracion(1,:) = aceleracion(2,:);
     end
 end
end 