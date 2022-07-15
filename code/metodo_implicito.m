%r = (m-1/2)*dr
%z = (n-1)*dz
%p = (p-1)*dt
%M = numero de nodos radiales;
%N = numero de nodos en el eje z;
%P = numero de nodos temporales;
%tc = tiempo total;
% dt = tc/(P-1);
% dr = 2*R/(2*M-1);
% dz = L/(N-1);
%F_o = constante de Fourier
%Bi = numero de Biot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculo numero de Fourier y numero de Biot:
% alpha=k/p*c
% Bi=(h*dr)/k
% F_o=alpha*dt/(dr^2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%T_o = Temperatura inicial del cuerpo
%T_amb = Temperatura del ambiente
%R = radio del cilindro
%L = longitud del cilindro

function [T]=metodo_implicito(F_o, Bi, T_o, T_amb, R, L, M, P)

N = 1+(L/R)*(M-0.5);
%T: Matriz de temperatura para el analisis transitorio.
T = zeros(M,N,P);
%T_v: Matriz simbolica usada para guardar los valores de las ecuaciones
%resueltas
T_v = sym('T_%d_%d', [M N]);
%Estableciendo la temperatura inicial del cilindro
T(:,:,1)=ones(M,N)*T_o;
%Inicializando una matriz para guardar las ecuaciones
EQN = sym('eqn',[M N]);
%vector auxiliar para ordenamiento de variables antes de
%resolver sistema de ecuaciones
VAR = sym('var', [1 M*N]);

for p=1:P-1
    for m=1:M
        m_a = m-1;
        for n=1:N
            %punto interior
            if (m<M && m>1 && n<N && n>1)
                EQN(m,n) = T_v(m,n)*(1+4*F_o)-F_o*(T_v(m,n-1)...
                +T_v(m,n+1)+((2*m_a)/(2*m_a+1))*T_v(m-1,n)+...
                ((2*m_a+2)/(2*m_a+1))*T_v(m+1,n)) == T(m,n,p); 
            end

            %punto exterior cara lateral
            if (m==M && n<N && n>1)
                EQN(m,n) = T_v(m,n)*(1+2*F_o+F_o*((8*m_a)/(4*m_a+1))+...
                F_o*Bi*((8*m_a+4)/(4*m_a+1)))-F_o*(((8*m_a)/(4*m_a+1))*T_v(m-1,n)...
                +T_v(m,n-1)+T_v(m,n+1)+Bi*((8*m_a+4)/(4*m_a+1))*T_amb)...
                == T(m,n,p);
            end

            %punto exterior base superior
            if (m>1 && m<M && n==N)
                EQN(m,n) = T_v(m,n)*(1+4*F_o+2*F_o*Bi)-2*F_o*...
                (T_v(m,n-1)+Bi*T_amb+((m_a+1)/(2*m_a+1))*T_v(m+1,n)...
                +(m_a/(2*m_a+1))*T_v(m-1,n))==T(m,n,p);
            end

            %punto exterior base inferior
            if (m>1 && m<M && n==1)
                EQN(m,n) = T_v(m,n)*(1+4*F_o+2*F_o*Bi)-2*F_o*...
                (T_v(m,n+1)+Bi*T_amb+((m_a+1)/(2*m_a+1))*T_v(m+1,n)...
                +(m_a/(2*m_a+1))*T_v(m-1,n))==T(m,n,p);
            end

            %punto exterior arista superior
            if (m==M && n==N) 
                EQN(m,n) = T_v(m,n)*(1+2*((8*m_a+1)/(4*m_a+1))*F_o+...
                2*((8*m_a+3)/(4*m_a+1))*F_o*Bi)-2*F_o*(T_v(m,n-1)+(4*m_a/(4*m_a+1))*...
                T_v(m-1,n)+Bi*((8*m_a+3)/(4*m_a+1))*T_amb)==T(m,n,p);
            end

            %punto exterior arista inferior
            if (m==M && n==1) 
                EQN(m,n) = T_v(m,n)*(1+2*((8*m_a+1)/(4*m_a+1))*F_o+...
                2*((8*m_a+3)/(4*m_a+1))*F_o*Bi)-2*F_o*(T_v(m,n+1)+(4*m_a/(4*m_a+1))*...
                T_v(m-1,n)+Bi*((8*m_a+3)/(4*m_a+1))*T_amb)==T(m,n,p);
            end
            
            %ecuaciones centro de cilindro
            %punto interior
            if (m==1 && n<N && n>1)
                EQN(m,n) = T_v(m,n)*(1+2*F_o+F_o*((2*m_a+2)/(2*m_a+1)))-F_o*(T_v(m,n-1)...
                +T_v(m,n+1)+ ((2*m_a+2)/(2*m_a+1))*T_v(m+1,n)) == T(m,n,p); 
            end
            %punto exterior base superior
            if (m==1 && n==N)
                EQN(m,n) = T_v(m,n)*(1+2*F_o+2*F_o*Bi+F_o*((2*m_a+2)/(2*m_a+1)))-2*F_o*...
                (T_v(m,n-1)+Bi*T_amb+((m_a+1)/(2*m_a+1))*T_v(m+1,n))==T(m,n,p);
            end

            %punto exterior base inferior
            if (m==1 && n==1)
                EQN(m,n) = T_v(m,n)*(1+2*F_o+2*F_o*Bi+F_o*((2*m_a+2)/(2*m_a+1)))-2*F_o*...
                (T_v(m,n+1)+Bi*T_amb+((m_a+1)/(2*m_a+1))*T_v(m+1,n))==T(m,n,p);
            end
            
        end
    end

    for m=1:M
        VAR(((m-1)*N+1):(N*m))= sort(symvar(T_v(m,:)));
    end
    S = solve(EQN,VAR);
    S1 = struct2cell(S);
    
    for m=1:M
        for n=1:N
            T(m,n,p+1) = S1{N*(m-1)+n};
        end
    end
end