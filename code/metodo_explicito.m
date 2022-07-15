function [T]=metodo_explicito(F_o, Bi, T_o, T_amb, R, L, M, P)

N = 1+(L/R)*(M-0.5);
%T: Matriz de temperatura para el analisis transitorio.
T = zeros(M,N,P);
%Estableciendo la temperatura inicial del cilindro
T(:,:,1)=ones(M,N)*T_o;


for p=1:P-1
    for m=1:M
        m_a = m-1;
        for n=1:N
            %punto interior
            if (m<M && m>1 && n<N && n>1)
                T(m,n,p+1)= T(m,n,p)*(1-4*F_o)+F_o*(T(m,n-1,p)...
                +T(m,n+1,p)+((2*m_a)/(2*m_a+1))*T(m-1,n,p)+...
                ((2*m_a+2)/(2*m_a+1))*T(m+1,n,p)); 
            end

            %punto exterior cara lateral
            if (m==M && n<N && n>1)
                T(m,n,p+1) = T(m,n,p)*(1-2*F_o-F_o*((8*m_a)/(4*m_a+1))-...
                F_o*Bi*((8*m_a+4)/(4*m_a+1)))+F_o*(((8*m_a)/(4*m_a+1))*T(m-1,n,p)...
                +T(m,n-1,p)+T(m,n+1,p)+Bi*((8*m_a+4)/(4*m_a+1))*T_amb);
            end

            %punto exterior base superior
            if (m>1 && m<M && n==N)
                T(m,n,p+1) = T(m,n,p)*(1-4*F_o-2*F_o*Bi)+2*F_o*...
                (T(m,n-1,p)+Bi*T_amb+((m_a+1)/(2*m_a+1))*T(m+1,n,p)...
                +(m_a/(2*m_a+1))*T(m-1,n,p));
            end

            %punto exterior base inferior
            if (m>1 && m<M && n==1)
                T(m,n,p+1) = T(m,n,p)*(1-4*F_o-2*F_o*Bi)+2*F_o*...
                (T(m,n+1,p)+Bi*T_amb+((m_a+1)/(2*m_a+1))*T(m+1,n,p)...
                +(m_a/(2*m_a+1))*T(m-1,n,p));
            end

            %punto exterior arista superior
            if (m==M && n==N) 
                T(m,n,p+1) = T(m,n,p)*(1-2*((8*m_a+1)/(4*m_a+1))*F_o-...
                2*((8*m_a+3)/(4*m_a+1))*F_o*Bi)+2*F_o*(T(m,n-1,p)+(4*m_a/(4*m_a+1))*...
                T(m-1,n,p)+Bi*((8*m_a+3)/(4*m_a+1))*T_amb);
            end

            %punto exterior arista inferior
            if (m==M && n==1) 
                T(m,n,p+1) = T(m,n,p)*(1-2*((8*m_a+1)/(4*m_a+1))*F_o-...
                2*((8*m_a+3)/(4*m_a+1))*F_o*Bi)+2*F_o*(T(m,n+1,p)+(4*m_a/(4*m_a+1))*...
                T(m-1,n,p)+Bi*((8*m_a+3)/(4*m_a+1))*T_amb);
            end
            
            %ecuaciones centro de cilindro
            %punto interior
            if (m==1 && n<N && n>1)
                T(m,n,p+1) = T(m,n,p)*(1-2*F_o-F_o*((2*m_a+2)/(2*m_a+1)))+F_o*(T(m,n-1,p)...
                +T(m,n+1,p)+ ((2*m_a+2)/(2*m_a+1))*T(m+1,n,p)); 
            end
            %punto exterior base superior
            if (m==1 && n==N)
                T(m,n,p+1) = T(m,n,p)*(1-2*F_o-2*F_o*Bi-F_o*((2*m_a+2)/(2*m_a+1)))+2*F_o*...
                (T(m,n-1,p)+Bi*T_amb+((m_a+1)/(2*m_a+1))*T(m+1,n,p));
            end

            %punto exterior base inferior
            if (m==1 && n==1)
                T(m,n,p+1) = T(m,n,p)*(1-2*F_o-2*F_o*Bi-F_o*((2*m_a+2)/(2*m_a+1)))+2*F_o*...
                (T(m,n+1,p)+Bi*T_amb+((m_a+1)/(2*m_a+1))*T(m+1,n,p));
            end
            
        end
    end
end