function graf_dist_temp(R, L, M, T, p_e, T_amb, T_o)
%p_e = instante de tiempo
%N = numero de nodos en el eje z
%M = numero de nodos radiales
%dr = delta de r
%dz = delta de z
%r = parametro radial del cilindro
%z = parametro z del cilindro
%T = matriz de distribucion de temperatura
%T_amb = Temperatura del ambiente
%T_o = Temperatura inicial del cuerpo

N = 1+(L/R)*(M-0.5);
dr = 2*R/(2*M-1);
r = zeros(1,M);
dz = L/(N-1);
z = zeros(1,N);
for m = 1:M
    r(m) = (m-1/2)*dr;
end
for n = 1:N
    z(n) = (n-1)*dz;
end
T_a = T(:,:,p_e);
imshow(T_a, [T_amb T_o], 'XData', z, 'Ydata', r, 'InitialMagnification', 1600);
colormap(gca, jet(256));
colorbar(gca);
axis('on','image');
title('T(r,z)', 'Fontsize', 15);
xlabel('z', 'Fontsize', 15);
ylabel('r', 'Fontsize', 15);
