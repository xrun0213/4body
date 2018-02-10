# 4body
clc;clear

format long

k = 0.31624349300710735377439789256976;
r = 5/4; w = 1.722246782664867009854618407143*sqrt(1/(r^3));
h = 1e-4;

x0 = [r; 0; 0; w*r;
     k*r; 0; 0; w*k*r;
     -k*r; 0; 0; -w*k*r;
     -r; 0; 0; -w*r];

times = 1100;
 
x_euler = x0;
x_rk2 = x0;
x_rk4 = x0;
x_gps2 = x0;
x_midp = x0;

for i = 1:times
        
    %exact solution
    t(i) = i*h;
    tm1(i, :) = [r*cos(w*t(i)) r*sin(w*t(i))]; 
    tm2(i, :) = [k*r*cos(w*t(i)) k*r*sin(w*t(i))];
    tm3(i, :) = [-k*r*cos(w*t(i)) -k*r*sin(w*t(i))]; 
    tm4(i, :) = [-r*cos(w*t(i)) -r*sin(w*t(i))];
    
    %Euler
    X_Eu(:, i) = x_euler;
    
    f = dx(x_euler);
    x_euler = eu(f, h, x_euler);
    
    m1_Eu(i, :) = [x_euler(1) x_euler(3)];
    m2_Eu(i, :) = [x_euler(5) x_euler(7)];
    m3_Eu(i, :) = [x_euler(9) x_euler(11)];
    m4_Eu(i, :) = [x_euler(13) x_euler(15)];
    
    E_Eu(i,:) = [norm([m1_Eu(i,:)-tm1(i,:)]), norm([m2_Eu(i,:)-tm2(i,:)]), norm([m3_Eu(i,:)-tm3(i,:)]), norm([m4_Eu(i,:)-tm4(i,:)])];
    
    %RK2
    X_RK2(:, i) = x_rk2;
        
    f = dx(x_rk2);
    x_rk2 = rk2(f, h, x_rk2);
    
    m1_RK2(i, :) = [x_rk2(1) x_rk2(3)];
    m2_RK2(i, :) = [x_rk2(5) x_rk2(7)];
    m3_RK2(i, :) = [x_rk2(9) x_rk2(11)];
    m4_RK2(i, :) = [x_rk2(13) x_rk2(15)];
    
    E_RK2(i,:) = [norm([m1_RK2(i,:)-tm1(i,:)]), norm([m2_RK2(i,:)-tm2(i,:)]), norm([m3_RK2(i,:)-tm3(i,:)]), norm([m4_RK2(i,:)-tm4(i,:)])];
    
    %RK4
    X_RK4(:, i) = x_rk4;
    
    f = dx(x_rk4);
    x_rk4 = rk4(f, h, x_rk4);
    
    m1_RK4(i, :) = [x_rk4(1) x_rk4(3)];
    m2_RK4(i, :) = [x_rk4(5) x_rk4(7)];
    m3_RK4(i, :) = [x_rk4(9) x_rk4(11)];
    m4_RK4(i, :) = [x_rk4(13) x_rk4(15)];
    
    E_RK4(i,:) = [norm([m1_RK4(i,:)-tm1(i,:)]), norm([m2_RK4(i,:)-tm2(i,:)]), norm([m3_RK4(i,:)-tm3(i,:)]), norm([m4_RK4(i,:)-tm4(i,:)])];
    
    %GPS2
    X_GPS2(:, i) = x_gps2;
    
    f = dx(x_gps2);
    x_gps2 = gps2(f, h, x_gps2);
    
    m1_GPS2(i, :) = [x_gps2(1) x_gps2(3)];
    m2_GPS2(i, :) = [x_gps2(5) x_gps2(7)];
    m3_GPS2(i, :) = [x_gps2(9) x_gps2(11)];
    m4_GPS2(i, :) = [x_gps2(13) x_gps2(15)];
    
    E_GPS2(i,:) = [norm([m1_GPS2(i,:)-tm1(i,:)]), norm([m2_GPS2(i,:)-tm2(i,:)]), norm([m3_GPS2(i,:)-tm3(i,:)]), norm([m4_GPS2(i,:)-tm4(i,:)])];
    
    %MidPoint
    X_MidP(:, i) = x_midp;
    
    f = dx(x_midp);
    x_midp = midp(f, h, x_midp);
    
    m1_MidP(i, :) = [x_midp(1) x_midp(3)];
    m2_MidP(i, :) = [x_midp(5) x_midp(7)];
    m3_MidP(i, :) = [x_midp(9) x_midp(11)];
    m4_MidP(i, :) = [x_midp(13) x_midp(15)];
    
    E_MidP(i,:) = [norm([m1_MidP(i,:)-tm1(i,:)]), norm([m2_MidP(i,:)-tm2(i,:)]), norm([m3_MidP(i,:)-tm3(i,:)]), norm([m4_MidP(i,:)-tm4(i,:)])];
    
end
 
hold on
grid on

%error
plot(log10(1:i), log10(E_Eu(:, 1)), 'b.', 'Markersize', 10)
plot(log10(1:i), log10(E_RK2(:, 1)), 'g.', 'Markersize', 10)
plot(log10(1:i), log10(E_RK4(:, 1)), 'black.', 'Markersize', 10)
plot(log10(1:i), log10(E_GPS2(:, 1)), 'r.', 'Markersize', 10)
plot(log10(1:i), log10(E_MidP(:, 1)), 'm.', 'Markersize', 10)

%orbit
%plot(X_Eu(1,1:i), X_Eu(3,1:i), 'b.', 'Markersize', 10)

legend('Euler', 'RK2', 'RK4', 'GPS2', 'MidPoint', 2)
%axis square

title('collinear configuration errors')
