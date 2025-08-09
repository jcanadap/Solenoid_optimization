clear; clc;
close all;
Rin = 2; % Internal radius (mm)
wc = 0.6; % Conductive trace width (mm)
hc = 0.6; % Conductive trace height (mm)

mu0 = 1.26/1000000; % Vacuum permeability (H/m)
ro = 0.00032; % Electrifi resistivity (Ohm*m) cheked with 10-layered solenoid resistance (3 kOhm)
Imax = 35; % Maximum recommended current (mA)
Pmax = 180; % Maximum suppliable power (W)
C = 1; % Core multiplication factor
wi = 0.4; % Insulation width (mm)
hi = 0.05; % Insulation height (mm)
Rout = Rin:1:30; % Outside radius (mm)
L = 0.3:0.025:20; % Inductor length (mm)
z_offset = 0.7; % Distance from coil to probe (mm)

Vmax = Pmax/(Imax/1000);
Bo_out_of_plane = zeros(size(Rout,2),size(L,2));
R_out_of_plane = zeros(size(Rout,2),size(L,2));
V_out_of_plane = zeros(size(Rout,2),size(L,2));

for l=1:size(L,2)
    N_out = round(L(l)/(hc+hi));
    z = (z_offset+(0:(hc+hi):L(l)))/1000;
    for j=1:size(Rout,2)
        N_in = round((Rout(j)-Rin)/(wc+wi));
        radii = (Rin:(wc+wi):Rout(j))/1000;
        L_out_of_plane = 0;
        for k=1:N_out
            for i=1:N_in
                L_out_of_plane = L_out_of_plane + 2*pi*radii(i);
            end
        end
        R_out_of_plane(j,l) = ro*L_out_of_plane/(wc*hc/1000^2);
        V_out_of_plane(j,l) = R_out_of_plane(j,l)*Imax/1000;
        if V_out_of_plane(j,l) < Vmax
            for k=1:N_out
                for i=1:N_in
                    Bo_out_of_plane(j,l) = Bo_out_of_plane(j,l) + (((mu0/2)*(Imax/1000))*(radii(i)^2)/(z(k)^2+radii(i)^2)^(3/2))*10000;
                end
            end
        else
            I = Vmax/R_out_of_plane(j,l);
            V_out_of_plane(j,l) = Vmax;
            for k=1:N_out
                for i=1:N_in
                    Bo_out_of_plane(j,l) = Bo_out_of_plane(j,l) + (((mu0/2)*I)*(radii(i)^2)/(z(k)^2+radii(i)^2)^(3/2))*10000;
                end
            end
        end
    end
end
Bo_out_of_plane = C*Bo_out_of_plane;

figure
h = pcolor(L,Rout,Bo_out_of_plane);
set(h, 'LineStyle','none')
colorbar
ylabel('Inductor outer radius (mm)')
xlabel('Inductor length (mm)')
title(strcat('Maximum achievable magnetic field (G) Limit=',num2str(Vmax),' V'))


% In terms of Nloop and Nstack
Nloop = 1:90;
Nstack = 1:90;
Bo_out_of_plane = zeros(size(Nloop,2),size(Nstack,2));
R_out_of_plane = zeros(size(Nloop,2),size(Nstack,2));
V_out_of_plane = zeros(size(Nloop,2),size(Nstack,2));

for l=1:size(Nstack,2)
    N_out = Nstack(l);
    z = (z_offset+(0:(hc+hi):Nstack(l)*(hc+hi)))/1000;
    for j=1:size(Nloop,2)
        N_in = Nloop(j);
        radii = (Rin:(wc+wi):(Rin+(Nloop(j)-1)*(wc+wi)))/1000;
        L_out_of_plane = 0;
        for k=1:N_out
            for i=1:N_in
                L_out_of_plane = L_out_of_plane + 2*pi*radii(i);
            end
        end
        R_out_of_plane(j,l) = ro*L_out_of_plane/(wc*hc/1000^2);
        V_out_of_plane(j,l) = R_out_of_plane(j,l)*Imax/1000;
        if V_out_of_plane(j,l) < Vmax
            for k=1:N_out
                for i=1:N_in
                    Bo_out_of_plane(j,l) = Bo_out_of_plane(j,l) + (((mu0/2)*(Imax/1000))*(radii(i)^2)/(z(k)^2+radii(i)^2)^(3/2))*10000;
                end
            end
        else
            I = Vmax/R_out_of_plane(j,l);
            V_out_of_plane(j,l) = Vmax;
            for k=1:N_out
                for i=1:N_in
                    Bo_out_of_plane(j,l) = Bo_out_of_plane(j,l) + (((mu0/2)*I)*(radii(i)^2)/(z(k)^2+radii(i)^2)^(3/2))*10000;
                end
            end
        end
    end
end
Bo_out_of_plane = C*Bo_out_of_plane;

figure('Renderer', 'painters', 'Position', [400 300 400 300])
h = pcolor(Nstack,Nloop,Bo_out_of_plane);
set(h, 'LineStyle','none')
colorbar
xticks([10 30 50 70 90])
% title('Air-core copper-PLA-based solenoids')
ylabel('Number of In-Plane Loops (-)','fontweight','bold')
xlabel('Number of Stacked Layers (-)','fontweight','bold');

