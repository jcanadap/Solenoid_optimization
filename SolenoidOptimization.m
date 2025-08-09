clear; clc;
close all;

% Physical parameters
mu0 = 1.26/1000000; % Vacuum permeability (H/m)
ro = 0.00032; % Resistivity (Ohm*m)
Imax = 35; % Maximum recommended current (mA)
Pmax = 180; % Maximum suppliable power (W)

% Geometrical parameters
Rin = 2; % Inner radius (mm)
wc = 0.6; % Conductive trace width (mm)
hc = 0.6; % Conductive trace height (mm)
wi = 0.4; % Insulation width (mm)
hi = 0.05; % Insulation height (mm)
z_offset = 0.7; % Solenoid to probe distance (mm)

% Parameter sweep
% In terms of actual size
Rout = Rin:1:30; % Outer radius (mm)
L = 0.3:0.025:20; % Inductor length (mm)
% In terms of the number of loops and stacked layers
Nloop = 1:90; % Number of in-plane loops
Nstack = 1:90; % Number of stacked layers

% Implementation
Vmax = Pmax/(Imax/1000);

% In terms of actual size
Bo_out_of_plane = zeros(size(Rout,2),size(L,2));
R_out_of_plane = zeros(size(Rout,2),size(L,2));
V_out_of_plane = zeros(size(Rout,2),size(L,2));

for l=1:size(L,2)
    N_out = round(L(l)/(hc+hi));
    z = (z_offset+(0:(hc+hi):L(l)))/1000;
    for j=1:size(Rout,2)
        N_in = round((Rout(j)-Rin)/(wc+wi));
        radii = (Rin:(wc+wi):Rout(j))/1000;
        L_out_of_plane = N_out*sum(2*pi*radii);
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

figure
h = pcolor(L,Rout,Bo_out_of_plane);
set(h, 'LineStyle','none')
colorbar
ylabel('Inductor outer radius (mm)')
xlabel('Inductor length (mm)')
title(strcat('Maximum achievable magnetic field (G) Limit=',num2str(Vmax),' V'))


% In terms of Nloop and Nstack
Bo_out_of_plane = zeros(size(Nloop,2),size(Nstack,2));
R_out_of_plane = zeros(size(Nloop,2),size(Nstack,2));
V_out_of_plane = zeros(size(Nloop,2),size(Nstack,2));

for l=1:size(Nstack,2)
    N_out = Nstack(l);
    z = (z_offset+(0:(hc+hi):Nstack(l)*(hc+hi)))/1000;
    for j=1:size(Nloop,2)
        N_in = Nloop(j);
        radii = (Rin:(wc+wi):(Rin+(Nloop(j)-1)*(wc+wi)))/1000;
        L_out_of_plane = N_out*sum(2*pi*radii);
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

figure
% figure('Renderer', 'painters', 'Position', [400 300 400 300])
h = pcolor(Nstack,Nloop,Bo_out_of_plane);
set(h, 'LineStyle','none')
colorbar
xticks([10 30 50 70 90])
title(strcat('Maximum achievable magnetic field (G) Limit=',num2str(Vmax),' V'))
ylabel('Number of In-Plane Loops (-)','fontweight','bold');
xlabel('Number of Stacked Layers (-)','fontweight','bold');

