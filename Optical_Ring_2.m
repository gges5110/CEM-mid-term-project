clear all   
close all
%% DECLARATION OF MATERIAL
mu = 4*pi*10^-7;    epsilon = 1/(36*pi*10^9);    eta = sqrt(mu/epsilon);
c = 3E08;
delta_s = 10e-9;     delta_t = delta_s/sqrt(2)/c;   N_lambda=1.539*1e-6/delta_s;%1.5353 %1.55
%% DECLARATION OF DIMENSION
Nx = 601;   Ny = 701;   M = 16000;
L = round(1);
start_source = 32;
end_source = 74;
ring_origin_1 = [452.5 205.5];
ring_origin_2 = [452.5 458.5];
ring_radius_in = 77.5;      ring_radius_out = 122.5;
%% INITIALIZE
Ex = zeros(Nx-1,Ny);
Ey = zeros(Nx,Ny-1);
Hz = zeros(Nx-1,Ny-1);
D = ones(Nx,Ny).*delta_t./mu./delta_s;
Cx = ones(Nx-1,Ny).*delta_t./epsilon./delta_s;
Cy = ones(Nx,Ny-1).*delta_t./epsilon./delta_s;
%% ESTABLISH WAVEGUIDE
refrac_index = 6^0.5;
for i = 1:Nx
    for j = 1:Ny       
        if j > 31 && j < 75
            Cx(i,j) = delta_t/epsilon/delta_s/(refrac_index^2);
            Cy(i,j) = delta_t/epsilon/delta_s/(refrac_index^2);
        end         
        if j > 590 && j < 635
            Cx(i,j) = delta_t/epsilon/delta_s/(refrac_index^2);
            Cy(i,j) = delta_t/epsilon/delta_s/(refrac_index^2);
        end 
        if ((i + 0.5 - ring_origin_1(1))^2 + (j - ring_origin_1(2))^2) <= ring_radius_out^2 && ((i + 0.5 - ring_origin_1(1))^2 + (j - ring_origin_1(2))^2) > ring_radius_in^2
            Cx(i,j) = delta_t/epsilon/delta_s/refrac_index^2;
        end
        if ((i - ring_origin_1(1))^2 + (j + 0.5 - ring_origin_1(2))^2) <= ring_radius_out^2 && ((i - ring_origin_1(1))^2 + (j + 0.5 - ring_origin_1(2))^2) > ring_radius_in^2
            Cy(i,j) = delta_t/epsilon/delta_s/refrac_index^2;
        end      
        if ((i + 0.5 - ring_origin_2(1))^2 + (j - ring_origin_2(2))^2) <= ring_radius_out^2 && ((i + 0.5 - ring_origin_2(1))^2 + (j - ring_origin_2(2))^2) > ring_radius_in^2
            Cx(i,j) = delta_t/epsilon/delta_s/refrac_index^2;
        end
        if ((i - ring_origin_2(1))^2 + (j + 0.5 - ring_origin_2(2))^2) <= ring_radius_out^2 && ((i - ring_origin_2(1))^2 + (j + 0.5 - ring_origin_2(2))^2) > ring_radius_in^2
            Cy(i,j) = delta_t/epsilon/delta_s/refrac_index^2;
        end    
    end
end
factor = c*delta_t/delta_s;
%% UPDATE EQUATION
for m = 1:M
    %Update Hz
    Hz(1:Nx-1,1:Ny-1) = Hz(1:Nx-1,1:Ny-1) +  D(1:Nx-1,1:Ny-1).*(Ex(1:Nx-1,2:Ny) - Ex(1:Nx-1,1:Ny-1) + Ey(1:Nx-1,1:Ny-1) - Ey(2:Nx,1:Ny-1));                                                                                                         
    %Special update in scat/tot field boundary
    time_E =  (m-1)*delta_t-(L-1)*delta_s/c;    
    Ey_inc = .5*sin(((2*pi*(c/(delta_s*N_lambda))*(m-1)*delta_t)));
    Hz(L,start_source:end_source) = Hz(L,start_source:end_source) + D(L,start_source:end_source).*( ones(1,end_source - start_source + 1)* Ey_inc ); 
    %ABC
    Ey(1,1:Ny-1)= Ey(1,1:Ny-1)*(1-factor)+factor*Ey(2,1:Ny-1);
    Ey(Nx,1:Ny-1)=Ey(Nx,1:Ny-1)*(1-factor)+factor*Ey(Nx-1,1:Ny-1);
    Ex(1:Nx-1,1)=Ex(1:Nx-1,1)*(1-factor)+factor*Ex(1:Nx-1,2); 
    Ex(1:Nx-1,Ny)=Ex(1:Nx-1,Ny)*(1-factor)+factor*Ex(1:Nx-1,Ny-1);    
    
    %Update Ex and Ey
    Ex(1:Nx-1,2:Ny-1) = Ex(1:Nx-1,2:Ny-1) + Cx(1:Nx-1,2:Ny-1).*(Hz(1:Nx-1,2:Ny-1) - Hz(1:Nx-1,1:Ny-2));                                   
    Ey(2:Nx-1,1:Ny-1) = Ey(2:Nx-1,1:Ny-1) - Cy(2:Nx-1,1:Ny-1).*(Hz(2:Nx-1,1:Ny-1) - Hz(1:Nx-2,1:Ny-1));                                   
    %Special update in scat/tot field boundary
    time_H = (m-0.5)*delta_t-(L-0.5)*delta_s/c;
    Hz_inc = .5*sin(((2*pi*(c/(delta_s*N_lambda))*(m-1)*delta_t)))/eta;
    Ey(L,start_source:end_source) = Ey(L,start_source:end_source) +Cy(L,start_source:end_source).*( ones(1,end_source - start_source + 1)*Hz_inc );     
    %Plot
    h=imagesc(delta_s*(1:Ny)*10^6,(delta_s*(1:Nx)*10^6),Hz,[-.05,.05]);
    name = ['time step = ' num2str(m)];
    title(name);    
    pause(.01)
end
