clc;
clear;
close all;
format longE;

imported_data = importdata('D:\Code review\2025-External Flow Dynamics-Finite Volume Method Modeling\Jetting_velocity_profile.mat');
%imported_data = readmatrix('D:\Code review\2025-External Flow Dynamics-Finite Volume Method Modeling\Veolcity_profile(Interpolation).csv');
R_data = importdata('D:\Code review\2025-External Flow Dynamics-Finite Volume Method Modeling\CFD_validation_data.mat');

%%%%%%%%%%%%%%%%%%%% System configuration %%%%%%%%%%%%%%%%%%%%%
%%%% Material Properties %%%%
sigma=0.047;
mu=0.0157;
rho=1111.4;
l0=10^-6;
u0=1;
%%%% Material Properties %%%%

%%%% Solver Setting %%%%
dt_act=10^-9;
dz_act=0.352112*10^-6;
Rc_act=dz_act/10;
t_fin=110*10^-6; %t_fin=200*10^-6;
N=2091;
%%%% Solver Setting %%%%

%%%% System Construction %%%%
dt=dt_act*u0/l0;
dz=dz_act/l0;
Rc=Rc_act/l0;
Re=rho*l0*u0/mu;
We=rho*(u0^2)*l0/sigma;
z=[0:N-1]*dz;
%%%% System Construction %%%%
%%%%%%%%%%%%%%%%%%%% System configuration %%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% Solving %%%%%%%%%%%%%%%%%%%%%
%%%% Initialization %%%%
A=zeros(1,N-1); %%%%% continuity, staggered
R=zeros(1,N-1); %%%%% continuity, staggered
p=zeros(1,N-1); %%%%% continuity, staggered
u=zeros(1,N); %%%%% momentum, collocated
t=0;
%%%% Initialization %%%%

%%%% Transient Solving %%%%
while t <= t_fin*u0/l0
    A_old=A;
    u_old=u;
    R_old=sqrt(A/pi);
    e_max=1;
    e_iter=0;

    %%%% Internal Loop %%%%
    while (e_max>10^-4)&&(e_iter<2)
        A_ref=A;
        u_ref=u;

        %%%% Spatial Solving %%%%
        for i=2:N-2

            %%%% Area Solving %%%%
            if (R(i-1)>0)||(R(i)>0)||(R(i+1)>0)
                if R(i)~=0
                    Fe=pi/4*((R(i+1)+R(i))^2)*u(i+1);
                else
                    Fe=0;
                end
                if R(i-1)~=0
                    Fw=pi/4*((R(i)+R(i-1))^2)*u(i);
                else
                    Fw=0;
                end

                A(i)=A_old(i)-dt/dz*(Fe-Fw);

                if A(i)<0
                    A(i)=0;
                    R(i)=0;
                else
                    R(i)=sqrt(A(i)/pi);
                end
            end
            %%%% Area Solving %%%%

            %%%% Pressure Solving %%%%
            if R(i)>Rc
                if (R(i-1)>=Rc)&&(R(i)>=Rc)&&(R(i+1)>=Rc)
                    dRdz=(R(i+1)-R(i-1))/2/dz;
                    d2Rdz2=(R(i+1)-2*R(i)+R(i-1))/dz^2;
                    p(i)=(1/R(i)/(1+dRdz^2)^0.5-d2Rdz2/(1+dRdz^2)^1.5)/We;
                elseif (R(i-1)>=Rc)&&(R(i)>=Rc)&&(R(i+1)<Rc)
                    p(i)=p(i-1);
                elseif (R(i-1)<Rc)&&(R(i)>=Rc)&&(R(i+1)>=Rc)
                    p(i)=p(i+1);
                else
                    p(i)=0;
                end

                if i==2
                    p(i)=p(i+1);
                end
            end
            %%%% Pressure Solving %%%%

            %%%% Velocity Solving %%%%
            if i>=3
                if (R(i)>Rc)||(R(i-1)>Rc)
                    if (R(i)+R(i-1))/2>Rc
                        if R(i)<Rc
                            Fe=0;
                        else
                            Fe=A(i)*(u(i)^2-3/Re*(u(i+1)-u(i))/dz);
                        end

                        if R(i-1)<Rc
                            Fw=0;
                        else
                            Fw=A(i-1)*(u(i-1)^2-3/Re*(u(i)-u(i-1))/dz);
                        end

                        u(i)=0.3*u(i)+0.7*((u_old(i)*(R_old(i)+R_old(i-1))^2-4/pi*dt/dz*(Fe-Fw))/(R(i)+R(i-1))^2-dt/dz*(p(i)-p(i-1)));
                    elseif (R(i)<=Rc)&&(R(i-1)>=Rc) % Neumann boundary on droplet edge
                        u(i)=u(i-1);
                    elseif (R(i)<=Rc)&&(R(i+1)>=Rc) % Neumann boundary on droplet edge
                        u(i)=u(i+1);
                    end
                end
            end
            %%%% Velocity Solving %%%%

        end
        %%%% Spatial Solving %%%%

        %%%% Cleansing %%%%
        p_med=p;
        p_med(p_med==0)='';
        p_med=median(p_med);

        for i=2:N-2

            %%%% Mass Conservation %%%%
            if (R(i)<Rc)&&(R(i+1)>Rc)&&(0<R(i-1)<Rc)
                A(i+1)=A(i)+A(i+1);
                R(i+1)=sqrt(A(i+1)/pi);
                A(i)=0;
                R(i)=0;
                p(i)=0;

            elseif (R(i)<Rc)&&(0<R(i+1)<Rc)&&(R(i-1)>Rc)
                A(i-1)=A(i)+A(i-1);
                R(i-1)=sqrt(A(i-1)/pi);
                A(i)=0;
                R(i)=0;
                p(i)=0;

            elseif (R(i)<Rc)&&(R(i+1)>Rc)&&(R(i-1)>Rc)
                A(i+1)=A(i)/2+A(i+1);
                R(i+1)=sqrt(A(i+1)/pi);
                A(i-1)=A(i)/2+A(i-1);
                R(i-1)=sqrt(A(i-1)/pi);
                A(i)=0;
                R(i)=0;
                p(i)=0;

            elseif (R(i)~=0)&&(R(i+1)<Rc)&&(R(i-1)<Rc)
                A(i)=0;
                R(i)=0;
                p(i)=0;

            elseif p(i)>20*p_med
                A(i)=0;
                R(i)=0;
                p(i)=0;
            end
            %%%% Mass Conservation %%%%

            %%%% Velocity Cleansing %%%%
            if (i>=3)&&(R(i-1)==0)&&(R(i)==0)
                u(i)=0;
                u_old(i)=0;
            end
            %%%% Velocity Cleansing %%%%

        end
        %%%% Cleansing %%%%

        %%%% Boundary Conditions %%%%
        if floor(t/u0*l0/imported_data(2,1))<=size(imported_data,1)-1
            n=floor(t/u0*l0/imported_data(2,1))+1;
            r=t/u0*l0/imported_data(2,1)-floor(t/u0*l0/imported_data(2,1));

            if n==1
                u(1:2)=imported_data(n,3)*r/u0;
                R(1)=imported_data(n,2)*r/l0;
            else
                u(1:2)=(imported_data(n-1,3)+(imported_data(n,3)-imported_data(n-1,3))*r)/u0;
                R(1)=(imported_data(n-1,2)+(imported_data(n,2)-imported_data(n-1,2))*r)/l0;
            end

            A(1)=pi*R(1).^2;
        else
            u(1:2)=0;
            R(1)=0;
            A(1)=0;
        end
        %%%% Boundary Conditions %%%%

        %%%% Correction Analysis %%%%
        e_max=sqrt(sum(sqrt(((A_ref-A)./(A_ref+10^-10)).^2)+sum(((u_ref(2:end)-u(2:end))./(u_ref(2:end)+10^-10)).^2))); % Error calculation
        e_iter=e_iter+1;
        %%%% Correction Analysis %%%%

    end
    %%%% Internal Loop %%%%

    %%%% Time Stepping %%%%
    t=t+dt;
    %%%% Time Stepping %%%%

    %%%% Plotting %%%%
    if (mod(round(t/dt), 1000) == 0) || (t == dt)
        figure(1)

        z_um   = z * l0 * 1e6;
        R_um   = R * l0 * 1e6;

        idx    = min(floor(t/dt/1000)+1, size(R_data,1));
        R2_um  = R_data(idx, :) * l0 * 1e6;

        u_plot = u * 10;
        p_plot = p * 10;

        plot(z_um(1:end-1), R_um, 'k');
        hold on
        plot(z_um(1:end),   R2_um, '--r');
        plot(z_um,          u_plot, 'b');
        plot(z_um(1:end-1), p_plot, '--b');
        plot(z_um(1:end-1), -R_um, 'k');
        plot(z_um(1:end),   -R2_um, '--r');
        hold off

        title(['Validation of the transient analysis for a jetting droplet, (time = ' num2str(t/u0*l0) ' s)']);
        axis equal
        axis([min(z_um) max(z_um) -100 100])

        legend('Axi-symmetric slender jet model (1D)', ...
               'CFD validation data (2D)', ...
               'Velocity × 10', ...
               'Pressure × 10');
        legend('boxoff')

        xlabel('Travel distance (\mum)')
        ylabel('Droplet radius (\mum)')

        set(gcf, 'position', [20,20,1000,300])
        grid on
        drawnow;
    end
    %%%% Plotting %%%%

end
%%%% Transient Solving %%%%
%%%%%%%%%%%%%%%%%%%% Solving %%%%%%%%%%%%%%%%%%%%%