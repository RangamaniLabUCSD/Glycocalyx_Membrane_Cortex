clc;
clear;
close all;

%% Define parameters
%% membrane properties
kappa = 20;  %% bending rigidity, Unit: k_BT
sigma = 0.012;  %% membrane tension, Unit: k_BT/nm^2
lambda = 0.2;   %% Line tension, Unit: k_BT/nm
c_0 = 0.02;   %% spontaneous curveture, Unit: nm^(-1)
k = 0.025;  %% Unit: k_BT/nm^2
l_c = 5;  %% Unit: nm
rho_0 = 1*1e-4;   %% Unit: nm^-2
%% Polymer properties
a = 10;      %% effective monomer size,  Unit: nm
v = a^3;     %% Excluded volume parameter, Unit: nm^3
R_0 = 100;        %%   membrane patch radius,  Unit: nm
A_s = pi*R_0^2;   %%   membrane patch area,  Unit: nm^2

N_t = 20;
t = 0;
u = 0;
for N = N_t
    N
    u = u+1;
    R_G = (1/sqrt(6))*a*N^(3/5);  %% radius of gyration, Unit: nm
    R_F= a*N^(3/5);  %% Flory radius, Unit: nm

    xi_j = linspace(10,30,101);
    j = 0;
    for xi = xi_j
        xi
        j = j+1;
        S_R = xi^2;     
        rho_g = 1/xi^2;  
        rho_anchor = 1/xi^2;  
        N_p = A_s/S_R;   
        if xi<2*(1/sqrt(6))*a*N^(3/5)
            brush = 1;
        else
            brush = 0;
        end

        eta_span = 0.0:0.0005:1; %% fraction η which characterize the area ratio between a spherical cap and a sphere.
        i = 0;
        F_tot = 0;
        F_09_1_tot = 0;
        m = 0;

        for eta = eta_span
            i = i+1;
            H_flat = N*S_R^(-1/3)*(v*a^2/3)^(1/3);
            F_flat1 = N_p*N*(9/2)*(v/3)^(2/3)*(S_R*a)^(-2/3);
            F_flat2 = (9/2)*(pi*R_0^2/xi^2)*N*(v/(3*a*xi^2))^(2/3);
            if eta==0
                F(i) = ((9/2)*(pi*R_0^2/xi^2)*N*(v/(3*a*xi^2))^(2/3)+8*pi*kappa*(sqrt(eta)-R_0*c_0/4)^2+pi*sigma*eta*R_0^2+2*pi*lambda*R_0*sqrt(1-eta)...
                    +eta*pi*k*rho_0*R_0^4/6)/(pi*kappa);
            else
                R = R_0/(2*sqrt(eta));
                eta_c = (l_c/R_0)^2;
                H_spherical = R*(1+5/3*H_flat/R)^(3/5)-R;
                if 0<eta && eta<=eta_c
                    F(i) = ((9/2)*(pi*R_0^2/xi^2)*(R_0/(2*sqrt(eta)))*(3*v^(1/2)/(xi*a^2))^(2/3)*((1+10*N*sqrt(eta)/3/R_0*(v*a^2/3/xi^2)^(1/3))^(1/5)-1)...
                    + 8*pi*kappa*(sqrt(eta)-R_0*c_0/4)^2 + pi*sigma*eta*R_0^2 + 2*pi*lambda*R_0*sqrt(1-eta)...
                    + eta*pi*k*rho_0*R_0^4/6)/(pi*kappa);
                else
                    F(i) = ((9/2)*(pi*R_0^2/xi^2)*(R_0/(2*sqrt(eta)))*(3*v^(1/2)/(xi*a^2))^(2/3)*((1+10*N*sqrt(eta)/3/R_0*(v*a^2/3/xi^2)^(1/3))^(1/5)-1)...
                    + 8*pi*kappa*(sqrt(eta)-R_0*c_0/4)^2 + pi*sigma*eta*R_0^2 + 2*pi*lambda*R_0*sqrt(1-eta)...
                    + k/2*pi*rho_0*l_c^2*R_0^2*(1-l_c/(R_0*sqrt(eta))) + 1/6*pi*k*rho_0*R_0*l_c^3/sqrt(eta))/(pi*kappa);
                end
                F_tot = F_tot + exp(-F(i));
            end

            if eta>=1 && eta<=1
                m = m+1;
                if eta==0
                    F_09_1(m) = ((9/2)*(pi*R_0^2/xi^2)*N*(v/(3*a*xi^2))^(2/3)+8*pi*kappa*(sqrt(eta)-R_0*c_0/4)^2+pi*sigma*eta*R_0^2+2*pi*lambda*R_0*sqrt(1-eta)...
                        +eta*pi*k*rho_0*R_0^4/6)/(pi*kappa);
                else
                    R = R_0/(2*sqrt(eta));
                    eta_c = (l_c/R_0)^2;
                    if 0<eta && eta<=eta_c
                        F_09_1(m) = ((9/2)*(pi*R_0^2/xi^2)*(R_0/(2*sqrt(eta)))*(3*v^(1/2)/(xi*a^2))^(2/3)*((1+10*N*sqrt(eta)/3/R_0*(v*a^2/3/xi^2)^(1/3))^(1/5)-1)...
                            + 8*pi*kappa*(sqrt(eta)-R_0*c_0/4)^2 + pi*sigma*eta*R_0^2 + 2*pi*lambda*R_0*sqrt(1-eta)...
                            + eta*pi*k*rho_0*R_0^4/6)/(pi*kappa);
                    else
                        F_09_1(m) = ((9/2)*(pi*R_0^2/xi^2)*(R_0/(2*sqrt(eta)))*(3*v^(1/2)/(xi*a^2))^(2/3)*((1+10*N*sqrt(eta)/3/R_0*(v*a^2/3/xi^2)^(1/3))^(1/5)-1)...
                            + 8*pi*kappa*(sqrt(eta)-R_0*c_0/4)^2 + pi*sigma*eta*R_0^2 + 2*pi*lambda*R_0*sqrt(1-eta)...
                            + k/2*pi*rho_0*l_c^2*R_0^2*(1-l_c/(R_0*sqrt(eta))) + 1/6*pi*k*rho_0*R_0*l_c^3/sqrt(eta))/(pi*kappa);
                    end
                    F_09_1_tot = F_09_1_tot + exp(-F_09_1(m));
                end
            end
        F = F';
        [y,index_min] = min(F);
        eta_min(j) = eta_span(index_min);
        E{j} = F;
        leg{j}=strcat('\xi=',num2str(xi),'nm');
        end
        t = t+1;
        Results(t,:) = [xi,rho_g,N_p,N,eta_span(index_min),lambda,c_0,rho_0,kappa,sigma,k,l_c,a,R_0];
    end
end

subplot(2,2,1);
for j=1:length(xi_j)
    plot(eta_span,E{j},'linewidth',2)
    hold on
end
hold off
legend(leg)
xlabel('\eta');
ylabel('\itF_{tot}/\pi\kappa')

subplot(2,2,2);
plot(Results(:,2),eta_min,'linewidth',2)
xlabel('\it\rho');
ylabel('\it\eta_{min}')

subplot(2,2,3);
% plot(N_j,P_001,'--o',N_j,P_191,'r','linewidth',2)
scatter(Results(:,3),Results(:,4),[],Results(:,5),'filled')
u=colorbar;
set(u,'FontName','Times New Roman','FontSize',15,'linewidth',2.5,'FontWeight','bold');
set(get(u,'title'),'string','\eta');
% u.Label.String = '\eta';
set(gca,'xscale','log');
xlim([30 330]);
xlabel('\itN_p');
ylabel('\itN')
box on
