%% Alfa = 0
% define the mean chord function

np=300;
k=1;
%np=160;
    n=np(k)+1;% # of nodes
    delta_c=1/np(k); %length of each panel
    
    x=linspace(0,1,n)';
    n=length(x);
    z=zeros(n,1);
    % sum of the three last digits of my student number (240)=6
    %Reference Airfoil NACA 2408
    m=2/100;
    p=4/10;
    
    for i=1:n
        if x(i)<=p
            z(i)=m/(p^2)*(2*p*x(i)-(x(i))^2);
        else
            z(i)=m/((1-p)^2)*(1-2*p+2*p*x(i)-(x(i))^2);
        end
    end
    
    % plot mean camber line 
    % figure(1)
    % plot(x,z)
    % axis equal
    % axis on
    % grid on
    % title('NACA 2408 mean camber line')
    % xlabel('$x/c$',Interpreter='latex')
    % ylabel('$z/c$',Interpreter='latex')
    
    
    
    m=zeros(np(k),1); %slope of each panel
    for i=1:np(k)
            %central difference
            m(i)=(z(i+1)-z(i))/(delta_c);
    end
    
    alpha_i=atan(-m); %relative angle of each panel
    
    % z coordinate of each vortex
    z_i=zeros(np(k),1);
    for i=1:np(k)
        z_i(i)=(z(i)+z(i+1))/2; %average of two nodes
    end
    
    
    % sum of the three last digits of my student number (240)=6
    %Reference Airfoil NACA 2408
    m=2/100;
    p=4/10;
    %local variables for the position of vortices and zero normal velocity
    x_vort=zeros(np(k),1);
    x_norm=zeros(np(k),1);
    z_vort=zeros(np(k),1);
    z_norm=zeros(np(k),1);
    for i=1:np(k)
        x_vort(i)=delta_c*((i-1)+1/4);%quarter chord
        x_norm(i)=delta_c*((i-1)+3/4);%three-quarter cord
        %evaluate on the camber line function
        z_vort(i)=chord(m,p,x_vort(i));
        z_norm(i)=chord(m,p,x_norm(i));  
    end
    
    %find the coefficients of the matrix
    A=zeros(np(k),np(k));
    for i=1:n-1
        for j=1:n-1
            x_d=x_norm(i)-x_vort(j);
            z_d=z_norm(i)-z_vort(j);
            r_sq=x_d^2+z_d^2;
            u=z_d/(2*pi*r_sq);
            w=-x_d/(2*pi*r_sq);
            %dot product
            A(i,j)=(u*sin(alpha_i(i)))+(w*cos(alpha_i(i)));
        end
    end
    

    c=linspace(0,1,np(k))';
    
    
    AoA=0;
    i=1;
            b=-sin(AoA(i)+alpha_i);%define b vector
            Gamma=A\b ;%solve matrix system
            L_p=sum(Gamma);%sum all circulations
            C_l=2*L_p; %after cancelling on the C_l formula
        
            %moment at LE
            M_LE=sum(Gamma.*(x_vort.*cos(AoA(i))+z_vort.*sin(AoA(i))));
            %chord normalized center of pressure
            x_cop=M_LE/L_p;
            %moment coefficient about the quarter-chord
            M_c4=-(x_cop-(1/4))*L_p;
            C_mc4=2*M_c4;
            Cp=2*Gamma/delta_c;
      
    
    % figure(3)
    % yyaxis left
    % plot(rad2deg(AoA),C_l,'*-')
    % ylabel('$C_l$',Interpreter='latex')
    % yyaxis right
    % plot(rad2deg(AoA),C_mc4,'^-')
    % ylabel('$C_{m, c/4}$',Interpreter='latex')
    % grid on
    % title('NACA 2408')
    % xlabel('$\alpha (^\circ)$',Interpreter='latex')
    % 
    % 




%% Panel Sweep

clear 
clc

np=[5,10,20,50,100,150,200];% numer of panels
dcl_da=zeros(1,length(np)); % derivative of the lift gradient

% For loop to iterate all calculations for different numbers of panels
for k=1:length(np)
    n=np(k)+1;% # of nodes
    delta_c=1/np(k); %length of each panel
   
    x=linspace(0,1,n)'; % x coordinate
    n=length(x);
    z=zeros(n,1); % z coordinate

    % sum of the three last digits of my student number (240)=6
    %Reference Airfoil NACA 2408
    m=2/100;
    p=4/10;

    % define the mean chord function
    for i=1:n
        if x(i)<=p
            z(i)=m/(p^2)*(2*p*x(i)-(x(i))^2);
        else
            z(i)=m/((1-p)^2)*(1-2*p+2*p*x(i)-(x(i))^2);
        end
    end
    
    %plot mean camber line 
    figure(1)
    plot(x,z)
    axis equal
    axis on
    grid on
    title('NACA 2408 mean camber line')
    xlabel('$x/c$',Interpreter='latex')
    ylabel('$z/c$',Interpreter='latex')
    
    
    m=zeros(np(k),1); %slope of each panel
    for i=1:np(k)
            %central difference
            m(i)=(z(i+1)-z(i))/(delta_c);
    end
    
    alpha_i=atan(-m); %relative angle of each panel
    
    % z coordinate of each vortex
    z_i=zeros(np(k),1);
    for i=1:np(k)
        z_i(i)=(z(i)+z(i+1))/2; %average of two nodes
    end
    
    m=2/100;
    p=4/10;
    %local variables for the position of vortices and zero normal velocity
    x_vort=zeros(np(k),1);
    x_norm=zeros(np(k),1);
    z_vort=zeros(np(k),1);
    z_norm=zeros(np(k),1);
    for i=1:np(k)
        x_vort(i)=delta_c*((i-1)+1/4);%quarter chord
        x_norm(i)=delta_c*((i-1)+3/4);%three-quarter cord
        %evaluate on the camber line function
        z_vort(i)=chord(m,p,x_vort(i));
        z_norm(i)=chord(m,p,x_norm(i));  
    end
    
    %Determine the coefficients of the matrix
    A=zeros(np(k),np(k));
    for i=1:n-1
        for j=1:n-1
            x_d=x_norm(i)-x_vort(j);
            z_d=z_norm(i)-z_vort(j);
            r_sq=x_d^2+z_d^2;
            u=z_d/(2*pi*r_sq);
            w=-x_d/(2*pi*r_sq);
            %dot product
            A(i,j)=(u*sin(alpha_i(i)))+(w*cos(alpha_i(i)));
        end
    end
   
    % AoA scan
    C_l=zeros(1,10); % Lift coeffs.
    C_mc4=zeros(1,10); % Moment coeffs.
    AoA=deg2rad(linspace(0,10,10));% angle of attack

        for i=1:length(AoA)
            b=-sin(AoA(i)+alpha_i);%define b vector
            Gamma=A\b ;%solve matrix system
            L_p=sum(Gamma);%sum all circulations
            C_l(i)=2*L_p; % Lift Coefficient (after cancelling on the C_l formula)
            %moment at LE
            M_LE=sum(Gamma.*(x_vort.*cos(AoA(i))+z_vort.*sin(AoA(i))));
            %chord normalized center of pressure
            x_cop=M_LE/L_p;
            %moment coefficient about the quarter-chord
            M_c4=-(x_cop-(1/4))*L_p;
            C_mc4(i)=2*M_c4;
            %pressure distribution
            Cp=2*Gamma/delta_c;
        end
  
    %computes the derivative of the lift gradient dC_l/da
    dcl_da(k)=diff(C_l)/diff(AoA);
end

%%
theo=linspace(2*pi,2*pi,7);
%plot(np,theo)
hold on
plot(np,dcl_da,'o-')
grid on
xlabel('Number of Panels')
ylabel('$dC_l/d\alpha$',Interpreter='latex')
%ylim([6.2 6.3])

%% Fixed Panel

clear 
clc

np=300;% numer of panels
dcl_da=zeros(1,length(np)); % derivative of the lift gradient

% For loop to iterate all calculations for different numbers of panels
for k=1:length(np)
    n=np(k)+1;% # of nodes
    delta_c=1/np(k); %length of each panel
   
    x=linspace(0,1,n)'; % x coordinate
    n=length(x);
    z=zeros(n,1); % z coordinate

    % sum of the three last digits of my student number (240)=6
    %Reference Airfoil NACA 2408
    m=2/100;
    p=4/10;

    % define the mean chord function
    for i=1:n
        if x(i)<=p
            z(i)=m/(p^2)*(2*p*x(i)-(x(i))^2);
        else
            z(i)=m/((1-p)^2)*(1-2*p+2*p*x(i)-(x(i))^2);
        end
    end
    
    %plot mean camber line 
    figure(1)
    plot(x,z)
    axis equal
    axis on
    grid on
    title('NACA 2408 mean camber line')
    xlabel('$x/c$',Interpreter='latex')
    ylabel('$z/c$',Interpreter='latex')
    
    
    m=zeros(np(k),1); %slope of each panel
    for i=1:np(k)
            %central difference
            m(i)=(z(i+1)-z(i))/(delta_c);
    end
    
    alpha_i=atan(-m); %relative angle of each panel
    
    % z coordinate of each vortex
    z_i=zeros(np(k),1);
    for i=1:np(k)
        z_i(i)=(z(i)+z(i+1))/2; %average of two nodes
    end
    
    m=2/100;
    p=4/10;
    %local variables for the position of vortices and zero normal velocity
    x_vort=zeros(np(k),1);
    x_norm=zeros(np(k),1);
    z_vort=zeros(np(k),1);
    z_norm=zeros(np(k),1);
    for i=1:np(k)
        x_vort(i)=delta_c*((i-1)+1/4);%quarter chord
        x_norm(i)=delta_c*((i-1)+3/4);%three-quarter cord
        %evaluate on the camber line function
        z_vort(i)=chord(m,p,x_vort(i));
        z_norm(i)=chord(m,p,x_norm(i));  
    end
    
    %Determine the coefficients of the matrix
    A=zeros(np(k),np(k));
    for i=1:n-1
        for j=1:n-1
            x_d=x_norm(i)-x_vort(j);
            z_d=z_norm(i)-z_vort(j);
            r_sq=x_d^2+z_d^2;
            u=z_d/(2*pi*r_sq);
            w=-x_d/(2*pi*r_sq);
            %dot product
            A(i,j)=(u*sin(alpha_i(i)))+(w*cos(alpha_i(i)));
        end
    end
   
    % AoA scan
    C_l=zeros(1,10); % Lift coeffs.
    C_mc4=zeros(1,10); % Moment coeffs.
    AoA=deg2rad(linspace(-5,10,20));% angle of attack

        for i=1:length(AoA)
            b=-sin(AoA(i)+alpha_i);%define b vector
            Gamma=A\b ;%solve matrix system
            L_p=sum(Gamma);%sum all circulations
            C_l(i)=2*L_p; % Lift Coefficient (after cancelling on the C_l formula)
            %moment at LE
            M_LE=-sum(Gamma.*(x_vort.*cos(AoA(i))+z_vort.*sin(AoA(i))));
            %chord normalized center of pressure
            %moment coefficient about the quarter-chord
            C_mc4(i)=2*M_LE+ C_l(i)/4;
            %pressure distribution
            Cp=2*Gamma/delta_c;
        end
  
    %computes the derivative of the lift gradient dC_l/da
    dcl_da(k)=diff(C_l)/diff(AoA);
end

%% Polar plot

% Experimental from Abott
alpha=[-4,-2,0,2,4,6,8,10];
Cl_exp=[-0.2,0,0.2,0.4,0.6,0.8,1,1.2];

AoA=rad2deg(AoA);
figure(2)
plot(AoA,C_l,'-o')
grid on
grid minor
hold on
xfoil=readtable('naca2408_xfoil.csv');
plot(xfoil.alpha,xfoil.CL,'-x')
plot(alpha,Cl_exp,'*',Color='magenta')
%ylim([-1 2])
xlim([-11 11])

xlabel('Angle of Attack [°]')
ylabel('$C_l$',Interpreter='latex')

xline(0, 'Color', 'k'); % Draw line for Y axis.
yline(0, 'Color', 'k'); % Draw line for X axis.
legend('Method of panels (n=300)', 'Xfoil (n=300)','Experimental data','','')
title('Lift Coefficient vs. AoA for NACA 2408')
%% Cm plot


% Experimental from Abott
alpha=[-4,-2,0,2,4,6,8,10];
Cm_exp=[-0.05,-0.05,-0.05,-0.05,-0.0495,-0.049,-0.0485,-0.048];

figure(3)
plot(AoA,C_mc4,'-o')
grid on
hold on
xfoil=readtable('naca2408_xfoil.csv');
plot(xfoil.alpha,xfoil.CM,'-x')
plot(alpha,Cm_exp,'*',Color='magenta')
ylim([-0.1,0])
xlim([-6,11])
legend('Method of panels (n=300)', 'Xfoil (n=300)','Experimental data')
xlabel('Angle of Attack [°]')
ylabel('$C_m$',Interpreter='latex')
title('$C_m$ vs. AoA for NACA 2408',Interpreter='latex')
%% Pressure plot

data=readtable("cp_2408.txt");
dcp=(-data.Var2(1:150)-flip(data.Var2(151:300)))/2;

plot(x(1:300),Cp,Color='red',LineWidth=1)
grid on
hold on
plot(data.Var1,-data.Var2,'--',Color='blue')
plot(data.Var1(1:150),dcp,Color='blue',LineWidth=1);
ylim([-0.5 0.5])
xlabel('$x/c$',Interpreter='latex')
ylabel('$\Delta C_p$',Interpreter='latex')
yline(0, 'Color', 'k'); % Draw line for X axis.
legend('Method of panels', 'Xfoil','Xfoil mean')
title('Pressure differential NACA 2408 (300 panels)')

%% Other airfoils

figure(1)
xfoil=readtable('naca2408_xfoil.csv');
plot(xfoil.alpha,xfoil.CL,'-x',LineWidth=1)
grid on
grid minor
hold on

xfoil=readtable('naca4408_xfoil.csv');
plot(xfoil.alpha,xfoil.CL,'-o',LineWidth=1)

xfoil=readtable('naca2208_xfoil.csv');
plot(xfoil.alpha,xfoil.CL,'-^',LineWidth=1)

%ylim([-1 2])
xlim([-11 11])

xlabel('Angle of Attack [°]')
ylabel('$C_l$',Interpreter='latex')

xline(0, 'Color', 'k'); % Draw line for Y axis.
yline(0, 'Color', 'k'); % Draw line for X axis.
legend('NACA 2408','NACA 4408 (increased camber)','NACA 2208 (decreased location of max. camber)')
title('Lift Coefficient vs. AoA')

%%

figure(1)
xfoil=readtable('naca2408_xfoil.csv');
plot(xfoil.alpha,xfoil.CM,'-x',LineWidth=1)
grid on
grid minor
hold on

xfoil=readtable('naca4408_xfoil.csv');
plot(xfoil.alpha,xfoil.CM,'-o',LineWidth=1)

xfoil=readtable('naca2208_xfoil.csv');
plot(xfoil.alpha,xfoil.CM,'-^',LineWidth=1)

ylim([-0.3 0])
xlim([-6 11])

xlabel('Angle of Attack [°]')
ylabel('$C_m$',Interpreter='latex')

xline(0, 'Color', 'k'); % Draw line for Y axis.
yline(0, 'Color', 'k'); % Draw line for X axis.
legend('NACA 2408','NACA 4408 (increased camber)','NACA 2208 (decreased location of max. camber)')
title('Moment Coefficient vs. AoA')