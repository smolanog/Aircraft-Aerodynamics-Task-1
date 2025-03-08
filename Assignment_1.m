%%

clear 
clc
% define the mean chord function

np=[10,25,50,100,300,500,1000];% numer of panels
%np=160;
avg=zeros(1,length(np));
for k=1:length(np)
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
    
    
    AoA=deg2rad(4);
    b=-sin(AoA+alpha_i);%define b vector
    Gamma=A\b ;%solve matrix system
    L_p=sum(Gamma);%sum all circulations
    c=linspace(0,1,np(k))';
    C_p=2*Gamma./delta_c; %after cancelling on the C_l formula
    
    
    
    % AoA scan
    C_l=zeros(1,10);
    C_mc4=zeros(1,10);
    AoA=deg2rad(linspace(0,10,10));% angle of attack
    %AoA=deg2rad(0);
        for i=1:length(AoA)
            b=-sin(AoA(i)+alpha_i);%define b vector
            Gamma=A\b ;%solve matrix system
            L_p=sum(Gamma);%sum all circulations
            C_l(i)=2*L_p; %after cancelling on the C_l formula
        
            %moment at LE
            M_LE=sum(Gamma.*(x_vort.*cos(AoA(i))+z_vort.*sin(AoA(i))));
            %chord normalized center of pressure
            x_cop=M_LE/L_p;
            %moment coefficient about the quarter-chord
            M_c4=-(x_cop-(1/4))*L_p;
            C_mc4(i)=2*M_c4;
            Cp=2*Gamma/delta_c;
        end
    
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

    % figure(4)
    % plot(c,Cp)
    % grid on
    
  
    %computes the derivative of the lift gradient dC_l/da
    avg(k)=diff(C_l)/diff(AoA);
end

%%
theo=linspace(2*pi,2*pi,7);
%plot(np,theo)
hold on
plot(np,avg,'o-')
grid on
xlabel('Number of Panels')
ylabel('$dC_l/d\alpha$',Interpreter='latex')
%ylim([6.2 6.3])



%% Comparison of results

%%

clear 
clc
% define the mean chord function

np=1000;
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

    figure(4)
    plot(c,Cp)
    grid on
    ylim([-2 1])
    









