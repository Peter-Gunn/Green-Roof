 function SSD = FD_Model(y)

%Boundaries for "str" and "fin" for monthly simulations
%     may = [2905:3648];
%     jun = [3649:4368];
%     jul = [4369:5112];
%     aug = [5113:5856];
%     sep = [5857:6576];

%Set temporal boundaries for calibration (see above)
%Here, only change the 4-digit number in line 13 to a starting value from
%above, i.e. only change "2905" and keep the "-1-48" in place. 
str = 2905-1-48; 
fin = 3648; 

%READ IN HOURLY DATA
num = readmatrix('hourly2016');
Troom = num(str:fin,7+1);
Tmat = num(str:fin,13+1);
sRad = num(str:fin,1+1);
rain = num(str:fin,2+1);
Tair = num(str:fin,3+1);
Tcan = num(str:fin,4+1);
Tsoil = num(str:fin,5+1);
TConc = num(str:fin,6+1);
VWC = num(str:fin,8+1);
km = num(str:fin,10+1);
Cm = num(str:fin,11+1);
Snow = num(str:fin,12+1);
Patm = num(str:fin,14+1);
nRad = num(str:fin,15+1);
RH = num(str:fin,16+1);
wind = num(str:fin,17+1);
TroomCR = num(str:fin,18+1);
TmatCR = num(str:fin,19+1);
Tsoil1 = num(str:fin,20+1);
Tsoil2 = num(str:fin,21+1);
Tsoil3 = num(str:fin,22+1);
VWC1 = num(str:fin,23+1);
VWC2 = num(str:fin,24+1);
VWC3 = num(str:fin,25+1);
rVWC = fillmissing(VWC,'linear');

Ca = 1*(10^-3); %MJ/kg K

%CONCRETE LAYER
Rc = y(7); Lc = 0.3; kc = Lc/Rc; pCc = y(8);

%CANOPY LAYER
% La = 0.1; ka = 5 ; pa = 1.225;  Cpa = 1005; Ra = La/ka; %K x 0.3 makes soil temp close ka = 5
Ra = y(1);  
La = 0.1;  
ka = La/Ra; 
pCa = y(2)*4800*580+(1-y(2))*1.184*1007;

%SOIL LAYER
Ls = 0.1;
Rs = y(4);
ks = Ls/Rs;
ps = 1101; 
pCs = y(3);

%INSULATION LAYER
Ri = y(6); Li = 0.1; ki = Li/Ri; pi = 45; Cpi = 1300;

%MODEL PROPERTIES
L = La+Ls+Li+Lc; %Total depth of roof envelope
numUnknowns = 120; %Total number of nodes

numBlocks = numUnknowns; 
dx = L/(numUnknowns); 
numNodes = numBlocks;  
Tcan0 = Tcan(2,1); 
Troom0 = Troom(2,1); 
Tmat0 = Tmat(2,1);
dt = 3600 ; N = numNodes;
Dur = length(Tcan)-1;
mDiv = 12;

%FOR POPULATING A MATRIX
fair = La/L; fsoil = Ls/L; fins = Li/L; fconc = Lc/L;
nair = fair*(numBlocks); nsoil = fsoil*(numBlocks); nins = fins*(numBlocks); nconc = fconc*(numBlocks);

%EXTERIOR CONVECTION COEFFICIENT
for hn = 2:length(Tsoil)   
    hext(hn,1) = (5.9+4.1*wind(hn,1)*((511+294)/(511+Tair(hn,1)+273)));
end


%INITIAL CONDITION - STEADY STATE TEMP. PROFILE
khm = (L)/(((La)/ka)+(Ls/ks)+(Li/ki)+(((Lc))/kc));
Qss = (khm/(L))*(Tcan0-Troom0); 
dTsa = (Qss*dx)/ka; dTss = (Qss*dx)/ks; dTsi = (Qss*dx)/ki; dTsc = (Qss*dx)/kc;

%FOR POPULATING A MATRIX pt.2
Mmain = transpose(ones(1,N-1));  Msup = transpose(ones(1,N-1));  Msub = transpose(ones(1,N-1));  

%INITIAL CONDITION - STEADY STATE TEMP. PROFILE pt.2
LocTCIint = round((La+Ls+Li)/dx)+1;
locTsoil = round((La+Ls*0.5)/dx);
locTmat = round((La+Ls)/dx);
locTcan = 1;
yx1 = transpose(linspace(dx,L-dx,N-2)); %For loop below
for n = 1:length(yx1)
    if n == 1
        b1(n,1) = Tcan0-dTsa;%+Qsurf(1,1)*(dx/ka)-dTsa;
    elseif yx1(n-1,1) < La
        b1(n,1) = b1(n-1,1)-dTsa;
    elseif yx1(n-1,1) < (La+Ls)
        b1(n,1) = b1(n-1,1)-dTss;
    elseif yx1(n-1,1) < (La+Ls+Li)
        b1(n,1) = b1(n-1,1)-dTsi;
    elseif yx1(n-1,1) < (La+Ls+Li+Lc)
        b1(n,1) = b1(n-1,1)-dTsc;
    end
end

%generating initial temperature profile of known and unknow portion of
%solution space
b12 = ones(N,1);% *Troom0; %FIXED TEMP

b12(1,1) = Tcan0;      %For Changing Boundaries
b12(2:N-1)=b1(:,1);      %For Changing Boundaries
b12(N,1) = b12(N-1,1)-dTsc;   %For Changing Boundaries

%ENABLE LINES 127-135 TO PLOT INITIAL CONDITION.
% yx2 = transpose(linspace(0,0.6,N+1)); %FOR TYPE I
% yx2 = transpose(linspace(0,0.55,N)); %FOR TYPE II
% figure(1)
% plot(b12,-yx2,'Linewidth',1,'Color',[1 0 0])
% title('Initial Temp. Profile')
% % xlim([18 38])
% ylabel('Depth (m)')
% xlabel('Temperature (°C)')
% grid on

%THERMAL CONDUCTIVITY VECTOR
for qq = 1:N-1
    if qq <= nair
        tk(qq,1) = ka;
        Lk(qq,1) = dx;
    elseif qq > nair && qq <= nair+nsoil
        tk(qq,1) = ks;
        Lk(qq,1) = dx;
    elseif qq > nsoil && qq <= nair+nsoil+nins
        tk(qq,1) = ki;
        Lk(qq,1) = dx;
    else
        tk(qq,1) = kc;
        Lk(qq,1) = dx;
    end
end

%HARMONIC MEAN OF THERMAL CONDUCTIVITY VECTOR
for ff = 1:N-1
    if ff == 1
        %Top node
        htk_sub(ff,1) = 0;
        htk_sup(ff,1) = dx/((dx/(2*tk(ff,1)))+(dx/(2*tk(ff+1,1))));
    elseif ff < N-1
        %Internal node
        htk_sub(ff,1) = dx/((dx/(2*tk(ff-1,1)))+(dx/(2*tk(ff,1))));
        htk_sup(ff,1) = dx/((dx/(2*tk(ff,1)))+(dx/(2*tk(ff+1,1))));
    else
        %2nd last node
        htk_sub(ff,1) = dx/((dx/(2*tk(ff,1)))+(dx/(2*tk(ff-1,1))));
        htk_sup(ff,1) = 0;
    end
end
%SMOOTHING AT INTERFACES
htk_sub(LocTCIint) = (0.5+0.0284)/2;
htk_sup(LocTCIint) = (0.5+0.0284)/2;


%SUB- and SUPER DIAGNOAL ELEMENTS OF A MATRIX, COEFFICIENTS "a" and "c"
for gg = 1:length(Msup)
    if gg == 1
        a(gg,1) = 0;%(-2)*(htk_sub(gg,1)*dt)/(pCa*dx^2);
        Msub(gg,1) = Msub(gg,1)*a(gg,1);
        
        c(gg,1) = (-2)*(htk_sup(gg,1)*dt)/(pCa*dx^2);
        Msup(gg,1) = Msup(gg,1)*c(gg,1);
    elseif gg > 1 && gg <= nair
        a(gg,1) = -(htk_sub(gg,1)*dt)/(pCa*dx^2);
        Msub(gg,1) = Msub(gg,1)*a(gg,1);
        
        c(gg,1) = -(htk_sup(gg,1)*dt)/(pCa*dx^2);
        Msup(gg,1) = Msup(gg,1)*c(gg,1);
    elseif gg > nair && gg <= nair+nsoil
        Cwat(gg,1) = 1000*4184*0.1*rVWC(1);
        a(gg,1) = -(htk_sub(gg,1)*dt)/((pCs+Cwat(gg,1)/0.1)*dx^2);
        Msub(gg,1) = Msub(gg,1)*a(gg,1);
        
        c(gg,1) = -(htk_sup(gg,1)*dt)/((pCs+Cwat(gg,1)/0.1)*dx^2);
        Msup(gg,1) = Msup(gg,1)*c(gg,1);
    elseif gg > nsoil && gg <= nair+nsoil+nins
        a(gg,1) = -(htk_sub(gg,1)*dt)/(pi*Cpi*dx^2);
        Msub(gg,1) = Msub(gg,1)*a(gg,1);
        
        c(gg,1) = -(htk_sup(gg,1)*dt)/(pi*Cpi*dx^2);
        Msup(gg,1) = Msup(gg,1)*c(gg,1);
    else
        a(gg,1) = -(htk_sub(gg,1)*dt)/(pCc*dx^2);
        Msub(gg,1) = Msub(gg,1)*a(gg,1);
        
        c(gg,1) = -(htk_sub(gg,1)*dt)/(pCc*dx^2);
        Msup(gg,1) = Msup(gg,1)*c(gg,1);
    end
end

%DIAGONAL ELEMENTS OF A MATRIX, COEFFICIENT "b"
for g = 1:length(Mmain)
    if g == 1
%         b(g,1) = (1+((2*(htk_sub(g,1)+htk_sup(g,1))*dt)/(pCa*dx^2)));
%         Mdiag(g,1) = Mmain(g,1)*b(g,1); %%%TYPEI BOUNDARY
                                                                da = (1+((2*htk_sup(g,1)*dt)/(pCa*dx^2)));
                                                                Mdiag(g,1) = Mmain(g,1)*da; %%%TYPEII BOUNDARY
    elseif  g > 1 && g <= nair
        b(g,1) = 1+(((htk_sub(g,1)+htk_sup(g,1))*dt)/(pCa*dx^2));
        Mdiag(g,1) = Mmain(g,1)*b(g,1);
        
    elseif g > nair && g <= (nair+nsoil)
        Cwat(g,1) = 1000*4184*0.1*rVWC(1);
        b(g,1) = 1+(((htk_sub(g,1)+htk_sup(g,1))*dt)/((pCs+Cwat(g,1)/0.1)*dx^2));
        Mdiag(g,1) = Mmain(g,1)*b(g,1);
        
    elseif g > nsoil && g <= (nair+nsoil+nins)
        b(g,1) = 1+(((htk_sub(g,1)+htk_sup(g,1))*dt)/(pi*Cpi*dx^2));
        Mdiag(g,1) = Mmain(g,1)*b(g,1);
        
    elseif g > (nair+nsoil+nins) && g < N-1
        b(g,1) = 1+(((htk_sub(g,1)+htk_sup(g,1))*dt)/(pCc*dx^2));
        Mdiag(g,1) = Mmain(g,1)*b(g,1);
    else
        b(g,1) = 1+((2*(htk_sub(g,1)+htk_sup(g,1))*dt)/(pCc*dx^2));
        Mdiag(g,1) = Mmain(g,1)*b(g,1);
    end
end

%POPULATING A MATRIX pt. 3
A = diag(Mdiag) + diag(Msup(1:N-2,1),1) + diag(Msub(2:N-1,1),-1);
Anew = A(1:N-1,1:N-1);
Anew2=Anew; %for reference
const1 = Msup(1,1);
const2 = Msub(N-1,1);
const3 = 2*dt/(pCa*dx);

%POPULATING B MATRIX pt.1
solnspace = zeros(N-1,Dur);
Bnew = zeros(N-1,1);

for m = 1:width(solnspace)
        %FOR INITIAL TIMESTEP
        if m == 1               
            for k = 1:N-1
                Bnew(k,1) = b12(k,1); 
                Bnewplot = Bnew;
            end
%           %EXTERNAL SURFACE BOUNDARY CONDITION
            Bnew(1,1) = Bnew(1,1)+const3*(nRad(m+1,1)+ ((Tair(m+1,1)-Bnew(1,1))*y(5)*hext(m+1,1)) +...
               (((Tair(m+1,1)-Bnew(1,1))*y(5)*hext(m+1,1))/(Ca*(Patm(m+1,1)/10)*(Bnew(1,1)-Tair(m+1,1))...
               ./(0.622*(2.50-(2.36*10^-3)*Bnew(1,1))*((0.611*exp((17.3*Bnew(1,1))...
               ./((Bnew(1,1)+273.3))))-((RH(m+1,1)./100)*(0.611*exp((17.3*Bnew(1,1))...
               ./((Bnew(1,1)+273.3)))))))))); %%%TYPEII (top)
           
            %INTERNAL SURFACE BOUNDARY CONDITION
            Bnew(N-1,1) = Bnew(N-1,1)-const2*Troom(m+1,1); %%%TYPEI (bottom)

            %FOR TESTING        
%             Bnew(1,1) = Bnew(1,1) - ea*Qsurf; %%%TYPEII (top)
%             Bnew(1,1) = Bnew(1,1)-const1*10; %%%TYPEI (top)
%             Bnew(N-1,1) = Bnew(N-1,1)-const2*Tbot; %%%TYPEI (bottom)
            
            %THOMAS ALGORITHM - TRIDIAGONAL MATRIX SOLVER
            for j = 1:length(Anew)
                    if j == 1        
                        Bnew(j,:) = Bnew(j,j)./Anew(j,j);
                        Anew(j,:) = Anew(j,:)./Anew(j,1);
                    else  
                        %multiplying previous eqn by current row's "a"
                        Bnew(j-1,:) = Bnew(j-1,:)*Anew(j,j-1);
                        Anew(j-1,:) = Anew(j-1,:)*Anew(j,j-1);        
                        %subtracting modified previous row from current row
                        Bnew(j,:) = Bnew(j,:) - Bnew(j-1,:);
                        Anew(j,:) = Anew(j,:) - Anew(j-1,:);
                        %dividing through both sides by coeff. on x_q
                        Bnew(j,:) = Bnew(j,:)./Anew(j,j);
                        Anew(j,:) = Anew(j,:)/Anew(j,j);      
                    end
            end                   
            for p = length(Anew):-1:1
                if p == length(Anew)
                    x(p,1) = Bnew(p,1)/Anew(p,p);
                else
                    x(p,1) = (Bnew(p,1) - Anew(p,p+1)*x(p+1,1))/Anew(p,p);
                end
            end
        %PUTTING RESULTS INTO SOLUTION MATRIX
        solnspace(:,m) = x;
        
           
        else  
        %FOR SUBSEQUENT TIMESTEPS
                    clear tk
                    clear Lk
                    clear htk_sub
                    clear htk_sup
                    clear Mmain
                    clear Msup
                    clear Msub
                    clear Mdiag
                    clear b
        Mmain = transpose(ones(1,N-1));  Msup = transpose(ones(1,N-1));  Msub = transpose(ones(1,N-1));  
            %THERMAL CONDUCTIVITY VECTOR
            for qq = 1:N-1
                if qq <= nair
                    tk(qq,1) = ka;
                    Lk(qq,1) = dx;
                elseif qq > nair && qq <= nair+nsoil
                    tk(qq,1) = ks;
%                     tk(qq,1) = 1.4018*rVWC(m+1,1)+0.285;
                    Lk(qq,1) = dx;
                elseif qq > nsoil && qq <= nair+nsoil+nins
                    tk(qq,1) = ki;
                    Lk(qq,1) = dx;
                else
                    tk(qq,1) = kc;
                    Lk(qq,1) = dx;
                end
            end
            clear qq
            %HARMONIC MEAN OF THERMAL CONDUCTIVITY VECTOR
            for ff = 1:N-1
                if ff == 1
                    %Top node
                    htk_sub(ff,1) = 0;
                    htk_sup(ff,1) = dx/((dx/(2*tk(ff,1)))+(dx/(2*tk(ff+1,1))));
                elseif ff < N-1
                    %Internal node
                    htk_sub(ff,1) = dx/((dx/(2*tk(ff-1,1)))+(dx/(2*tk(ff,1))));
                    htk_sup(ff,1) = dx/((dx/(2*tk(ff,1)))+(dx/(2*tk(ff+1,1))));
                else
                    %2nd last node
                    htk_sub(ff,1) = dx/((dx/(2*tk(ff,1)))+(dx/(2*tk(ff-1,1))));
                    htk_sup(ff,1) = 0;
                end
            end
            %SMOOTHING AT INTERFACES
            htk_sub(LocTCIint) = (0.5+0.0284)/2;
            htk_sup(LocTCIint) = (0.5+0.0284)/2;    
            
            clear ff
            clear Cwat

            %SUB- AND SUPER DIAGNOAL ELEMENTS OF A MATRIX, Coefficients "a" and "c"
            for gg = 1:length(Msup)
                if gg == 1
                    a(gg,1) = 0;%(-2)*(htk_sub(gg,1)*dt)/(pCa*dx^2);
                    Msub(gg,1) = Msub(gg,1)*a(gg,1);
        
                    c(gg,1) = (-2)*(htk_sup(gg,1)*dt)/(pCa*dx^2);
                    Msup(gg,1) = Msup(gg,1)*c(gg,1);
                elseif gg > 1 && gg <= nair
                    a(gg,1) = -(htk_sub(gg,1)*dt)/(pCa*dx^2);
                    Msub(gg,1) = Msub(gg,1)*a(gg,1);

                    c(gg,1) = -(htk_sup(gg,1)*dt)/(pCa*dx^2);
                    Msup(gg,1) = Msup(gg,1)*c(gg,1);
                elseif gg > nair && gg <= nair+nsoil
                    Cwat(gg,1) = 1000*4184*0.1*rVWC(m+1,1);
                    a(gg,1) = -(htk_sub(gg,1)*dt)/((pCs+Cwat(gg,1)/0.1)*dx^2);
                    Msub(gg,1) = Msub(gg,1)*a(gg,1);

                    c(gg,1) = -(htk_sup(gg,1)*dt)/((pCs+Cwat(gg,1)/0.1)*dx^2);
                    Msup(gg,1) = Msup(gg,1)*c(gg,1);
                elseif gg > nsoil && gg <= nair+nsoil+nins
                    a(gg,1) = -(htk_sub(gg,1)*dt)/(pi*Cpi*dx^2);
                    Msub(gg,1) = Msub(gg,1)*a(gg,1);

                    c(gg,1) = -(htk_sup(gg,1)*dt)/(pi*Cpi*dx^2);
                    Msup(gg,1) = Msup(gg,1)*c(gg,1);
                else
                    a(gg,1) = -(htk_sub(gg,1)*dt)/(pCc*dx^2);
                    Msub(gg,1) = Msub(gg,1)*a(gg,1);

                    c(gg,1) = -(htk_sub(gg,1)*dt)/(pCc*dx^2);
                    Msup(gg,1) = Msup(gg,1)*c(gg,1);
                end
            end
            clear gg
            
            %DIAGONAL ELEMENTS OF A MATRIX, COEFFICIENT "b"
            for g = 1:length(Mmain)
                 if g == 1
%                            b(g,1) = (1+((2*(htk_sub(g,1)+htk_sup(g,1))*dt)/(pCa*dx^2)));
%                            Mdiag(g,1) = Mmain(g,1)*b(g,1); %%%TYPEI
                                                                da = (1+((2*htk_sup(g,1)*dt)/(pCa*dx^2)));
                                                                Mdiag(g,1) = Mmain(g,1)*da; %%%TYPEII
                elseif  g > 1 && g <= nair
                    b(g,1) = 1+(((htk_sub(g,1)+htk_sup(g,1))*dt)/(pCa*dx^2));
                    Mdiag(g,1) = Mmain(g,1)*b(g,1);

                elseif g > nair && g <= (nair+nsoil)
                    Cwat(g,1) = 1000*4184*0.1*rVWC(m+1,1);
                    b(g,1) = 1+(((htk_sub(g,1)+htk_sup(g,1))*dt)/((pCs+Cwat(g,1)/0.1)*dx^2));
                    Mdiag(g,1) = Mmain(g,1)*b(g,1);

                elseif g > nsoil && g <= (nair+nsoil+nins)
                    b(g,1) = 1+(((htk_sub(g,1)+htk_sup(g,1))*dt)/(pi*Cpi*dx^2));
                    Mdiag(g,1) = Mmain(g,1)*b(g,1);
                    
                elseif g > (nair+nsoil+nins) && g < N-1
                    b(g,1) = 1+(((htk_sub(g,1)+htk_sup(g,1))*dt)/(pCc*dx^2));
                    Mdiag(g,1) = Mmain(g,1)*b(g,1);
                 else
                    b(g,1) = 1+((2*(htk_sub(g,1)+htk_sup(g,1))*dt)/(pCc*dx^2));
                    Mdiag(g,1) = Mmain(g,1)*b(g,1);
                end
            end
            clear g

            const1 = Msup(1,1);
            const2 = Msub(N-1,1);
            const3 = 2*dt/(pCa*dx);
            Bnew = solnspace(:,m-1);

                %EXTERNAL SURFACE BOUNDARY CONDITION
                Bnew(1,1) = Bnew(1,1) + const3*(nRad(m+1,1)+ ((Tair(m+1,1)-Bnew(1,1))*y(5)*hext(m+1,1)) +...
                (((Tair(m+1,1)-Bnew(1,1))*y(5)*hext(m+1,1))/(Ca*(Patm(m+1,1)/10)*(Bnew(1,1)-Tair(m+1,1))...
                ./(0.622*(2.50-(2.36*10^-3)*Bnew(1,1))*((0.611*exp((17.3*Bnew(1,1))...
               ./((Bnew(1,1)+273.3))))-((RH(m+1,1)./100)*(0.611*exp((17.3*Bnew(1,1))...
                ./((Bnew(1,1)+273.3))))))))));  %%%TYPEII (top)
            
                %INTERNAL SURFACE BOUNDARY CONDITION
                Bnew(N-1,1) = Bnew(N-1,1)-const2*Troom(m+1,1); %%%TYPEI (bottom)
            
             %FOR TESTING
             %Bnew(1,1) = Bnew(1,1) - ea*Qsurf; %%%TYPEII (top)
             %Bnew(1,1) = Bnew(1,1)-const1*Ttop; %%%TYPEI (top)
             %Bnew(N-1,1) = Bnew(N-1,1)-const2*Tbot; %%%TYPEI (bottom)
            
            %POPULATING A MATRIX
            A = diag(Mdiag) + diag(Msup(1:N-2,1),1) + diag(Msub(2:N-1,1),-1);
            Anew = A(1:N-1,1:N-1);
            
            %THOMAS ALGORITHM - TRIDIAGONAL MATRIX SOLVER
            for j = 1:length(Anew)
                    if j == 1        
                        Bnew(j,:) = Bnew(j,j)./Anew(j,j);
                        Anew(j,:) = Anew(j,:)./Anew(j,1);
                    else  
                        %multiplying previous eqn by current row's "a"
                        Bnew(j-1,:) = Bnew(j-1,:)*Anew(j,j-1);
                        Anew(j-1,:) = Anew(j-1,:)*Anew(j,j-1);        
                        %subtracting modified previous row from current row
                        Bnew(j,:) = Bnew(j,:) - Bnew(j-1,:);
                        Anew(j,:) = Anew(j,:) - Anew(j-1,:);
                        %dividing through both sides by coeff. on x_q
                        Bnew(j,:) = Bnew(j,:)./Anew(j,j);
                        Anew(j,:) = Anew(j,:)/Anew(j,j);      
                    end
            end         
            for p = length(Anew):-1:1
                if p == length(Anew)
                    x(p,1) = Bnew(p,1)/Anew(p,p);
                else
                    x(p,1) = (Bnew(p,1) - Anew(p,p+1)*x(p+1,1))/Anew(p,p);
                end
            end
            %PUTTING RESULTS INTO SOLUTION MATRIX
            solnspace(:,m) = x;
        end
    end

%OBJECTIVE FUNCTION - SUM OF SQUARED DIFFERENCES
MesTcant = Tcan(2:length(Tcan));
MesTsoilt = Tsoil(2:length(Tsoil));
MesTmatt = Tmat(2:length(Tmat));
ModTcant = transpose(solnspace(locTcan,:));
ModTsoilt = transpose(solnspace(locTsoil,:));
ModTmatt = transpose(solnspace(locTmat,:));
for pm = 48:length(MesTsoilt)
    sq3(pm,1) = (((ModTsoilt(pm,1)-MesTsoilt(pm,1))^2)...
        +(3*(ModTmatt(pm,1)-MesTmatt(pm,1))^2)...
        +((ModTcant(pm,1)-MesTcant(pm,1))^2));
end
SSD = sum(sq3);
end  

