function SSD = RC_Model(x)

%Boundaries for "str" and "fin" for monthly simulations
%     may = [2905:3648];
%     jun = [3649:4368];
%     jul = [4369:5112];
%     aug = [5113:5856];
%     sep = [5857:6576];

%Set temporal boundaries for calibration (see above)
%Here, only change the 4-digit number in line 12, i.e. only change "2905"
%and keep the "-1" in place. 
str = 2905-1; 
fin = 3648; 

%READING IN HOURLY DATA
% [num,~,~] = xlsread('hourly2016.xlsx');

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

%EXTERIOR CONVECTION COEFFICIENT
for hn = 2:length(Tsoil)
    Rs(hn,1) = 0.1/(1.4018*rVWC(hn,1)+0.285);
    hext(hn,1) = (5.9+4.1*wind(hn,1)*((511+294)/(511+Tair(hn,1)+273)));
end

Cins = 45*1300*0.1;

    for i = 2:length(Tsoil)
        %INITIAL TIMESTEP
        if i == 2
            Tcan0 = Tcan(i);
            ModTcant(i,1) = Tcan0;
            MesTcant(i,1) = Tcan0;
            
            Tsoil0 = Tsoil(i);
            ModTsoilt(i,1) = Tsoil0;
            MesTsoilt(i,1) = Tsoil0;
            
            Troom0 = Troom(i);
            ModTroomt(i,1) = Troom0;
            MesTroomt(i,1) = Troom0;
            
            Tmat0 = Tmat(i);
            ModTmatt(i,1) = Tmat0;
            MesTmatt(i,1) = Tmat0;
            
            TConc0 = TConc(i,1);
            MesTConct(i,1) = TConc0;
            ModTConct(i,1) = TConc0;
            
            Tins0 = (Tmat0+TConc0)/2;                    
            MesTinst(i,1) = Tins0;
            ModTinst(i,1) = Tins0;
        else 

            %TIMESTEPS AFTER 1 HOUR:
            %Canopy Node
            dModTcan(i,1) = (((Tair(i-1,1)-ModTcant(i-1,1))*(x(5)*hext(i-1,1))) + ...
                ((ModTsoilt(i-1,1)-ModTcant(i-1,1))/((x(4)*0.5+x(1)))) +...
                nRad(i-1,1) + (((Tair(i-1,1)-ModTcant(i-1,1))*(x(5)*hext(i-1,1)))/(Ca*(Patm(i-1,1)/10)*(ModTcant(i-1,1)-Tair(i-1,1))./(0.622*(2.50-(2.36*10^-3)*ModTcant(i-1,1))*((0.611*exp((17.3*ModTcant(i-1,1))./((ModTcant(i-1,1)+273.3))))-((RH(i-1,1)./100)*(0.611*exp((17.3*ModTcant(i-1,1))./((ModTcant(i-1,1)+273.3))))))))) )...
                *(3600/(x(2)*582*4800*0.1+(1-x(2))*1.184*1007*0.1));        

            %Soil Node
             dModTsoil(i,1) = (((ModTcant(i-1,1)-ModTsoilt(i-1,1))/(x(4)*0.5+x(1))) +...
                ((ModTmatt(i-1,1)-ModTsoilt(i-1,1))/(x(4)*0.5)))*(3600/(x(3)+1000*4184*rVWC(i-1,1)*0.1));                                
            
            %Mat Node
            ModTmatt(i,1) = 2*((ModTsoilt(i-1,1)/x(4))+(ModTinst(i-1,1)/x(6)))*((2/x(4)+(2/x(6)))^-1);

            %Ins Node
            dModTins(i,1) = (((ModTmatt(i-1,1)-ModTinst(i-1,1))/(x(6)*0.5)) +...
                ((ModTConct(i-1,1)-ModTinst(i-1,1))/(x(6)*0.5+x(7)*0.5)))*(3600/(Cins));
           
            %Conc node
            dModTConc(i,1) = (((ModTinst(i-1,1)-ModTConct(i-1,1))/(x(6)*0.5+x(7)*0.5)) +...
                ((MesTroomt(i-1,1)-ModTConct(i-1,1))/(x(7)*0.5)))*(3600/(x(8)));


              
 
            %For plotting purposes         
            ModTcant(i,1) = ModTcant(i-1,1) + dModTcan(i,1);
            dMesTcan(i,1) = Tcan(i,1)-Tcan(i-1,1);
            MesTcant(i,1) = MesTcant(i-1,1) + dMesTcan(i,1);
            
            ModTsoilt(i,1) = ModTsoilt(i-1,1) + dModTsoil(i,1);
            dMesTsoil(i,1) = Tsoil(i,1)-Tsoil(i-1,1);
            MesTsoilt(i,1) = MesTsoilt(i-1,1) + dMesTsoil(i,1);
            
            dMesTmat(i,1) = Tmat(i,1)-Tmat(i-1,1);
            MesTmatt(i,1) = MesTmatt(i-1,1) + dMesTmat(i,1);
            
            ModTinst(i,1) = ModTinst(i-1,1) + dModTins(i,1);
        
            ModTConct(i,1) = ModTConct(i-1,1) + dModTConc(i,1);
            
            dMesTroom(i,1) = Troom(i,1)-Troom(i-1,1);
            MesTroomt(i,1) = MesTroomt(i-1,1) + dMesTroom(i,1);

        end
    end
    
    


%OBJECTIVE FUNCTION, Sum of Squared Differences
sq = zeros(length(MesTcant),1);

for pm = 2:length(MesTcant)     
    sq(pm,1) = ((ModTcant(pm,1)-MesTcant(pm,1))^2)...
        + ((ModTsoilt(pm,1)-MesTsoilt(pm,1))^2) +...
        (3*(ModTmatt(pm,1)-MesTmatt(pm,1))^2);
end

SSD = sum(sq);
end

