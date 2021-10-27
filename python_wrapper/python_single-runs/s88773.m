%Code for running the model for s88773

%Reading in input parameters
input_parameter_file="s88773_input.mat";
load(input_parameter_file,'a','b');

%***** Parameters *******
k = @(~,state) b+a*state.u; %W/mK - 
% (previous fit): 0.045*state.u + 3000./sqrt(abs(state.u)), (de koker fit):
% 30+0.0288t, (good linear fit): 120+0.02t
d=12800; %kg/m^3 - determined from Smith et al. 2018 isentrope (previous at 11000) 190 GPa for MgO
d2=12800; %kg/m^3 - density of iron in LiF experiment for 220 GPa
c=500; %J/kg K
t0=2000; %K/counts
w=0.2; % width (um)
mesh_dens = 0.1e-6; % spacing of mesh points (um)

%MgO
km=100; %W/mK used to be 5
dm=5800; %kg/m^3 (4200) from Jin et al. 2010
cm=800; %J/kg K (800)
t0m=2000; %K - old value 1000

tlist = [15:0.2:40]*10^-9;

%**** Creating geometry ****

ratio = 7900/d; %ambient density divided by measured density
inc = 0.1; %amount measuring into window
finc = 0.1; %measurement into iron
lwin=2.0; % length of window (um)

lfe_long=ratio*1.07;
lfe2_long=ratio*2.22;
lfe3_long=ratio*3.07;

sfe_long=(lfe_long-finc)*10^-6; % Measurement spot in iron
smgo_long=(lfe_long+inc)*10^-6; % Measurement spot in mgo

sfe2_long=(lfe2_long-finc)*10^-6; % Measurement spot in iron
smgo2_long=(lfe2_long+inc)*10^-6; % Measurement spot in mgo

sfe3_long=(lfe3_long-finc)*10^-6; % Measurement spot in iron
smgo3_long=(lfe3_long+inc)*10^-6; % Measurement spot in mgo


gap = 0.2;
w2=w+gap;
w3=2*w+gap;
w4=2*w+2*gap;
w5=3*w+2*gap;

w6=3*w+3*gap;
w7=4*w+3*gap;
w8=4*w+4*gap;
w9=5*w+4*gap;
w10=5*w+5*gap;
w11=6*w+5*gap;

w12=6*w+6*gap;
w13=7*w+6*gap;
w14=7*w+7*gap;
w15=8*w+7*gap;
w16=8*w+8*gap;
w17=9*w+8*gap;

w18=9*w+9*gap;
w19=10*w+9*gap;
w20=10*w+10*gap;
w21=11*w+10*gap;
w22=11*w+11*gap;
w23=12*w+11*gap;

r_73 = [3 4 0 lfe_long*10^-6 lfe_long*10^-6 0 0 0 w*10^-6 w*10^-6]';
rmgo_73 = [3 4 lfe_long*10^-6 (lfe_long+lwin)*10^-6 (lfe_long+lwin)*10^-6 lfe_long*10^-6 0 0 w*10^-6 w*10^-6]';

r2_73 = [3 4 0 lfe2_long*10^-6 lfe2_long*10^-6 0 w2*10^-6 w2*10^-6 w3*10^-6 w3*10^-6]';
rmgo2_73 = [3 4 lfe2_long*10^-6 (lfe2_long+lwin)*10^-6 (lfe2_long+lwin)*10^-6 lfe2_long*10^-6 w2*10^-6 w2*10^-6 w3*10^-6 w3*10^-6]';

r3_73 = [3 4 0 lfe3_long*10^-6 lfe3_long*10^-6 0 w4*10^-6 w4*10^-6 w5*10^-6 w5*10^-6]';
rmgo3_73 = [3 4 lfe3_long*10^-6 (lfe3_long+lwin)*10^-6 (lfe3_long+lwin)*10^-6 lfe3_long*10^-6 w4*10^-6 w4*10^-6 w5*10^-6 w5*10^-6]';



thermalmodelT = createpde('thermal','transient');

gdm2 = [r_73 rmgo_73 r2_73 rmgo2_73 r3_73 rmgo3_73];
sf2 = 'r_73+rmgo_73+r2_73+rmgo2_73+r3_73+rmgo3_73';
ns2=char('r_73','rmgo_73','r2_73','rmgo2_73','r3_73','rmgo3_73');
ns2=ns2';
g2 = decsg(gdm2,sf2,ns2);
geometryFromEdges(thermalmodelT,g2);

%iron
thermalProperties(thermalmodelT,'ThermalConductivity',k,...
                                'MassDensity',d,...
                                'SpecificHeat',c,'Face',4);                                                  
thermalProperties(thermalmodelT,'ThermalConductivity',k,...
                                'MassDensity',d,...
                                'SpecificHeat',c,'Face',5);                     
thermalProperties(thermalmodelT,'ThermalConductivity',k,...
                                'MassDensity',d,...
                                'SpecificHeat',c,'Face',6);
                            
%MgO
thermalProperties(thermalmodelT,'ThermalConductivity',km,...
                                'MassDensity',dm,...
                                'SpecificHeat',cm,'Face',1);
thermalProperties(thermalmodelT,'ThermalConductivity',km,...
                                'MassDensity',dm,...
                                'SpecificHeat',cm,'Face',2);
thermalProperties(thermalmodelT,'ThermalConductivity',km,...
                                'MassDensity',dm,...
                                'SpecificHeat',cm,'Face',3)


%***Solve***
square=true;

if square==true
    thermalIC(thermalmodelT,t0);             
    thermalBC(thermalmodelT,'Edge',8,'Temperature',@BC_s88773); %8, 12
    thermalBC(thermalmodelT,'Edge',9,'Temperature',@BC_s88773); %9, 13
    thermalBC(thermalmodelT,'Edge',10,'Temperature',@BC_s88773); %10, 14
    
else
    thermalIC(thermalmodelT,t0);             
    thermalBC(thermalmodelT,'Edge',51,'Temperature',@BC_s88773); %8, 12
    thermalBC(thermalmodelT,'Edge',52,'Temperature',@BC_s88773); %9, 13
    thermalBC(thermalmodelT,'Edge',53,'Temperature',@BC_s88773); %10, 14
    
end

R1 = solve(thermalmodelT,tlist);
T1 = R1.Temperature;

getClosestNode = @(p,x,y) min((p(1,:) - x).^2 + (p(2,:) - y).^2);
[~,nid1_73] = getClosestNode( msh.Nodes, sfe_long, (w*10^-6)/2 );
[~,nid2_73] = getClosestNode( msh.Nodes, sfe2_long, (1.5*w+gap)*10^-6);
[~,nid3_73] = getClosestNode( msh.Nodes, sfe3_long, (2.5*w+2*gap)*10^-6);

%******    Plotting    ******

%SOP data for s88773%
start1=179;
start2=259;
start3=310;
stop1=265;
stop2=415;
stop3=706;
stop_final=750;
opts = detectImportOptions('../s88773_sop_lineouts_TP.xlsx','Sheet','June2021');
m = readtable('../s88773_sop_lineouts_TP.xlsx',opts,'Sheet','June2021');
t_data1=m.Time(start1:stop1)*10^-9;
t_data2=m.Time(start2:stop2)*10^-9;
t_data3=m.Time(start3:stop3)*10^-9;
t_data1end=m.Time(stop1:stop_final)*10^-9;
t_data2end=m.Time(stop2:stop_final)*10^-9;
t_data2end=m.Time(stop3:stop_final)*10^-9;
temp1_data=m.step1_corrected(start1:stop1);
temp2_data=m.step2_corrected(start2:stop2);
temp3_data=m.step3_corrected(start3:stop3);
xw = 1; %throughput correction based on slit width and shot number
aeta=25.5; %1/delta_t of SOP data, was 46
a0=481000;% Based on ND filter of 0 (481000) 
t0a=1.909;% Based on ND filter of 0 (1.909) 

reference = 0.55; %calculatd based on reflectivity data of ambient iron assuming perfect window transparency

temp1=real(11605*t0a./(log(1+((1-reference)*a0./(aeta*temp1_data)))));
temp2=real(11605*t0a./(log(1+((1-reference)*a0./(aeta*temp2_data)))));
temp3=real(11605*t0a./(log(1+((1-reference)*a0./(aeta*temp3_data)))));

figure
ax1=axes('Position',[0.13,0.11,0.85,0.69]);
nexttile
plot(t_data1*10^9,temp1*10^-3,'o',t_data2*10^9,temp2*10^-3,'o',t_data3*10^9,temp3*10^-3,'o',...
    tlist*10^9,T1(nid1_73,:)*10^-3,'-',tlist*10^9,T1(nid2_73,:)*10^-3,'-',tlist*10^9,T1(nid3_73,:)*10^-3,'-','Linewidth',2)
legend('Data 1.07 um','Data 2.3 um','Data 3.23 um','Model 1.07 um','Model 2.3 um', 'Model 3.23 um')
grid on
if square==true
    title '1D Model With square Input: s88773 (Fe-LiF)';
else
    title '1D Model With exp Input: s88773 (Fe-LiF)';
end
xlabel 'Time (nanoseconds)'
ylabel 'Temperature (K)'
%xticklabels({'15','20','25','30','35','40','45','50','55'})
%ylim([0 3e4]);

T11=T1(nid1_73,:);
T12=T1(nid2_73,:);
T13=T1(nid3_73,:);

save('73_output.mat','tlist','T11','T12','T13','t_data1','t_data2','t_data3','temp1','temp2','temp3')