%Main script for analyzing multiple shots simultaneously.
%s88773 MgO, s88776 LiF, s88780 MgO short


%****     2um     *****
input_parameter_file="inputs/1d_combined_input_matrix.mat";
load(input_parameter_file,'a','b');

%***** Parameters *******
k = @(~,state) b+a*state.u; %W/mK - 
% (previous fit): 0.045*state.u + 3000./sqrt(abs(state.u)), (de koker fit):
% 30+0.0288t, (good linear fit): 120+0.02t
d=12500; %kg/m^3 - determined from Smith et al. 2018 isentrope (previous at 11000)
c=500; %J/kg K
t0=2000; %K/counts
w=0.2; % width (um)
mesh_dens = 0.1e-6; % spacing of mesh points (um)

%MgO
km=100; %W/mK used to be 5
dm=5800; %kg/m^3 (4200) from Jin et al. 2010
cm=800; %J/kg K (800)
t0m=1000; %K - old value 1000

%LiF
klif=50; %W/mK used to be 5
dlif=4700; %kg/m^3 
clif=1600; %J/kg K (800)
t0lif=1000; %K - old value 1000

tlist = [20:0.5:50]*10^-9;


%**** Creating geometry ****

ratio = 7900/d; %ambient density divided by measured density
inc = 0.1; %amount measuring into window
finc = 0.1; %measurement into iron
lwin=2.0; % length of window (um)


lfe_long=ratio*1.07; % Length of iron (um) 0.89 0.72 0.77
lfe2_long=ratio*2.3; % Length of iron (um) 1.57
lfe3_long=ratio*3.23; % Length of iron (um) 1.57

lfe_short=ratio*0.45; % Length of iron (um) 0.89 0.72 0.77
lfe2_short=ratio*1.68; % Length of iron (um) 1.57
lfe3_short=ratio*2.61; % Length of iron (um) 1.57

sfe_long=(lfe_long-finc)*10^-6; % Measurement spot in iron
smgo_long=(lfe_long+inc)*10^-6; % Measurement spot in mgo

sfe2_long=(lfe2_long-finc)*10^-6; % Measurement spot in iron
smgo2_long=(lfe2_long+inc)*10^-6; % Measurement spot in mgo

sfe3_long=(lfe3_long-finc)*10^-6; % Measurement spot in iron
smgo3_long=(lfe3_long+inc)*10^-6; % Measurement spot in mgo


sfe_short=(lfe_short-finc)*10^-6; % Measurement spot in iron
smgo_short=(lfe_short+inc)*10^-6; % Measurement spot in mgo

sfe2_short=(lfe2_short-finc)*10^-6; % Measurement spot in iron
smgo2_short=(lfe2_short+inc)*10^-6; % Measurement spot in mgo

sfe3_short=(lfe3_short-finc)*10^-6; % Measurement spot in iron
smgo3_short=(lfe3_short+inc)*10^-6; % Measurement spot in mgo


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

r_73 = [3 4 0 lfe_long*10^-6 lfe_long*10^-6 0 0 0 w*10^-6 w*10^-6]';
rmgo_73 = [3 4 lfe_long*10^-6 (lfe_long+lwin)*10^-6 (lfe_long+lwin)*10^-6 lfe_long*10^-6 0 0 w*10^-6 w*10^-6]';

r2_73 = [3 4 0 lfe2_long*10^-6 lfe2_long*10^-6 0 w2*10^-6 w2*10^-6 w3*10^-6 w3*10^-6]';
rmgo2_73 = [3 4 lfe2_long*10^-6 (lfe2_long+lwin)*10^-6 (lfe2_long+lwin)*10^-6 lfe2_long*10^-6 w2*10^-6 w2*10^-6 w3*10^-6 w3*10^-6]';

r3_73 = [3 4 0 lfe3_long*10^-6 lfe3_long*10^-6 0 w4*10^-6 w4*10^-6 w5*10^-6 w5*10^-6]';
rmgo3_73 = [3 4 lfe3_long*10^-6 (lfe3_long+lwin)*10^-6 (lfe3_long+lwin)*10^-6 lfe3_long*10^-6 w4*10^-6 w4*10^-6 w5*10^-6 w5*10^-6]';


r_76 = [3 4 0 lfe_long*10^-6 lfe_long*10^-6 0 w6*10^-6 w6*10^-6 w7*10^-6 w7*10^-6]';
rlif = [3 4 lfe_long*10^-6 (lfe_long+lwin)*10^-6 (lfe_long+lwin)*10^-6 lfe_long*10^-6 w6*10^-6 w6*10^-6 w7*10^-6 w7*10^-6]';

r2_76 = [3 4 0 lfe2_long*10^-6 lfe2_long*10^-6 0 w8*10^-6 w8*10^-6 w9*10^-6 w9*10^-6]';
rlif2 = [3 4 lfe2_long*10^-6 (lfe2_long+lwin)*10^-6 (lfe2_long+lwin)*10^-6 lfe2_long*10^-6 w8*10^-6 w8*10^-6 w9*10^-6 w9*10^-6]';

r3_76 = [3 4 0 lfe3_long*10^-6 lfe3_long*10^-6 0 w10*10^-6 w10*10^-6 w11*10^-6 w11*10^-6]';
rlif3 = [3 4 lfe3_long*10^-6 (lfe3_long+lwin)*10^-6 (lfe3_long+lwin)*10^-6 lfe3_long*10^-6 w10*10^-6 w10*10^-6 w11*10^-6 w11*10^-6]';


r_80 = [3 4 0 lfe_short*10^-6 lfe_short*10^-6 0 w12*10^-6 w12*10^-6 w13*10^-6 w13*10^-6]';
rmgo_80 = [3 4 lfe_short*10^-6 (lfe_short+lwin)*10^-6 (lfe_short+lwin)*10^-6 lfe_short*10^-6 w12*10^-6 w12*10^-6 w13*10^-6 w13*10^-6]';

r2_80 = [3 4 0 lfe2_short*10^-6 lfe2_short*10^-6 0 w14*10^-6 w14*10^-6 w15*10^-6 w15*10^-6]';
rmgo2_80 = [3 4 lfe2_short*10^-6 (lfe2_short+lwin)*10^-6 (lfe2_short+lwin)*10^-6 lfe2_short*10^-6 w14*10^-6 w14*10^-6 w15*10^-6 w15*10^-6]';

r3_80 = [3 4 0 lfe3_short*10^-6 lfe3_short*10^-6 0 w16*10^-6 w16*10^-6 w17*10^-6 w17*10^-6]';
rmgo3_80 = [3 4 lfe3_short*10^-6 (lfe3_short+lwin)*10^-6 (lfe3_short+lwin)*10^-6 lfe3_short*10^-6 w16*10^-6 w16*10^-6 w17*10^-6 w17*10^-6]';



thermalmodelT = createpde('thermal','transient');

gdm2 = [r_73 rmgo_73 r2_73 rmgo2_73 r3_73 rmgo3_73 r_76 rlif r2_76 rlif2 r3_76 rlif3 r_80 rmgo_80 r2_80 rmgo2_80 r3_80 rmgo3_80];
sf2 = 'r_73+rmgo_73+r2_73+rmgo2_73+r3_73+rmgo3_73+r_76+rlif+r2_76+rlif2+r3_76+rlif3+r_80+rmgo_80+r2_80+rmgo2_80+r3_80+rmgo3_80';
ns2=char('r_73','rmgo_73','r2_73','rmgo2_73','r3_73','rmgo3_73', 'r_76', 'rlif', 'r2_76', 'rlif2', 'r3_76', 'rlif3', 'r_80', 'rmgo_80', 'r2_80', 'rmgo2_80', 'r3_80', 'rmgo3_80');
ns2=ns2';
g2 = decsg(gdm2,sf2,ns2);
geometryFromEdges(thermalmodelT,g2);

%iron
thermalProperties(thermalmodelT,'ThermalConductivity',k,...
                                'MassDensity',d,...
                                'SpecificHeat',c,'Face',4);                                                  
thermalProperties(thermalmodelT,'ThermalConductivity',k,...
                                'MassDensity',d,...
                                'SpecificHeat',c,'Face',6);                     
thermalProperties(thermalmodelT,'ThermalConductivity',k,...
                                'MassDensity',d,...
                                'SpecificHeat',c,'Face',8);

thermalProperties(thermalmodelT,'ThermalConductivity',k,...
                                'MassDensity',d,...
                                'SpecificHeat',c,'Face',10);                                                  
thermalProperties(thermalmodelT,'ThermalConductivity',k,...
                                'MassDensity',d,...
                                'SpecificHeat',c,'Face',12);                     
thermalProperties(thermalmodelT,'ThermalConductivity',k,...
                                'MassDensity',d,...
                                'SpecificHeat',c,'Face',14);
                            
thermalProperties(thermalmodelT,'ThermalConductivity',k,...
                                'MassDensity',d,...
                                'SpecificHeat',c,'Face',15);                                                  
thermalProperties(thermalmodelT,'ThermalConductivity',k,...
                                'MassDensity',d,...
                                'SpecificHeat',c,'Face',16);                     
thermalProperties(thermalmodelT,'ThermalConductivity',k,...
                                'MassDensity',d,...
                                'SpecificHeat',c,'Face',17);                            
                            
%MgO
thermalProperties(thermalmodelT,'ThermalConductivity',km,...
                                'MassDensity',dm,...
                                'SpecificHeat',cm,'Face',1);
thermalProperties(thermalmodelT,'ThermalConductivity',km,...
                                'MassDensity',dm,...
                                'SpecificHeat',cm,'Face',2);
thermalProperties(thermalmodelT,'ThermalConductivity',km,...
                                'MassDensity',dm,...
                                'SpecificHeat',cm,'Face',3);
                            
thermalProperties(thermalmodelT,'ThermalConductivity',km,...
                                'MassDensity',dm,...
                                'SpecificHeat',cm,'Face',5);
thermalProperties(thermalmodelT,'ThermalConductivity',km,...
                                'MassDensity',dm,...
                                'SpecificHeat',cm,'Face',7);
thermalProperties(thermalmodelT,'ThermalConductivity',km,...
                                'MassDensity',dm,...
                                'SpecificHeat',cm,'Face',9);
                            
                            
%LiF
thermalProperties(thermalmodelT,'ThermalConductivity',klif,...
                                'MassDensity',dlif,...
                                'SpecificHeat',clif,'Face',11);
thermalProperties(thermalmodelT,'ThermalConductivity',klif,...
                                'MassDensity',dlif,...
                                'SpecificHeat',clif,'Face',13);
thermalProperties(thermalmodelT,'ThermalConductivity',klif,...
                                'MassDensity',dlif,...
                                'SpecificHeat',clif,'Face',18);

                          
msh = generateMesh(thermalmodelT,'Hmax',mesh_dens);

% figure 
% pdeplot(thermalmodelT); %
% axis equal
% title 'Block With Finite Element Mesh Displayed'
% figure
% pdegplot(thermalmodelT,'EdgeLabels','on','FaceLabels','on')


%***Solve***


thermalIC(thermalmodelT,t0);             
thermalBC(thermalmodelT,'Edge',43,'Temperature',@BC_external_exp); %8, 12
thermalBC(thermalmodelT,'Edge',44,'Temperature',@BC_external_exp); %9, 13
thermalBC(thermalmodelT,'Edge',45,'Temperature',@BC_external_exp); %10, 14

thermalBC(thermalmodelT,'Edge',46,'Temperature',@BC_external_exp); %8, 12
thermalBC(thermalmodelT,'Edge',47,'Temperature',@BC_external_exp); %9, 13
thermalBC(thermalmodelT,'Edge',48,'Temperature',@BC_external_exp); %10, 14

thermalBC(thermalmodelT,'Edge',49,'Temperature',@BC_external_exp); %8, 12
thermalBC(thermalmodelT,'Edge',50,'Temperature',@BC_external_exp); %9, 13
thermalBC(thermalmodelT,'Edge',51,'Temperature',@BC_external_exp); %10, 14



R1 = solve(thermalmodelT,tlist);
T1 = R1.Temperature;

getClosestNode = @(p,x,y) min((p(1,:) - x).^2 + (p(2,:) - y).^2);
[~,nid1_73] = getClosestNode( msh.Nodes, sfe_long, (w*10^-6)/2 );
[~,nid2_73] = getClosestNode( msh.Nodes, sfe2_long, (1.5*w+gap)*10^-6);
[~,nid3_73] = getClosestNode( msh.Nodes, sfe3_long, (2.5*w+2*gap)*10^-6);

[~,nid1_76] = getClosestNode( msh.Nodes, sfe_long, (3.5*w+3*gap)*10^-6 );
[~,nid2_76] = getClosestNode( msh.Nodes, sfe2_long, (4.5*w+4*gap)*10^-6);
[~,nid3_76] = getClosestNode( msh.Nodes, sfe3_long, (5.5*w+5*gap)*10^-6);

[~,nid1_80] = getClosestNode( msh.Nodes, sfe_short, (6.5*w+6*gap)*10^-6 );
[~,nid2_80] = getClosestNode( msh.Nodes, sfe2_short, (7.5*w+7*gap)*10^-6);
[~,nid3_80] = getClosestNode( msh.Nodes, sfe3_short, (8.5*w+8*gap)*10^-6);





%[~,nidfe1] = getClosestNode( msh.Nodes, sfe, (3*w*10^-6)/2);
% [~,nid1m] = getClosestNode( msh.Nodes, smgo, (w*10^-6)/2 );
% [~,nid2m] = getClosestNode( msh.Nodes, smgo2, (1.5*w+gap)*10^-6);
% [~,nid3m] = getClosestNode( msh.Nodes, smgo3, (2.5*w+2*gap)*10^-6);


%******    Plotting    ******

%SOP data for s88773%
start=50;
stop=950;
start1_73=150;
start2_73=195;
start3_73=220;
stop1_73=300;
stop2_73=950;
stop3_73=950;

opts = detectImportOptions('s88773_sop_lineouts_TP.xlsx','Sheet','June2021');
m = readtable('s88773_sop_lineouts_TP.xlsx',opts,'Sheet','June2021');
blank_t_data_73=m.Time*10^-9;
t_data_73=m.Time(start:stop)*10^-9;

temp1_data73=m.step1_corrected(start:stop);
temp2_data73=m.step2_corrected(start:stop);
temp3_data73=m.step3_corrected(start:stop);

% temp1_data73=m.step1_corrected(start1_73:stop1_73);
% temp2_data73=m.step2_corrected(start2_73:stop2_73);
% temp3_data73=m.step3_corrected(start3_73:stop3_73);

xw = 1; %throughput correction based on slit width and shot number
aeta=25.5; %1/delta_t of SOP data, was 46
a0=481000;% Based on ND filter of 0 (481000) (395900)
t0a=1.909;% Based on ND filter of 0 (1.909) (1.91)

%SOP data for s88776
start=50;
stop=750;
start1_76=50;
start2_76=50;
start3_76=50;
stop1_76=950;
stop2_76=950;
stop3_76=950;
opts = detectImportOptions('s88776_SOP_TP.xlsx');
m = readtable('s88776_SOP_TP.xlsx',opts);
blank_t_data_76=m.time*10^-9;
t_data_76=m.time(start:stop)*10^-9;
temp1_data76=m.step1_corrected(start:stop);
temp2_data76=m.step2_corrected(start:stop);
temp3_data76=m.step3_corrected(start:stop);
xw = 1; %throughput correction based on slit width and shot number
aeta=25.5; %1/delta_t of SOP data, was 46
a0_76=395900;% Based on ND filter of 0.1 (481000) (395900)
t0a_76=1.91;% Based on ND filter of 0.1 (1.909) (1.91)

%SOP data for s88780
start=50;
stop=750;
opts = detectImportOptions('s88780_SOP_TP.xlsx');
m = readtable('s88780_SOP_TP.xlsx',opts);
blank_t_data_80=m.time*10^-9;
t_data_80=m.time(start:stop)*10^-9;
temp1_data80=m.step1_corrected(start:stop);
temp2_data80=m.step2_corrected(start:stop);
temp3_data80=m.step3_corrected(start:stop);
xw = 1; %throughput correction based on slit width and shot number
aeta=25.5; %1/delta_t of SOP data, was 46
a0=481000;% Based on ND filter of 0 (481000) (395900)
t0a=1.909;% Based on ND filter of 0 (1.909) (1.91)

reference = 0.55; %calculatd based on reflectivity data of ambient iron assuming perfect window transparency

temp1_73=real(11605*t0a./(log(1+((1-reference)*a0./(aeta*temp1_data73)))));
temp2_73=real(11605*t0a./(log(1+((1-reference)*a0./(aeta*temp2_data73)))));
temp3_73=real(11605*t0a./(log(1+((1-reference)*a0./(aeta*temp3_data73)))));

temp1_76=real(11605*t0a_76./(log(1+((1-reference)*a0_76./(aeta*temp1_data76)))));
temp2_76=real(11605*t0a_76./(log(1+((1-reference)*a0_76./(aeta*temp2_data76)))));
temp3_76=real(11605*t0a_76./(log(1+((1-reference)*a0_76./(aeta*temp3_data76)))));

temp1_80=real(11605*t0a./(log(1+((1-reference)*a0./(aeta*temp1_data80)))));
temp2_80=real(11605*t0a./(log(1+((1-reference)*a0./(aeta*temp2_data80)))));
temp3_80=real(11605*t0a./(log(1+((1-reference)*a0./(aeta*temp3_data80)))));

figure
ax1=axes('Position',[0.13,0.11,0.85,0.69]);
nexttile
plot(t_data_73,temp1_73,'o',t_data_73,temp2_73,'o',t_data_73,temp3_73,'o',...
    tlist,T1(nid1_73,:),'-',tlist,T1(nid2_73,:),'-',tlist,T1(nid3_73,:),'-','Linewidth',2)
legend('Data 1 um','Data 2 um','Data 3 um','Fe 1um','Fe 2um', 'Fe 3um')
grid on
title '1D Model With Exp Input: s88773';
xlabel 'Time (nanoseconds)'
ylabel 'Temperature (K)'
%xticklabels({'15','20','25','30','35','40','45','50','55'})
ylim([0 3e4]);

figure
ax1=axes('Position',[0.13,0.11,0.85,0.69]);
plot(t_data_76,temp1_76,'o',t_data_76,temp2_76,'o',t_data_76,temp3_76,'o',...
    tlist,T1(nid1_76,:),'-',tlist,T1(nid2_76,:),'-',tlist,T1(nid3_76,:),'-','Linewidth',2)
legend('Data 1 um','Data 2 um','Data 3 um','Fe 1um','Fe 2um', 'Fe 3um')
grid on
title '1D Model With Exp Input: s88776';
xlabel 'Time (nanoseconds)'
ylabel 'Temperature (K)'
%xticklabels({'15','20','25','30','35','40','45','50','55'})
ylim([0 3e4]);

figure
ax1=axes('Position',[0.13,0.11,0.85,0.69]);
plot(t_data_80,temp1_80,'o',t_data_80,temp2_80,'o',t_data_80,temp3_80,'o',...
    tlist,T1(nid1_80,:),'-',tlist,T1(nid2_80,:),'-',tlist,T1(nid3_80,:),'-','Linewidth',2)
legend('Data 1 um','Data 2 um','Data 3 um','Fe 1um','Fe 2um', 'Fe 3um')
grid on
title '1D Model With Exp Input: s88780';
xlabel 'Time (nanoseconds)'
ylabel 'Temperature (K)'
%xticklabels({'15','20','25','30','35','40','45','50','55'})
ylim([0 3e4]);

T11=T1(nid1_73,:);
T12=T1(nid2_73,:);
T13=T1(nid3_73,:);
T21=T1(nid1_76,:);
T22=T1(nid2_76,:);
T23=T1(nid3_76,:);
T31=T1(nid1_80,:);
T32=T1(nid2_80,:);
T33=T1(nid3_80,:);

save('FEM_output/1d_combo_run_output.mat','tlist','T11','T12','T13','T21','T22','T23','T31','T32','T33','t_data_73','temp1_73','temp2_73','temp3_73','t_data_76','temp1_76','temp2_76','temp3_76','t_data_80','temp1_80','temp2_80','temp3_80')
