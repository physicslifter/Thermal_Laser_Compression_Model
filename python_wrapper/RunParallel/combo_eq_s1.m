
%s88773 MgO, s88776 LiF, s88780 MgO short


%****     2um     *****

input_parameter_file="1d_combined_input_matrixs1.mat";
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

%LiF
klif=50; %W/mK used to be 5
dlif=4700; %kg/m^3 
clif=1600; %J/kg K (800)
t0lif=2000; %K - old value 1000

tlist = [15:0.2:40]*10^-9;


%**** Creating geometry ****

ratio = 7900/d; %ambient density divided by measured density
inc = 0.1; %amount measuring into window
finc = 0.1; %measurement into iron
lwin=2.0; % length of window (um)

%Old values, lfe_long for 73 & 80, lfe_short for s88776
%lfe_long=ratio*1.07; % Length of iron (um) 0.89 0.72 0.77
%lfe2_long=ratio*2.3; % Length of iron (um) 1.57
%lfe3_long=ratio*3.23; % Length of iron (um) 1.57

%lfe_short=ratio*0.45; % Length of iron (um) 0.89 0.72 0.77
%lfe2_short=ratio*1.68; % Length of iron (um) 1.57
%lfe3_short=ratio*2.61; % Length of iron (um) 1.57

%lfe_83=ratio*0.76; % Length of iron (um) 0.89 0.72 0.77
%lfe2_83=ratio*2.76; % Length of iron (um) 1.57
%lfe3_83=ratio*4.76; % Length of iron (um) 1.57

%New values
%s88773 & s888776
lfe_long=ratio*1.07;
lfe2_long=ratio*2.22;
lfe3_long=ratio*3.07;

%s88780
lfe_short=ratio*0.45;
lfe2_short=ratio*1.7;
lfe3_short=ratio*2.8;

lfe_83=ratio*0.76;
lfe2_83=ratio*2.76;
lfe3_83=ratio*4.76;


sfe_long=(lfe_long-finc)*10^-6; % Measurement spot in iron
smgo_long=(lfe_long+inc)*10^-6; % Measurement spot in mgo

sfe2_long=(lfe2_long-finc)*10^-6; % Measurement spot in iron
smgo2_long=(lfe2_long+inc)*10^-6; % Measurement spot in mgo

sfe3_long=(lfe3_long-finc)*10^-6; % Measurement spot in iron
smgo3_long=(lfe3_long+inc)*10^-6; % Measurement spot in mgo


sfe_short=(lfe_short-finc)*10^-6; % Measurement spot in iron
smgo1_short=(lfe_short+inc)*10^-6; % Measurement spot in mgo

sfe2_short=(lfe2_short-finc)*10^-6; % Measurement spot in iron
smgo2_short=(lfe2_short+inc)*10^-6; % Measurement spot in mgo

sfe3_short=(lfe3_short-finc)*10^-6; % Measurement spot in iron
smgo3_short=(lfe3_short+inc)*10^-6; % Measurement spot in mgo


sfe_83=(lfe_83-finc)*10^-6; % Measurement spot in iron
smgo1_83=(lfe_83+inc)*10^-6; % Measurement spot in mgo

sfe2_83=(lfe2_83-finc)*10^-6; % Measurement spot in iron
smgo2_83=(lfe2_83+inc)*10^-6; % Measurement spot in mgo

sfe3_83=(lfe3_83-finc)*10^-6; % Measurement spot in iron
smgo3_83=(lfe3_83+inc)*10^-6; % Measurement spot in mgo


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


r_83 = [3 4 0 lfe_83*10^-6 lfe_83*10^-6 0 w18*10^-6 w18*10^-6 w19*10^-6 w19*10^-6]';
rmgo_83 = [3 4 lfe_83*10^-6 (lfe_83+lwin)*10^-6 (lfe_83+lwin)*10^-6 lfe_83*10^-6 w18*10^-6 w18*10^-6 w19*10^-6 w19*10^-6]';

r2_83 = [3 4 0 lfe2_83*10^-6 lfe2_83*10^-6 0 w20*10^-6 w20*10^-6 w21*10^-6 w21*10^-6]';
rmgo2_83 = [3 4 lfe2_83*10^-6 (lfe2_83+lwin)*10^-6 (lfe2_83+lwin)*10^-6 lfe2_83*10^-6 w20*10^-6 w20*10^-6 w21*10^-6 w21*10^-6]';

r3_83 = [3 4 0 lfe3_83*10^-6 lfe3_83*10^-6 0 w22*10^-6 w22*10^-6 w23*10^-6 w23*10^-6]';
rmgo3_83 = [3 4 lfe3_83*10^-6 (lfe3_83+lwin)*10^-6 (lfe3_83+lwin)*10^-6 lfe3_83*10^-6 w22*10^-6 w22*10^-6 w23*10^-6 w23*10^-6]';



thermalmodelT = createpde('thermal','transient');

gdm2 = [r_73 rmgo_73 r2_73 rmgo2_73 r3_73 rmgo3_73 r_76 rlif r2_76 rlif2 r3_76 rlif3 r_80 rmgo_80 r2_80 rmgo2_80 r3_80 rmgo3_80 r_83 rmgo_83 r2_83 rmgo2_83 r3_83 rmgo3_83];
sf2 = 'r_73+rmgo_73+r2_73+rmgo2_73+r3_73+rmgo3_73+r_76+rlif+r2_76+rlif2+r3_76+rlif3+r_80+rmgo_80+r2_80+rmgo2_80+r3_80+rmgo3_80+r_83+rmgo_83+r2_83+rmgo2_83+r3_83+rmgo3_83';
ns2=char('r_73','rmgo_73','r2_73','rmgo2_73','r3_73','rmgo3_73', 'r_76', ...
    'rlif', 'r2_76', 'rlif2', 'r3_76', 'rlif3', 'r_80', 'rmgo_80', 'r2_80',...
    'rmgo2_80', 'r3_80', 'rmgo3_80',  'r_83', 'rmgo_83', 'r2_83',...
    'rmgo2_83', 'r3_83', 'rmgo3_83');
ns2=ns2';
g2 = decsg(gdm2,sf2,ns2);
geometryFromEdges(thermalmodelT,g2);

%iron
thermalProperties(thermalmodelT,'ThermalConductivity',k,...
                                'MassDensity',d,...
                                'SpecificHeat',c,'Face',7);                                                  
thermalProperties(thermalmodelT,'ThermalConductivity',k,...
                                'MassDensity',d,...
                                'SpecificHeat',c,'Face',9);                     
thermalProperties(thermalmodelT,'ThermalConductivity',k,...
                                'MassDensity',d,...
                                'SpecificHeat',c,'Face',11);

thermalProperties(thermalmodelT,'ThermalConductivity',k,...
                                'MassDensity',d,...
                                'SpecificHeat',c,'Face',18);                                                  
thermalProperties(thermalmodelT,'ThermalConductivity',k,...
                                'MassDensity',d,...
                                'SpecificHeat',c,'Face',19);                     
thermalProperties(thermalmodelT,'ThermalConductivity',k,...
                                'MassDensity',d,...
                                'SpecificHeat',c,'Face',20);
                            
thermalProperties(thermalmodelT,'ThermalConductivity',k,...
                                'MassDensity',d,...
                                'SpecificHeat',c,'Face',13);                                                  
thermalProperties(thermalmodelT,'ThermalConductivity',k,...
                                'MassDensity',d,...
                                'SpecificHeat',c,'Face',15);                     
thermalProperties(thermalmodelT,'ThermalConductivity',k,...
                                'MassDensity',d,...
                                'SpecificHeat',c,'Face',17);
                            
thermalProperties(thermalmodelT,'ThermalConductivity',k,...
                                'MassDensity',d,...
                                'SpecificHeat',c,'Face',21);                                                  
thermalProperties(thermalmodelT,'ThermalConductivity',k,...
                                'MassDensity',d,...
                                'SpecificHeat',c,'Face',22);                     
thermalProperties(thermalmodelT,'ThermalConductivity',k,...
                                'MassDensity',d,...
                                'SpecificHeat',c,'Face',23);   
   
                            
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
                                'SpecificHeat',cm,'Face',4);
thermalProperties(thermalmodelT,'ThermalConductivity',km,...
                                'MassDensity',dm,...
                                'SpecificHeat',cm,'Face',5);
thermalProperties(thermalmodelT,'ThermalConductivity',km,...
                                'MassDensity',dm,...
                                'SpecificHeat',cm,'Face',6);
                            
thermalProperties(thermalmodelT,'ThermalConductivity',km,...
                                'MassDensity',dm,...
                                'SpecificHeat',cm,'Face',8);
thermalProperties(thermalmodelT,'ThermalConductivity',km,...
                                'MassDensity',dm,...
                                'SpecificHeat',cm,'Face',10);
thermalProperties(thermalmodelT,'ThermalConductivity',km,...
                                'MassDensity',dm,...
                                'SpecificHeat',cm,'Face',12);
                            
                            
%LiF
thermalProperties(thermalmodelT,'ThermalConductivity',klif,...
                                'MassDensity',dlif,...
                                'SpecificHeat',clif,'Face',14);
thermalProperties(thermalmodelT,'ThermalConductivity',klif,...
                                'MassDensity',dlif,...
                                'SpecificHeat',clif,'Face',16);
thermalProperties(thermalmodelT,'ThermalConductivity',klif,...
                                'MassDensity',dlif,...
                                'SpecificHeat',clif,'Face',24);

                          
msh = generateMesh(thermalmodelT,'Hmax',mesh_dens);

figure 
pdeplot(thermalmodelT); %
axis equal
title 'Block With Finite Element Mesh Displayed'
figure
pdegplot(thermalmodelT,'EdgeLabels','on','FaceLabels','on')


%***Solve***
square=true;

if square==true
    thermalIC(thermalmodelT,t0);             
    thermalBC(thermalmodelT,'Edge',51,'Temperature',@BCs1); %8, 12
    thermalBC(thermalmodelT,'Edge',52,'Temperature',@BCs1); %9, 13
    thermalBC(thermalmodelT,'Edge',53,'Temperature',@BCs1); %10, 14

    thermalBC(thermalmodelT,'Edge',54,'Temperature',@BCs1); %8, 12
    thermalBC(thermalmodelT,'Edge',55,'Temperature',@BCs1); %9, 13
    thermalBC(thermalmodelT,'Edge',56,'Temperature',@BCs1); %10, 14

    thermalBC(thermalmodelT,'Edge',57,'Temperature',@BCs1); %8, 12
    thermalBC(thermalmodelT,'Edge',58,'Temperature',@BCs1); %9, 13
    thermalBC(thermalmodelT,'Edge',59,'Temperature',@BCs1); %10, 14
    
    thermalBC(thermalmodelT,'Edge',60,'Temperature',@BCs1); %8, 12
    thermalBC(thermalmodelT,'Edge',61,'Temperature',@BCs1); %9, 13
    thermalBC(thermalmodelT,'Edge',62,'Temperature',@BCs1); %10, 14
else
    thermalIC(thermalmodelT,t0);             
    thermalBC(thermalmodelT,'Edge',51,'Temperature',@BCs1); %8, 12
    thermalBC(thermalmodelT,'Edge',52,'Temperature',@BCs1); %9, 13
    thermalBC(thermalmodelT,'Edge',53,'Temperature',@BCs1); %10, 14

    thermalBC(thermalmodelT,'Edge',54,'Temperature',@BCs1); %8, 12
    thermalBC(thermalmodelT,'Edge',55,'Temperature',@BCs1); %9, 13
    thermalBC(thermalmodelT,'Edge',56,'Temperature',@BCs1); %10, 14

    thermalBC(thermalmodelT,'Edge',57,'Temperature',@BCs1); %8, 12
    thermalBC(thermalmodelT,'Edge',58,'Temperature',@BCs1); %9, 13
    thermalBC(thermalmodelT,'Edge',59,'Temperature',@BCs1); %10, 14
    
    thermalBC(thermalmodelT,'Edge',60,'Temperature',@BCs1); %8, 12
    thermalBC(thermalmodelT,'Edge',61,'Temperature',@BCs1); %9, 13
    thermalBC(thermalmodelT,'Edge',62,'Temperature',@BCs1); %10, 14
end

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

[~,nid1_83] = getClosestNode( msh.Nodes, sfe_83, (9.5*w+9*gap)*10^-6 );
[~,nid2_83] = getClosestNode( msh.Nodes, sfe2_83, (10.5*w+10*gap)*10^-6);
[~,nid3_83] = getClosestNode( msh.Nodes, sfe3_83, (11.5*w+11*gap)*10^-6);

[~,nidmgo1_80] = getClosestNode( msh.Nodes, smgo1_short, (6.5*w+6*gap)*10^-6 );
[~,nidmgo2_80] = getClosestNode( msh.Nodes, smgo2_short, (7.5*w+6*gap)*10^-6 );



%[~,nidfe1] = getClosestNode( msh.Nodes, sfe, (3*w*10^-6)/2);
% [~,nid1m] = getClosestNode( msh.Nodes, smgo, (w*10^-6)/2 );
% [~,nid2m] = getClosestNode( msh.Nodes, smgo2, (1.5*w+gap)*10^-6);
% [~,nid3m] = getClosestNode( msh.Nodes, smgo3, (2.5*w+2*gap)*10^-6);


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
t_data1_73=m.Time(start1:stop1)*10^-9;
t_data2_73=m.Time(start2:stop2)*10^-9;
t_data3_73=m.Time(start3:stop3)*10^-9;
t_data1end_73=m.Time(stop1:stop_final)*10^-9;
t_data2end_73=m.Time(stop2:stop_final)*10^-9;
t_data2end_73=m.Time(stop3:stop_final)*10^-9;
temp1_data73=m.step1_corrected(start1:stop1);
temp2_data73=m.step2_corrected(start2:stop2);
temp3_data73=m.step3_corrected(start3:stop3);
xw = 1; %throughput correction based on slit width and shot number
aeta=25.5; %1/delta_t of SOP data, was 46
a0=481000;% Based on ND filter of 0 (481000) 
t0a=1.909;% Based on ND filter of 0 (1.909) 

%SOP data for s88776
start1=88;
start2=93;
start3=246;
stop1=172;
stop2=418;
stop3=468;
stop_final=750;
opts = detectImportOptions('../s88776_SOP_TP.xlsx');
m = readtable('../s88776_SOP_TP.xlsx',opts);
t_data1_76=m.time(start1:stop1)*10^-9;
t_data2_76=m.time(start2:stop2)*10^-9;
t_data3_76=m.time(start3:stop3)*10^-9;
t_data1end_76=m.time(stop1:stop_final)*10^-9;
t_data2end_76=m.time(stop2:stop_final)*10^-9;
t_data2end_76=m.time(stop3:stop_final)*10^-9;
temp1_data76=m.step1_corrected(start1:stop1);
temp2_data76=m.step2_corrected(start2:stop2);
temp3_data76=m.step3_corrected(start3:stop3);
xw = 1; %throughput correction based on slit width and shot number
aeta=25.5; %1/delta_t of SOP data, was 46
a0_76=481000;% Based on ND filter of 0 (481000) (395900)
t0a_76=1.909;% Based on ND filter of 0 (1.909) (1.91)

%SOP data for s88780
start1=158;
start2=262;
start3=377;
stop1=190;
stop2=370;
stop3=500;
stop_final=750;
opts = detectImportOptions('../s88780_SOP_TP.xlsx');
m = readtable('../s88780_SOP_TP.xlsx',opts);
t_data1_80=m.time(start1:stop1)*10^-9;
t_data2_80=m.time(start2:stop2)*10^-9;
t_data3_80=m.time(start3:stop3)*10^-9;
t_data1end_80=m.time(stop1:stop_final)*10^-9;
t_data2end_80=m.time(stop2:stop_final)*10^-9;
t_data2end_80=m.time(stop3:stop_final)*10^-9;
temp1_data80=m.step1_corrected(start1:stop1);
temp2_data80=m.step2_corrected(start2:stop2);
temp3_data80=m.step3_corrected(start3:stop3);
xw = 1; %throughput correction based on slit width and shot number
aeta=25.5; %1/delta_t of SOP data, was 46
a0=481000;% Based on ND filter of 0 (481000) (395900)
t0a=1.909;% Based on ND filter of 0 (1.909) 


%SOP data for s86483
start1=342;
start2=582;
start3=377;
stop1=412;
stop2=732;
stop3=800;
stop_final=750;
opts = detectImportOptions('../s86483_SOP_TP.xlsx');
m = readtable('../s86483_SOP_TP.xlsx',opts);
t_data1_83=m.time(start1:stop1)*10^-9;
t_data2_83=m.time(start2:stop2)*10^-9;
t_data3_83=m.time(start3:stop3)*10^-9;
t_data1end_83=m.time(stop1:stop_final)*10^-9;
t_data2end_83=m.time(stop2:stop_final)*10^-9;
t_data2end_83=m.time(stop3:stop_final)*10^-9;
temp1_data83=m.step1_corrected(start1:stop1);
temp2_data83=m.step2_corrected(start2:stop2);
temp3_data83=m.step3_corrected(start3:stop3);
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

temp1_83=real(11605*t0a./(log(1+((1-reference)*a0./(aeta*temp1_data83)))));
temp2_83=real(11605*t0a./(log(1+((1-reference)*a0./(aeta*temp2_data83)))));
temp3_83=real(11605*t0a./(log(1+((1-reference)*a0./(aeta*temp3_data83)))));

figure
ax1=axes('Position',[0.13,0.11,0.85,0.69]);
nexttile
plot(t_data1_73*10^9,temp1_73*10^-3,'o',t_data2_73*10^9,temp2_73*10^-3,'o',t_data3_73*10^9,temp3_73*10^-3,'o',...
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

figure
ax1=axes('Position',[0.13,0.11,0.85,0.69]);
plot(t_data1_76*10^9,temp1_76*10^-3,'o',t_data2_76*10^9,temp2_76*10^-3,'o',t_data3_76*10^9,temp3_76*10^-3,'o',...
    tlist*10^9,T1(nid1_76,:)*10^-3,'-',tlist*10^9,T1(nid2_76,:)*10^-3,'-',tlist*10^9,T1(nid3_76,:)*10^-3,'-','Linewidth',2)
legend('Data 1.07 um','Data 2.3 um','Data 3.23 um','Model 1.07 um','Model 2.3 um', 'Model 3.23 um')
grid on
if square==true
    title '1D Model With square Input: s88776 (Fe-LiF)';
else
    title '1D Model With exp Input: s88776 (Fe-LiF)';
end
xlabel 'Time (nanoseconds)'
ylabel 'Temperature (K)'
%xticklabels({'15','20','25','30','35','40','45','50','55'})
%ylim([0 3e4]);


figure
ax1=axes('Position',[0.13,0.11,0.85,0.69]);
plot(t_data1_80*10^9,temp1_80*10^-3,'o',t_data2_80*10^9,temp2_80*10^-3,'o',t_data3_80*10^9,temp3_80*10^-3,'o',...
    tlist*10^9,T1(nid1_80,:)*10^-3,'-',tlist*10^9,T1(nid2_80,:)*10^-3,'-',tlist*10^9,T1(nid3_80,:)*10^-3,'-',...
    'Linewidth',2)
%tlist*10^9,T1(nidmgo1_80,:)*10^-3,'--', tlist*10^9,T1(nidmgo2_80,:)*10^-3,'--',
% plot(t_data1end_80,temp1end_data80,'o')
% alpha(0.2)
legend('Data 0.45 um','Data 1.68 um','Data 2.61 um','Model 0.45 um','Model 1.68 um', 'Model 2.61 um')
%,'Model MgO 1','Model MgO 2')
grid on
if square==true
    title '1D Model With square Input: s88780 (Fe-MgO)';
else
    title '1D Model With exp Input: s88780 (Fe-MgO)';
end
xlabel 'Time (nanoseconds)'
ylabel 'Temperature (K*10^3)'
%xticklabels({'10','15','20','25','30','35','40','45'})
ylim([0 25]);
xlim([18 35]);


figure
ax1=axes('Position',[0.13,0.11,0.85,0.69]);
plot(t_data1_83,temp1_83,'o',t_data2_83,temp2_83,'o',t_data3_83,temp3_83,'o',...
    tlist,T1(nid1_83,:),'-',tlist,T1(nid2_83,:),'-',tlist,T1(nid3_83,:),'-','Linewidth',2)
legend('Data 0.76 um','Data 2.76 um','Data 4.76 um','Model 0.76 um','Model 2.76 um', 'Model 4.76 um')
grid on
if square==true
    title '1D Model With square Input: s86483 (Fe-MgO)';
else
    title '1D Model With exp Input: s86483 (Fe-MgO)';
end
xlabel 'Time (nanoseconds)'
ylabel 'Temperature (K)'
%xticklabels({'15','20','25','30','35','40','45','50','55'})
ylim([0 3e4]);

% figure
% pdeplot(thermalmodelT,'XYData',T1(:,end),'Contour','on','ColorMap','hot');
% figure
% ax1=axes('Position',[0.13,0.11,0.85,0.69]);
% plot(t_data,temp1_corrected,'o',t_data,temp2_corrected,'o',t_data,temp3_corrected,'o',...
%     tlist,T1(nid1m,:),'-',tlist,T1(nid2m,:),'-',tlist,T1(nid3m,:),'-','Linewidth',2)
% legend('Data 1 um','Data 2 um','Data 3 um','MgO interface 1um','MgO interface 2um', 'MgO interface 3um')
% grid on
% title '2D Model: MgO side of interface';
% xlabel 'Time (nanoseconds)'
% ylabel 'Temperature (K)'
% %xticklabels({'15','20','25','30','35','40','45','50','55'})
% ylim([0 3e4]);

T11=T1(nid1_73,:);
T12=T1(nid2_73,:);
T13=T1(nid3_73,:);
T21=T1(nid1_76,:);
T22=T1(nid2_76,:);
T23=T1(nid3_76,:);
T31=T1(nid1_80,:);
T32=T1(nid2_80,:);
T33=T1(nid3_80,:);
T41=T1(nid1_83,:);
T42=T1(nid2_83,:);
T43=T1(nid3_83,:);

save('combos1.mat','tlist','T11','T12','T13','T21','T22','T23','T31','T32','T33','T41','T42','T43','t_data1_73','t_data2_73','t_data3_73','temp1_73','temp2_73','temp3_73','t_data1_76','t_data2_76','t_data3_76','temp1_76','temp2_76','temp3_76','t_data1_80','t_data2_80','t_data3_80','temp1_80','temp2_80','temp3_80','t_data1_83','t_data2_83','t_data3_83','temp1_83','temp2_83','temp3_83')