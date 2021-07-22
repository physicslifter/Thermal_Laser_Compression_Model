%Main script for analyzing shot s88773
%two different mesh setups with iron and MgO blocks in contact at ~1 um and
%~2 um (corresponds to the ~2 and ~3 um steps)

%****     2um     *****
%load input parameters
load('inputs/input_matrix.mat','a','b')

%load input data
load('fitting_data.mat','a0','aeta','init_fit','m','mr','opts','optsr','p','ref2','ref2_fit','ref2_fourier','start','stop','t0a','t_data','t_ref','temp1_corrected','temp1_data','temp2_corrected','temp2_data','temp3_corrected','temp3_data','xw')

%***** Parameters *******
k = @(~,state) b+a*state.u; %W/mK - 
% (previous fit): 0.045*state.u + 3000./sqrt(abs(state.u)), (de koker fit):
% 30+0.0288t, (good linear fit): 120+0.02t
d=12500; %kg/m^3 - determined from Smith et al. 2018 isentrope (previous at 11000)
c=500; %J/kg K
t0=2000; %K/counts
w=0.2; % width (um)
mesh_dens = 0.1e-6; % spacing of mesh points (um)


km=100; %W/mK used to be 5
dm=5800; %kg/m^3 (4200) from Jin et al. 2010
cm=800; %J/kg K (800)
t0m=1000; %K - old value 1000

tlist = [18:0.5:40]*10^-9;


%**** Creating geometry ****

ratio = 7900/d; %ambient density divided by measured density
inc = 0.1; %amount measuring into window

lfe=ratio*1.07; % Length of iron (um) 0.89 0.72 0.77
lmgo=2; % length of mgo (um)

lfe2=ratio*2.3; % Length of iron (um) 1.57
lmgo2=2; % length of mgo (um)

lfe3=ratio*3.23; % Length of iron (um) 1.57
lmgo3=2; % length of mgo (um)

finc = 0.1;
sfe=(lfe-finc)*10^-6; % Measurement spot in iron
smgo=(lfe+inc)*10^-6; % Measurement spot in mgo

sfe2=(lfe2-finc)*10^-6; % Measurement spot in iron
smgo2=(lfe2+inc)*10^-6; % Measurement spot in mgo

sfe3=(lfe3-finc)*10^-6; % Measurement spot in iron
smgo3=(lfe3+inc)*10^-6; % Measurement spot in mgo

gap = 0.2;
w2=w+gap;
w3=2*w+gap;
w4=2*w+2*gap;
w5=3*w+2*gap;

r = [3 4 0 lfe*10^-6 lfe*10^-6 0 0 0 w*10^-6 w*10^-6]';
rmgo = [3 4 lfe*10^-6 (lfe+lmgo)*10^-6 (lfe+lmgo)*10^-6 lfe*10^-6 0 0 w*10^-6 w*10^-6]';

r2 = [3 4 0 lfe2*10^-6 lfe2*10^-6 0 w2*10^-6 w2*10^-6 w3*10^-6 w3*10^-6]';
rmgo2 = [3 4 lfe2*10^-6 (lfe2+lmgo2)*10^-6 (lfe2+lmgo2)*10^-6 lfe2*10^-6 w2*10^-6 w2*10^-6 w3*10^-6 w3*10^-6]';

r3 = [3 4 0 lfe3*10^-6 lfe3*10^-6 0 w4*10^-6 w4*10^-6 w5*10^-6 w5*10^-6]';
rmgo3 = [3 4 lfe3*10^-6 (lfe3+lmgo3)*10^-6 (lfe3+lmgo3)*10^-6 lfe3*10^-6 w4*10^-6 w4*10^-6 w5*10^-6 w5*10^-6]';

thermalmodelT = createpde('thermal','transient');

gdm2 = [r rmgo r2 rmgo2 r3 rmgo3];
sf2 = 'r+rmgo+r2+rmgo2+r3+rmgo3';
ns2=char('r','rmgo','r2','rmgo2','r3','rmgo3');
ns2=ns2';
g2 = decsg(gdm2,sf2,ns2);
geometryFromEdges(thermalmodelT,g2);

thermalProperties(thermalmodelT,'ThermalConductivity',k,...
                                'MassDensity',d,...
                                'SpecificHeat',c,'Face',4);                                                  
thermalProperties(thermalmodelT,'ThermalConductivity',k,...
                                'MassDensity',d,...
                                'SpecificHeat',c,'Face',5);                     
thermalProperties(thermalmodelT,'ThermalConductivity',k,...
                                'MassDensity',d,...
                                'SpecificHeat',c,'Face',6);

thermalProperties(thermalmodelT,'ThermalConductivity',km,...
                                'MassDensity',dm,...
                                'SpecificHeat',cm,'Face',1);
thermalProperties(thermalmodelT,'ThermalConductivity',km,...
                                'MassDensity',dm,...
                                'SpecificHeat',cm,'Face',2);
thermalProperties(thermalmodelT,'ThermalConductivity',km,...
                                'MassDensity',dm,...
                                'SpecificHeat',cm,'Face',3);

                          
msh = generateMesh(thermalmodelT,'Hmax',mesh_dens);

% figure 
% pdeplot(thermalmodelT); %
% axis equal
% title 'Block With Finite Element Mesh Displayed'
% figure
% pdegplot(thermalmodelT,'EdgeLabels','on','FaceLabels','on')


%***Solve***


thermalIC(thermalmodelT,t0);             
thermalBC(thermalmodelT,'Edge',8,'Temperature',@BC_external_exp); %8, 12
thermalBC(thermalmodelT,'Edge',9,'Temperature',@BC_external_exp); %9, 13
thermalBC(thermalmodelT,'Edge',10,'Temperature',@BC_external_exp); %10, 14



R1 = solve(thermalmodelT,tlist);
T1 = R1.Temperature;

getClosestNode = @(p,x,y) min((p(1,:) - x).^2 + (p(2,:) - y).^2);
[~,nid1] = getClosestNode( msh.Nodes, sfe, (w*10^-6)/2 );
[~,nid2] = getClosestNode( msh.Nodes, sfe2, (1.5*w+gap)*10^-6);
[~,nid3] = getClosestNode( msh.Nodes, sfe3, (2.5*w+2*gap)*10^-6);
%[~,nidfe1] = getClosestNode( msh.Nodes, sfe, (3*w*10^-6)/2);

[~,nid1m] = getClosestNode( msh.Nodes, smgo, (w*10^-6)/2 );
[~,nid2m] = getClosestNode( msh.Nodes, smgo2, (1.5*w+gap)*10^-6);
[~,nid3m] = getClosestNode( msh.Nodes, smgo3, (2.5*w+2*gap)*10^-6);


%******    Plotting    ******



figure
ax1=axes('Position',[0.13,0.11,0.85,0.69]);
plot(t_data,temp1_corrected,'o',t_data,temp2_corrected,'o',t_data,temp3_corrected,'o',...
    tlist,T1(nid1,:),'-',tlist,T1(nid2,:),'-',tlist,T1(nid3,:),'-','Linewidth',2)
legend('Data 1 um','Data 2 um','Data 3 um','Fe 1um','Fe 2um', 'Fe 3um')
grid on
title '1D Model With Exp Input: MgO';
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


