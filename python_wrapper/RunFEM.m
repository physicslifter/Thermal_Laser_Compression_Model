function [] = RunFEM(input_parameter_file,run_name)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% This code takes the raw data from an excel sheet and converts to temperature.
% This version uses refelctivity information. This code should be run
% before heat equation code.

load(input_parameter_file,'sop_lineouts_data','reflectivity_data','xw','aeta','a0','t0a','d','c','t0','w','km','dm','cm','t0m','peak_temp')

%SOP data%
start=150;
stop=890;
opts = detectImportOptions(sop_lineouts_data);
m = readtable(sop_lineouts_data,opts);
t_data=m.Time_ns_(start:stop)*10^-9;
temp1_data=m.x1_07um(start:stop);
temp2_data=m.x2_3um(start:stop);
temp3_data=m.x3_23um(start:stop);

%Reflectivity data
optsr = detectImportOptions(reflectivity_data);
mr = readtable(reflectivity_data,optsr);
t_ref = mr.Time_ns_*10^-9;

ref2=mr.Reflectivity_1;
ref2_fourier=fit(t_ref,ref2,'fourier7');
ref2_fit=ref2_fourier(t_data);
plot(t_ref, ref2, 'o', t_data, ref2_fit)


%Converting from counts to temp
temp1_corrected=real(11605*t0a./(log(1+((1-ref2_fit)*a0./(aeta*temp1_data)))));
temp2_corrected=real(11605*t0a./(log(1+((1-ref2_fit)*a0./(aeta*temp2_data)))));
temp3_corrected=real(11605*t0a./(log(1+((1-ref2_fit)*a0./(aeta*temp3_data)))));


% The following code is for the "internal boundary condition" method. This
% this should probably change to interpolation
[p] = polyfit(t_data,temp1_corrected,9);
init_fit = polyval(p,t_data);

%This loop writes out the polynomial so I can copy and paste it into the
%HeatedBlock Function as the boundary condition
eqn='';
for c = p
    eqn=append(eqn,num2str(c),'*t^',num2str(length(p)-find(p==c)),' + ');
end
%eqn
%figure
%plot(t_data,temp1_corrected,'o',t_data,init_fit,'-')
%legend('data','fit')

% Main script for analyzing shot s88773
% Full 2D model of experiment
% utilizes square wave temp input

%****     2um     *****

%***** Parameters *******
%Sample
k = @(~,state) 30+0.035*state.u; %W/mK linear temperature dependent model for thermal conductivity
% (previous fit): 0.045*state.u + 3000./sqrt(abs(state.u)), (de koker fit):
% 30+0.0288t, (good linear fit): 120+0.02t

tlist = [18:0.5:35]*10^-9;


%**** Creating geometry ****

ratio = 7900/d; %ambient density divided by measured density
inc = 0.05; %amount measuring into window

lfe=ratio*1.07; % Length of iron (um) 0.89 0.72 0.77
lmgo=2; % length of mgo (um)

lfe2=ratio*2.3; % Length of iron (um) 1.57
lmgo2=2; % length of mgo (um)

lfe3=ratio*3.23; % Length of iron (um) 1.57
lmgo3=2; % length of mgo (um)

sfe=(lfe-0.01)*10^-6; % Measurement spot in iron
smgo=(lfe+inc)*10^-6; % Measurement spot in mgo

sfe2=(lfe2-0.01)*10^-6; % Measurement spot in iron
smgo2=(lfe2+inc)*10^-6; % Measurement spot in mgo

sfe3=(lfe3-0.01)*10^-6; % Measurement spot in iron
smgo3=(lfe3+inc)*10^-6; % Measurement spot in mgo


r = [3 4 0 lfe*10^-6 lfe*10^-6 0 0 0 w*10^-6 w*10^-6]';
rmgo = [3 4 lfe*10^-6 (lfe+lmgo)*10^-6 (lfe+lmgo)*10^-6 lfe*10^-6 0 0 w*10^-6 w*10^-6]';

r2 = [3 4 0 lfe2*10^-6 lfe2*10^-6 0 w*10^-6 w*10^-6 2*w*10^-6 2*w*10^-6]';
rmgo2 = [3 4 lfe2*10^-6 (lfe2+lmgo2)*10^-6 (lfe2+lmgo2)*10^-6 lfe2*10^-6 w*10^-6 w*10^-6 2*w*10^-6 2*w*10^-6]';

r3 = [3 4 0 lfe3*10^-6 lfe3*10^-6 0 2*w*10^-6 2*w*10^-6 3*w*10^-6 3*w*10^-6]';
rmgo3 = [3 4 lfe3*10^-6 (lfe3+lmgo3)*10^-6 (lfe3+lmgo3)*10^-6 lfe3*10^-6 2*w*10^-6 2*w*10^-6 3*w*10^-6 3*w*10^-6]';

thermalmodelT = createpde('thermal','transient');

gdm2 = [r rmgo r2 rmgo2 r3 rmgo3];
sf2 = 'r+rmgo+r2+rmgo2+r3+rmgo3';
ns2=char('r','rmgo','r2','rmgo2','r3','rmgo3');
ns2=ns2';
g2 = decsg(gdm2,sf2,ns2);
geometryFromEdges(thermalmodelT,g2);

thermalProperties(thermalmodelT,'ThermalConductivity',k,...
                                'MassDensity',d,...
                                'SpecificHeat',c,'Face',2);                           
                            
thermalProperties(thermalmodelT,'ThermalConductivity',k,...
                                'MassDensity',d,...
                                'SpecificHeat',c,'Face',4);
                                                        
thermalProperties(thermalmodelT,'ThermalConductivity',k,...
                                'MassDensity',d,...
                                'SpecificHeat',c,'Face',5);

thermalProperties(thermalmodelT,'ThermalConductivity',km,...
                                'MassDensity',dm,...
                                'SpecificHeat',cm,'Face',1);
thermalProperties(thermalmodelT,'ThermalConductivity',km,...
                                'MassDensity',dm,...
                                'SpecificHeat',cm,'Face',3);
thermalProperties(thermalmodelT,'ThermalConductivity',km,...
                                'MassDensity',dm,...
                                'SpecificHeat',cm,'Face',6);

                          
msh = generateMesh(thermalmodelT,'Hmax',0.05e-6);
%figure 
%pdeplot(thermalmodelT); %
%axis equal
%title 'Block With Finite Element Mesh Displayed'
%figure
%pdegplot(thermalmodelT,'EdgeLabels','on','FaceLabels','on')


%***Solve***

% This is where we call our temperature input function
thermalIC(thermalmodelT,t0);             
thermalBC(thermalmodelT,'Edge',12,'Temperature',@transientBCHeatedBlock_square);
thermalBC(thermalmodelT,'Edge',13,'Temperature',@transientBCHeatedBlock_square);
thermalBC(thermalmodelT,'Edge',14,'Temperature',@transientBCHeatedBlock_square);


R1 = solve(thermalmodelT,tlist);
T1 = R1.Temperature;

getClosestNode = @(p,x,y) min((p(1,:) - x).^2 + (p(2,:) - y).^2);
[~,nid1] = getClosestNode( msh.Nodes, sfe, (w*10^-6)/2 );
[~,nid2] = getClosestNode( msh.Nodes, sfe2, (3*w*10^-6)/2 );
[~,nid3] = getClosestNode( msh.Nodes, sfe3, (5*w*10^-6)/2 );
%[~,nidfe1] = getClosestNode( msh.Nodes, sfe, (3*w*10^-6)/2);

[~,nid1m] = getClosestNode( msh.Nodes, smgo, (w*10^-6)/2 );
[~,nid2m] = getClosestNode( msh.Nodes, smgo2, (3*w*10^-6)/2 );
[~,nid3m] = getClosestNode( msh.Nodes, smgo3, (5*w*10^-6)/2 );

%Naming some important variables
T11=T1(nid1,:);
T12=T1(nid2,:);
T13=T1(nid3,:);

%******    Plotting    ******



% figure
% ax1=axes('Position',[0.13,0.11,0.85,0.69]);
% plot(t_data,temp1_corrected,'o',t_data,temp2_corrected,'o',t_data,temp3_corrected,'o',...
%     tlist,T1(nid1,:),'-',tlist,T1(nid2,:),'-',tlist,T1(nid3,:),'-','Linewidth',2)
% legend('Data 1 um','Data 2 um','Data 3 um','Fe interface 1um','Fe interface 2um', 'Fe interface 3um')
% grid on
% title '2D Model: MgO side of interface';
% xlabel 'Time (nanoseconds)'
% ylabel 'Temperature (K)'
% xticklabels({'15','20','25','30','35','40','45','50','55'})
% ylim([0 3e4]);
% 
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

%Save variables in .mat file
mydate=date;
mytime=string(clock);
filename=strcat('FEM_output/',run_name,'.mat');
save(filename,'tlist','T11','T12','T13','t_data','temp1_corrected', 'temp2_corrected','temp3_corrected')
end

