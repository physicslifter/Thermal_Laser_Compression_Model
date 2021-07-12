% This code takes the raw data from an excel sheet and converts to temperature.
% This version uses refelctivity information. This code should be run
% before heat equation code.

clear

%SOP data%
start=150;
stop=890;

%Pat: altered the code here to read to the proper filepaths for the git
%repository. (The files are in the "Data" subfolder)
opts = detectImportOptions('../Data/s88773_sop_lineouts_TP.xlsx');
m = readtable('../Data/s88773_sop_lineouts_TP.xlsx',opts);

t_data=m.Time_ns_(start:stop)*10^-9;

%Pat: changed from Um to um in the variable names, in order to get the
%filenames to mesh
temp1_data=m.x1_07um(start:stop);
temp2_data=m.x2_3um(start:stop);
temp3_data=m.x3_23um(start:stop);

xw = 1; %throughput correction based on slit width and shot number
aeta=25.5; %1/delta_t of SOP data, was 46
a0=481000;% Based on ND filter
t0a=1.909;% Based on ND filter

%Reflectivity data
optsr = detectImportOptions('../Data/Reflectivity_s88773.xlsx');
mr = readtable('../Data/Reflectivity_s88773.xlsx',optsr);
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
[p] = polyfit(t_data,temp1_corrected,9)
init_fit = polyval(p,t_data);

%This loop writes out the polynomial so I can copy and paste it into the
%HeatedBlock Function as the boundary condition
eqn='';
for c = p
    eqn=append(eqn,num2str(c),'*t^',num2str(length(p)-find(p==c)),' + ');
end
eqn
figure
plot(t_data,temp1_corrected,'o',t_data,init_fit,'-')
legend('data','fit')

