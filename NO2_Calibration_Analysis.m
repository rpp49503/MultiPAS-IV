%% Inputs

Num_Backgrounds=2;

cross_sec=[0;1.4696E-19;3.6142E-21;6.9597E-22]; % Literature cross-section values (532, 662, and 780 nm respectively)

% Setting variable name for filepath for easier access (less user input). Enter date as "yyyymmdd".
file_path = "Z:\Groups\Smith_G\Ryan Poland\MultiPAS-IV\NO2 Calibration\Data\20220513";

%% Tabulate text files

files=dir(file_path + '\*.txt');
files_2=fullfile(file_path,{files.name});

% Loop to tabulate all data from folder into new cell
Datainput={length(files_2)};
for n=1:length(files_2)
    Datainput{n}=readtable(files_2(n));
end

%% Number Densities

NO2_flows=Datainput{1}; % Create new table for NO2 flow data for easier access (may be unnecessary)
d_f = NO2_flows.NO2Flow ./ (NO2_flows.NO2Flow + NO2_flows.N2Flow); % ratio of flow for NO2
NO2_cell = (NO2_flows.TankConc_(1)*1E-6)*((6.022141e23*(NO2_flows.Pressure(1)*0.000986923))/((NO2_flows.Temperature(1)+273.15)*82.057338))*d_f; % Concentration in chamber

%% Theoretical(Burrows) and CRD Absorptions

theo_abs_662=NO2_cell*cross_sec(3)*1E8; % Calculate theoretical absorption at 662nm in Mm^-1

% Tabulate tau values for all runs
tau_t=[length(Datainput)];
for n=2:length(Datainput)
    tau=Datainput{1,n};
    tau_t(:,n)=tau.tau_sec(3);
end

Taulabels={files.name}; % Designate variable names for tau_tt columns
tau_tt=array2table(tau_t,"VariableNames",Taulabels); % Re-shape tau values into table for easier indexing

% Pull signals from background runs
s_0_t=[];
for z=2:Num_Backgrounds+1
    s_I0=Datainput{1,z};
    s_0_t(:,n)=s_I0.mic_mV;
end

% Pull and compare tau0 values from background runs
for k=z
    I0_t(:,k)=tau_tt.("I0_"+k+".txt");
end

tau0_diff=abs(1-(min(I0_t)/max(I0_t)));

% Alert if difference in tau0 values is greater than 1%
if tau0_diff > 0.01
    fig=uifigure;
    uialert(fig,'Tau values not in agreement (greater than 1%)','Error')
    s_0=
else ;tau0=(tau_tt.("I0_1.txt")+tau_tt.("I0_2.txt"))/2; % Average tau0 if in agreement within 1%
    s_0=sum(s_0_t,2)/2;
end

% Compute CRD Absorption from Tau values
CRD_abs=(1/2.99792458E8)*((1./tau_t)-(1/tau0))*1E6;
CRD_abs_t=array2table(CRD_abs,"VariableNames",Taulabels); % Tabulate CRD absorptions for easier viewing (non-essential)
CRD_abs=CRD_abs(4:length(CRD_abs));

%% Signal to Power Calculations

% Pulling signal and power from "Datainput"
signal_t=[];
power_t=[];
for i=1:4
    for n=4:length(Datainput)
    signal=Datainput{1,n};
    signal_t(i,n)=signal.mic_mV(i);
    power_t(i,n)=signal.pd_mV(i);
    end
end

power_tt=array2table(power_t,"VariableNames",Taulabels); % Table for viewing power measurements (non-essential)
signal_tt=array2table(signal_t,"VariableNames",Taulabels); % Table for viewing signal measurements (non-essential)

% Correct all signals for background (s_0 defined in section above)
s_s_0=[];
for n=2:4
    s_s_0(n,:)=signal_t(n,4:length(signal_t))-s_0(n);
end

% Calculating S/P and correcting for cross-sections at each wavelength
s_p=[];
for z=2:4
    s_p(z,:)=((s_s_0(z,:)/1000)./(power_t(z,4:length(power_t))/1000))*(cross_sec(3)/cross_sec(z));
end

%% Calibration Plot

x=CRD_abs; % define 'x' for shortened notation

figure(1)

scatter(x,s_p(2,:),'filled','MarkerFaceColor',[0 1 0]) % Plotting absorption vs. s/p
mdl_532=fitlm(x,s_p(2,:)); % Create linear model for 532 nm
hold on
plot(x,mdl_532.Fitted,'g','LineWidth',1) % Plot linear model from "fitted" data.

% Repeat for other lines
scatter(theo_abs_662,s_p(3,:),'filled','MarkerFaceColor',[0.6 0.6 0.6]) % 662 nm from flows
mdl_662_flow=fitlm(x,s_p(3,:));
plot(x,mdl_662_flow.Fitted,'Color',[0.6 0.6 0.6],'LineWidth',1)

scatter(x,s_p(3,:),'filled','MarkerFaceColor' ,'r') % 662 nm
mdl_662=fitlm(x,s_p(3,:));
plot(x,mdl_662.Fitted,'r','LineWidth',1)

scatter(x,s_p(4,:),'filled','MarkerFaceColor',[0.5 0 0]) % 780 nm
mdl_780=fitlm(x,s_p(4,:));
plot(x,mdl_780.Fitted,'Color',[0.5 0 0],'LineWidth',1)

xlabel('Absorption_{662 nm} (Mm^{-1})')
ylabel('(Mic Signal(V) / Power(W)) * (\sigma 662 nm / \sigma x nm)')
title('NO_2 Calibration Plot')
legend('532 nm','','662 nm (flow)','','662 nm','','780 nm','Location','northwest')

hold off

%% Pulling calibration constants into new file

% Extracting slopes of linear fit lines from 'fitlm' models
calibration_532=mdl_532.Coefficients.Estimate(2);
calibration_662_flow=mdl_662_flow.Coefficients.Estimate(2);
calibration_662=mdl_662.Coefficients.Estimate(2);
calibration_780=mdl_780.Coefficients.Estimate(2);

% Create table of constants and write as text file to dated folder
calibration_t=table(calibration_532,calibration_662_flow,calibration_662,calibration_780);
writetable(calibration_t,file_path + '\calibration_constants','Delimiter','\t');


