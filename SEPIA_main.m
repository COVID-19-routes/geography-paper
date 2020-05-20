% Code to run the scenarios shown in Figure 2 of the preprint:
% The geography of COVID-19 spread in Italy and implications for the relaxation of confinement measures
% link: https://www.medrxiv.org/content/10.1101/2020.04.30.20083568v2.full.pdf+html
% 
% Enrico Bertuzzo et al., 20/05/2020; mail to enrico.bertuzzo@unive.it
clear all
close all

addpath('./input_data')
addpath('./functions')

%% READ POPULATION DATA

load Population_data.mat

% Content

% N    : Population of each Italian provice
% n    : number of Italian provinces
% n_reg: number of Italian regions
% p    : mobile population fraction for each Italian province
% q    : n*n mobility matrix. Element qij represents the fraction of the
%        mobile population living in province i that moves to province j
% prov_name  : name of the provinces
% reg_name   : name of the regions
% prov_IDprov: Census ID of the province
% prov_IDreg : Census ID of the region each province belongs to
% prov2reg   : n_reg*n matrix. Element ij is equal to 1 if province j 
%              belongs to region i

% Transfer data into V structure
V=Population_data;
clear Population_data

%% READ EPIDEMIOLOGICAL DATA
load('Epi_data.mat')

% Content

% Date     : date of epidemiologica data. Matlab datenum format
% prov_Hnew: n*length(Date) Matrix, reconstructed daily new hospitalized
%            cases for each province
% prov_Hnew: n_reg*length(Date) Matrix, reconstructed daily new 
%            hospitalized cases for each region

% Transfer data into V structure
V=append_struct(V,Epi_data);
clear Epi_data

%% READ MOBILITY REDUCTION DATA

load MobilityReduction_data.mat

%Content

% Information are extracted from preprint
% https://www.medrxiv.org/content/10.1101/2020.03.22.20039933v2 and
% represents the reduction in extraprovice mobility during the lockdown
% in Italy.

% tmob: Date (Matlab datenum format) for which values of mobility are 
%       provided
% mob : n*length(tmob) matrix. Element ij represents the ratio of
%       extra-province mobility of provice i at time tmob(j) with respect 
%       to the pre-COVID mobility
%
% Timeserire of mobility is then obtained through linear interpolation of
% the values provided in tmob and mob.

% Transfer data into V structure
V=append_struct(V,MobilityReduction_data);
clear MobilityReduction_data

%% READ MODEL PARAMETERS

load Parameters.mat

% Content

% PAR_names        : Parameter names (see preprint text)
% Posteriror_sample: Posterior samples of PAR_names parameters (note that 
%                    some parameters have a fixed values, see main text)
% nPAR_model       : number of structural model parameters (first 
%                    nPAR_model parameters of PAR_names)
% x0               : initial conditions (1 exposed individual in the 
%                    province of Lodi)
% seeding          : indexes of province for which initial conditions are 
%                    also estimated

% Additional fixed parameters
% gammaA_over_gammaQ: ratio \gamma_A/\gamma_Q
% gammaQ_over_gammaH: ratio \gamma_Q/\gamma_H
% zeta: fraction of cases that are quarantined (see preprint text)

% Transfer data into V structure
V=append_struct(V,Parameters);
clear Parameters


%% SIMULATION SET UP

NSample=100; %number of sample from  posterior distribution of parameters
NNoise4sample=20; %number of sample from the negative binomial distribution for each posterior sample
NReal=NSample*NNoise4sample; %total number of realizations
quant=[0.025 0.25 0.5 0.75 0.975]; %quantiles of interests

V.time_model_final=datenum('6/30/2020')+1; %simulation horizon
betaInc=1;  %increase in trasmission after May 3, 2020 (Figure 2 of the preprint is produced using value 1, 1.2 and 1.4)

% Apply increase in trasmission in all provinces
V.ProvBetaInc=ones(V.n,1);
V.ProvBetaInc=betaInc;

% Define schedule of variation for beta transmission parameters
V.tbeta=zeros(1,9);
V.tbeta(1)=datenum('1/01/2020');
V.tbeta(2)=datenum('2/24/2020');
V.tbeta(3)=datenum('2/26/2020');
V.tbeta(4)=datenum('3/08/2020');
V.tbeta(5)=datenum('3/11/2020');
V.tbeta(6)=datenum('3/22/2020');
V.tbeta(7)=datenum('5/04/2020');
V.tbeta(8)=datenum('5/07/2020');
V.tbeta(9)=datenum('12/31/2020');

% String command to be used to define beta variations
V.beta_string='[on, on, betaP1P0*on, betaP1P0*on, betaP1P0*betaP2P1*on, betaP1P0*betaP2P1*betaP3P2.*on, betaP1P0*betaP2P1*betaP3P2.*on, betaP1P0*betaP2P1*betaP3P2.*on.*V.ProvBetaInc, betaP1P0*betaP2P1*betaP3P2.*on.*V.ProvBetaInc ]';

% Time_resample: time at which output are provided
time_resample=V.Date(1):V.time_model_final-1;

% Define colorscheme
V.col=get(groot,'defaultAxesColorOrder');

%% MAIN

% Preallocation of variables
reg_Hnew_real=zeros(V.n_reg,length(time_resample)-1,NReal); %simulated new hospitalized cases for each region (1st dimension), day (2nd dimension) and realization (third dimension)

cont_real=0; %inizialization of realization counter
for cont_sample=1:NSample
    display([' Posterior sample ',num2str(cont_sample),' of ',num2str(NSample)]);
    % Run model for this posterior sample
    PAR_real=V.Posterior_sample(ceil(rand*length(V.Posterior_sample)),:); %select one random posterior sample
    [x,V]=SEPIA(PAR_real,V); %run model
    
    prov_cumH=interp1(V.time_model,x(:,10*V.n+1:11*V.n),time_resample)'; %cumulative hospitalized cases for each province (1st dimension) and day (2nd dimension) of time_resample
    prov_Hnew=diff(prov_cumH,1,2); %new hospitalized cases for each province (1st dimension) and day (2nd dimension) of time_resample
    
    % Add negative binomial error 
    omega=PAR_real(12); %parameter of the NB1 parametrization of the negative binomial distribution
    % Compute r and p parametr of the negative binomila distribution
    r=prov_Hnew/(omega-1);
    p=1/omega;
    for cont_noise=1:NNoise4sample
        cont_real=cont_real+1; %update counter of realizations   
        
        % Add error at the province scale and upscale results at the
        % regional scale
        reg_Hnew_real(:,:,cont_real)=V.prov2reg*nbinrnd(r,p);    
    end
end

% Compute relevant quantiles at the regional level
regHnew_quantile=quantile(reg_Hnew_real,quant,3);
% Compute relevant quantiles at the country level
countryHnew_quantile=quantile(squeeze(sum(reg_Hnew_real)),quant,2);
  
%% EDITING RESULTS

RegFig1=[3 8 1 5 11]; %regions to be displayed in the first figure (Figure 2 of the preprint)
RegFig2=[2 4 6 7 9 10 12:20]; %regions to be displayed in the second figure  (Figure S1 of the preprint)

% figure: model fit against the data of hospitalized individuals for the 
%         whole Italy and the most affected regions.
figure
set(gcf,'color','w')
subplot(3,2,1)
patch([time_resample(2:end) flip(time_resample(2:end))]',[countryHnew_quantile(:,1); flip(countryHnew_quantile(:,5))],V.col(1,:),'edgecolor','n','facealpha',0.15,'Handlevisibility','off')
hold on
patch([time_resample(2:end) flip(time_resample(2:end))]',[countryHnew_quantile(:,2); flip(countryHnew_quantile(:,4))],V.col(1,:),'edgecolor','n','facealpha',0.25,'Handlevisibility','off')
plot(V.Date(2:end),sum(V.reg_Hnew),'o','color',V.col(2,:),'markersize',5,'markerfacecolor','w','linewidth',1.5,'displayname','New Hospitalized')
plot(time_resample(2:end),countryHnew_quantile(:,3),'color',V.col(1,:),'linewidth',1.5)
box off
set(gca,'Xlim',[time_resample(1) time_resample(end)])
datetick('x','keeplimits','keepticks')
ylabel('New Hospitalized Cases')
title('Italy')

for i=1:length(RegFig1)
    r=RegFig1(i);
    subplot(3,2,i+1)
    patch([time_resample(2:end) flip(time_resample(2:end))]',...
        [squeeze(regHnew_quantile(r,:,1))'; flip(squeeze(regHnew_quantile(r,:,5)))'],V.col(1,:),'edgecolor','n','facealpha',0.15,'Handlevisibility','off')
    hold on
    patch([time_resample(2:end) flip(time_resample(2:end))]',...
        [squeeze(regHnew_quantile(r,:,2))'; flip(squeeze(regHnew_quantile(r,:,4)))'],V.col(1,:),'edgecolor','n','facealpha',0.25,'Handlevisibility','off')
    plot(V.Date(2:end),V.reg_Hnew(r,:),'o','color',V.col(2,:),'markersize',5,'markerfacecolor','w','linewidth',1.5,'displayname','New Hospitalized')
    plot(time_resample(2:end),squeeze(regHnew_quantile(r,:,3)),'color',V.col(1,:),'linewidth',1.5)
    box off
    set(gca,'Xlim',[time_resample(1) time_resample(end)])
    datetick('x','keeplimits','keepticks')
    box off
    title(V.reg_name(r))    
end

% figure: model fit against the data of hospitalized individuals for the 
%         less affected regions.
figure
set(gcf,'color','w')
for i=1:length(RegFig2)
    r=RegFig2(i);
    subplot(5,3,i)
    patch([time_resample(2:end) flip(time_resample(2:end))]',...
        [squeeze(regHnew_quantile(r,:,1))'; flip(squeeze(regHnew_quantile(r,:,5)))'],V.col(1,:),'edgecolor','n','facealpha',0.15,'Handlevisibility','off')
    hold on
    patch([time_resample(2:end) flip(time_resample(2:end))]',...
        [squeeze(regHnew_quantile(r,:,2))'; flip(squeeze(regHnew_quantile(r,:,4)))'],V.col(1,:),'edgecolor','n','facealpha',0.25,'Handlevisibility','off')
    plot(V.Date(2:end),V.reg_Hnew(r,:),'o','color',V.col(2,:),'markersize',5,'markerfacecolor','w','linewidth',1.5,'displayname','New Hospitalized')
    plot(time_resample(2:end),squeeze(regHnew_quantile(r,:,3)),'color',V.col(1,:),'linewidth',1.5)
    box off
    set(gca,'Xlim',[time_resample(1) time_resample(end)])
    datetick('x','keeplimits','keepticks')
    box off
    title(V.reg_name(r))
end
