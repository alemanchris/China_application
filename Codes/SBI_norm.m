%{
DESCRIPTION
SBI_norm calls functions:
-smooth_ts: Smooths series, using a cubic spline
-norm_func: Performs normalization CtoT, generates policy estimates and
figures

INPUT: 
-Outcome Variable TxN Matrix, where T(time observations e.g year) N(number of regions e.g US States)
-Population by region TxN
-Time units Tx1 (Years, Months Days, etc)
-Names by regions

OUTPUT: 
-psi_esti: Nx3 Mapping parameters 
-pol_esti: Nx1 Policy Estimate in log difference
-tp_norm1:  Nx1 Location of policy date after normalization
-fp_norm:  Nx1 Location of the first point in C to time of T
-l_opt:    Nx1 Lenght of the overlap interval 
-flag :    Nx1 1:Mapping didnt converge 2:Mapping delivers non-sence
-Figures: 5 sets per region (both in .png and .eps formats)
    fig_1: Mapping function
    fig_2: Before Normalization
    fig_3: After Normalization
    fig_4: Zoom Overlap Interval
    fig_5: Zoom Overlap Interval Policy Effect in log difference
Figures have the option to have the mapping parameters and flag in the
title, for this set opt.title_coeff=1

REMARKS: 
-Trimming(anchor) is optional, to trim select par.cut>0
-Local minizer is implemented, thus no need for penalization when observations are not mapped.

Before you run the this routine, make sure you have the following:
-At least two regions and their respective population
-Date of NATIONWIDE policy implementation. (par.tp)

Left to implement in this routine:
3. Mapping all non-leading regions to C, to get asociated fig of pol_effects
4. 2^N -1 for distribution of policy effects
5. Smoothing with Chevichev Regression

% Left to implement:
-Bootstrapping figures
-fixing flips
-smoother comparison
-Automatize smoother
-Automatize for 1 region and search of reference region. 
-Automatize 2^N-1

%}

function SBI_norm(in_data,in_tp,in_stp)
%clear all
close all
clc
%% Set paths and "names"
path = cd;
par.outpath = [path,'/output'];
par.inpath  = [path,'/input'];

%load([par.inpath,'/example_data']); 

% Set labels:
par.outcome_name = 'XD';
par.time_name = 'Day';

% Set date of policy implementation: (must be in the same units as par.time)
%***************
par.tp = in_tp;
par.stp = in_stp;      % Step when the frequency is higher than 1
par.nm = 2;           % 2: Linear mapping 3: Quadratic
%***************

%% Data Treatment
% Extract Data
data.outcome = in_data.rate;
data.pop = in_data.popth;
data.names_st = ['REG01';'REG02';'REG03';'REG04';'REG05';'REG06';'REG07'];
par.time = in_data.time;       % This could be in Years, or Dates

par.ns = size(data.outcome,2);       % Number of regions 
par.nsp = par.ns+1;                  % Number of regions + (1) Benchmark
par.nt = size(data.outcome,1);       % Number of time obeservations

% Create Population Weights:
data.pop_weight = NaN(par.nt,par.ns);
for i = 1:par.nt      
   data.pop_weight(i,:) = data.pop(i,:)./sum(data.pop(i,:),2);    
end
%% Options for smoothing:
opt.smooth = 1;                % 0: No smoothing 1: Smoothing All(pre-post) 2: Smooth only-pre
opt.cheb_nodes = 1;            % (only for B-splines)1: use cheby nodes 2: Evenly Spaced 0: Automatic Matlab
opt.maw = 9;                   % Moving average window.
par.m = round(0.9*par.nt);     % Mumber of collocation points: min: n+1,
par.n = 7;                     % Degree of the polinomial;
opt.type_sth = 0;  
%{
 0: Interpolation smoother CSAPS 
 1: Chebyshev nodes with Cheby Regression
 2: B-Splines 
 3: Moving average
%}
%% Options for SBI Normalization:

opt.ilog = 1;           % 1.Multiplicative on the logs, 0.Multiplicative on the levels
opt.restr = 0;          % Restrict one level shifter (Deprecated)
opt.manual_guess = 2;   % 1: Pick Coefficient Initial Guess Manually 2: Polinomial analitical guess
opt.wght = 0;           % Weighted estimation
opt.imeth = 'linear';   % Type of interpolation
par.cut = 0;            % Trimming, keep it at 0
opt.sgues = 0;          % Use fmincon with sgues
par.min_match = 5;      % Minimum number of matched points
opt.lead_def = 1;       % Default leading region, from internal rankin 0: pick own or construct from more than one
opt.graph_smooth = 1;   % Graph smoothed functions
%% Boothstrapping params
opt.bb = 1;         % 1: Bootstrap every region
par.nb = 500;
par.wdth = 4;       % length of bootstrap interval in Block Bootstrapping

%% Options for solvers 
opt.sol = optimset('Display','off','TolFun',1.0e-08,'TolX',1.0e-08);
opt.sol_fmin = optimoptions('fmincon','Display','off','FunctionTolerance',1.0e-10,'OptimalityTolerance',1.0e-10);


% Generate some parameters
par.t0 = par.time(1);
par.t1 = par.time(end);
par.sp = 0.0008;%0.003;           % Smoothing parameter csaps
par.nm = 2;                       % 2: Linear mapping 3: Quadratic


% Parameters for figures
par.minyear = par.t0-1;     % Min year to graph
par.maxyear = par.t1;       % Max year to graph
%par.tu = par.tp-par.t0+1;   % Position of policy date
I = par.time<=par.tp;
par.tu = sum(I);   % Position of policy date
par.vis_ind = 1;
par.visible = 'on';        % 'on' if you want to see the graphs as they come
opt.title_coeff = 1;        % Add coefficients and flag to the title o fig1
opt.lw = 1.5;               % Linewidth
opt.lw2 = 1.7;
opt.fz2 = 19;% font size
opt.fz1 = 13;% font size
opt.fz3 = 23;% font size

% Initialize output vars

olpCT = NaN(par.nsp,1);
pol_estiCT = NaN(par.nsp,1);
LS_CT = NaN(par.nsp,1);
psi_estiCT = NaN(par.nsp,par.nm+1);
tp_normCT  = NaN(par.nsp,1);
fp_norm  = NaN(par.nsp,1);
eflagCT  = NaN(par.nsp,1);

olpTC = NaN(par.nsp,1);
pol_estiTC = NaN(par.nsp,1);
LS_TC = NaN(par.nsp,1);
psi_estiTC = NaN(par.nsp,par.nm+1);
tp_normTC  = NaN(par.nsp,1);
eflagTC  = NaN(par.nsp,1);


for i = 1:par.ns
    par.fapp = ['reg_',num2str(i),'sm_',num2str(opt.smooth),'tp_',num2str(opt.type_sth),'wght_',num2str(opt.wght)];
    %% Generate Treatment region Rest of XX 
    data.C = data.outcome(:,i);
    if i ==par.ns
               I = 1:1:(i-1);
    elseif  i ==1
               I = (i+1):1:par.ns;
    else
               I = [1:1:(i-1),(i+1):1:par.ns];
    end
    
    data.T = nansum(data.outcome(:,I).*(data.pop_weight(:,I)./sum(data.pop_weight(:,I),2)),2);
    par.Tname = char('RoXX');
    par.Cname = char(data.names_st(i,:));
    %% Smoothing
    [datas]=smooth_ts(data,par,opt);
    %% Perform SBI
    data.T = datas.T;
    data.C = datas.C;
    data.TO = datas.TO;
    data.CO = datas.CO;
    
    [out_r] = norm_func(data,par,opt);
    
    %% Extract and Save Results

    olpCT(i,1) = out_r.olpCT;
    olpTC(i,1) = out_r.olpTC;
    pol_estiCT(i,1) = out_r.pol_estiCT;
    pol_estiTC(i,1) = out_r.pol_estiTC;    
    psi_estiCT(i,:) = out_r.psi_estiCT;       
    psi_estiTC(i,:) = out_r.psi_estiTC;
    LS_CT(i,1) = out_r.LS_CT;       
    LS_TC(i,1) = out_r.LS_TC;
    tp_normCT(i,1)  = out_r.tp_normCT;
    tp_normTC(i,1)  = out_r.tp_normTC;
    fp_norm(i,1)  = out_r.fp_norm;
    eflagCT(i,1)  = out_r.flagCT;
    eflagTC(i,1)  = out_r.flagTC;
    
  
  
    
end

% Identify the first four leading regions:
[B,index_s]=sortrows(olpCT(1:end-1,1));  
datanames_aux = cellstr(data.names_st);
names_sort = datanames_aux(index_s);
disp(['Leading Region 1:',char(names_sort(end))])
disp(['Leading Region 2:',char(names_sort(end-1))])
disp(['Leading Region 3:',char(names_sort(end-2))])
disp(['Leading Region 4:',char(names_sort(end-3))])

% Pick (Generate) benchmark Result:
% For example a top 3 or bottom 3:
    % Two options: 
    %{
    Default: Pick the leading region from above
    I = index_s(end);
    Non default Generate a leading region yourself
    %}
    if opt.lead_def==0
       I = [2,6,4];
       par.Tname = char('RoXX');
       par.Cname = char('Top XX');
    else
       I = index_s(end);
       par.Tname = char('RoXX');
       par.Cname = char('Top XX');
    end
    
    i = par.nsp; % Index for the benchmark 
    
    data.C = nansum(data.outcome(:,I).*(data.pop_weight(:,I)./sum(data.pop_weight(:,I),2)),2);
    aux2 = 1:par.ns;
    aux1 = ismember(aux2,I);
    
    II =  aux2(not(aux1));   
    data.T = nansum(data.outcome(:,II).*(data.pop_weight(:,II)./sum(data.pop_weight(:,II),2)),2);
    
    data.oriC = data.C;
    data.oriT = data.T;
    %% Smoothing
    [datas]=smooth_ts(data,par,opt);
    %% Perform SBI
    data.T = datas.T;
    data.C = datas.C;
    data.TO = datas.TO;
    data.CO = datas.CO;
    
        % Graph Smooth
    if opt.graph_smooth==1
        figure(1000)
        hold on
        plot(par.time,datas.T,'r-')
        plot(par.time,datas.C,'b-')
        plot(par.time,datas.TO,'rx')
        plot(par.time,datas.CO,'bo')
    end
    
    [out_r] = norm_func(data,par,opt);
    
    %% Extract and Save Results

    olpCT(i,1) = out_r.olpCT;
    olpTC(i,1) = out_r.olpTC;
    pol_estiCT(i,1) = out_r.pol_estiCT;
    pol_estiTC(i,1) = out_r.pol_estiTC;    
    psi_estiCT(i,:) = out_r.psi_estiCT;       
    psi_estiTC(i,:) = out_r.psi_estiTC;
    LS_CT(i,1) = out_r.LS_CT;       
    LS_TC(i,1) = out_r.LS_TC;
    tp_normCT(i,1)  = out_r.tp_normCT;
    tp_normTC(i,1)  = out_r.tp_normTC;
    fp_norm(i,1)  = out_r.fp_norm;
    eflagCT(i,1)  = out_r.flagCT;
    eflagTC(i,1)  = out_r.flagTC;
    
    out_r.Csmooth = datas.C;
    out_r.Tsmooth = datas.T;
    out_r.resid_C = zeros(par.nt,1); % for when there are zero deaths
    out_r.resid_T = zeros(par.nt,1); % for when there are zero deaths
    I =  datas.CO>0;
    J =  datas.C>0;
    II = logical(I.*J);
    out_r.resid_C(II) = log(datas.CO(II))-log(datas.C(II));
    I =  datas.TO>0;
    J =  datas.T>0;
    II = logical(I.*J);
    out_r.resid_T(II) = log(datas.TO(II))-log(datas.T(II));
    
  %% Perform Block Bootstrapping:
    flag_out = NaN(par.nb,1);           % Outlier indicator 1: outlier, 0: healthy iteration   
    phi_bd = -99999*ones(par.nb,3);           % Coefficients
    olp_bd = -99999*ones(par.nb,1);                  % Overlap interval
    pol_estiCT = -99999*ones(par.nb,1);         % Gamma
    LS_CT = -99999*ones(par.nb,1);                   % LS
    fvalCT = -99999*ones(par.nb,1);                 % Norm Min Error
    pdiff_CT = -99999*ones(par.nt,par.nb);      
    CO = -99999*ones(par.nt,par.nb);                 % C Series
    CO_norm = -99999*ones(par.nt,par.nb);       % C Series
    TO = -99999*ones(par.nt,par.nb);                 % T Series
    time_d = -99999*ones(par.nt,par.nb);           % Time
    tp_normCT = -99999*ones(par.nb,1);     % Norm policy date
    
    if opt.bb==1
    par.vis_ind = 0; 
    rng(123,'twister')
    for i = 1:par.nb
        disp(i)
        %{
        rng(123+i)
        bd = (randperm(par.nt))'; 
        cmeC_bd = cmeC(bd);
        cmeT_bd = cmeT(bd);
        %}

        rv_draw_pre = zeros(par.nt,1);
        rv_draw = zeros(par.nt,1);

        t0 = 1; % The region that stards further. 
        nr = par.tp-t0+1;  % Number of prepolicy days
        nx = ceil(nr/par.wdth); % Number of bootstrap intervals in prepolicy days ex 5
        ind = randi([0 nr-par.wdth],nx,1); % Generate nx=5 draws, including 0?
        tvec = [t0:par.tp]';
        for xc=1:nx % Roll on boothstrap intervals
            for tc=1:par.wdth % Roll on values inside an interval
                tcc = (xc-1)*par.wdth+tc+(t0-1); % 
                if ( tcc>par.tp )
                    break
                end
                rv_draw_pre(tcc) = tvec(ind(xc)+tc);
            end
        end
        % after policy reform:
        nr = par.nt-(par.tp+1)+1;
        nx = ceil(nr/par.wdth);
        ind = randi([0 nr-par.wdth],nx,1);
        tvec = [par.tp+1:par.nt]';
        for xc=1:nx
            for tc=1:par.wdth
                tcc = (xc-1)*par.wdth+tc+par.tp;
                if ( tcc>par.nt )
                    break
                end
                rv_draw(tcc) = tvec(ind(xc)+tc);
            end
        end

        % build error draw for both regions
        cmeC_bd = ones(par.nt,1);
        cmeT_bd = ones(par.nt,1);
        % Pre-policy
        cmeC_bd(t0:par.tp) =  out_r.resid_C(rv_draw_pre(t0:par.tp));
        cmeT_bd(t0:par.tp) =  out_r.resid_T(rv_draw_pre(t0:par.tp));
        % Post-Policy
        cmeC_bd(par.tp+1:par.nt) = out_r.resid_C(rv_draw(par.tp+1:par.nt));
        cmeT_bd(par.tp+1:par.nt) = out_r.resid_T(rv_draw(par.tp+1:par.nt));           

        data.C = out_r.Csmooth.*exp(cmeC_bd);%  Dashed Blue    
        data.T = out_r.Tsmooth.*exp(cmeT_bd);%  Dashed Red
        
        % Smooth 
            [datas]=smooth_ts(data,par,opt);
        %% Perform SBI
        data.T = datas.T;
        data.C = datas.C;
        data.TO = datas.TO;
        data.CO = datas.CO;


        [out] = norm_func(data,par,opt);

        

        % Save:
        flag_out(i) = out.flag_out;           % Coefficients
        if i ==100
            dert_stop = 1;
        end
        if flag_out(i)==0
            disp('++++++++++++++++++++++++++++++++++++++++++++++')
            disp('Iter Sucess')
            disp('++++++++++++++++++++++++++++++++++++++++++++++')

            phi_bd(i,:) = out.psi_estiCT;           % Coefficients
            olp_bd(i) = out.olpCT;                  % Overlap interval
            pol_estiCT(i) = out.pol_estiCT;         % Gamma
            LS_CT(i) = out.LS_CT;                   % LS
            fvalCT(i) = out.fvalCT;                 % Norm Min Error
            pdiff_CT(:,i) = out.pdiff_CT;  
            CO(:,i) = out.CO;                 % C Series
            CO_norm(:,i) = out.CO_norm;       % C Series
            TO(:,i) = out.TO;                 % T Series
            time_d(:,i) = out.time;           % Time
            tp_normCT(i) = out.tp_normCT;     % Norm policy date
        end
    end
    % Get CI's throught the percentile method

    save(['nb500 ',num2str(opt.type_sth)])
    
    % 
    %% Get the CI's
    lev = 0.1; % 90% 
    % Phi's
    [mean_phi,med_phi,UCI_phi,BCI_phi,ntest_phi,nobs_phi,Edist_phi]=getCI(phi_bd,lev,flag_out);
    % Policy Effects
    [mean_pe,med_pe,UCI_pe,BCI_pe,ntest_pe,nobs_pe,Edist_pe]=getCI(pol_estiCT,lev,flag_out);
    [mean_pet,med_pet,UCI_pet,BCI_pet,ntest_pet,nobs_pet,Edist_pet]=getCI(pol_esti_true,lev,flag_out);
    % Cum gamma
    [mean_cg,med_cg,UCI_cg,BCI_cg,ntest_cg,nobs_cg,Edist_cg]=getCI(pdiff_CT',lev,flag_out);
    [mean_cgt,med_cgt,UCI_cgt,BCI_cgt,ntest_cgt,nobs_cgt,Edist_cgt]=getCI(pdiff_True',lev,flag_out);
    % Series 
    [mean_cd,med_cd,UCI_cd,BCI_cd,ntest_cd,nobs_cd,Edist_cd]=getCI(CO',lev,flag_out);
    [mean_td,med_td,UCI_td,BCI_td,ntest_td,nobs_td,Edist_td]=getCI(TO',lev,flag_out);
    % policy date
    [mean_tpn,med_tpn,UCI_tpn,BCI_tpn,ntest_tpn,nobs_tpn,Edist_tpn]=getCI(tp_normCT,lev,flag_out);
    % Cut for day around the median Overlap interval
    flag_out_alt = flag_out;
    I = tp_normCT>med_tpn+0.5;
    II = tp_normCT<med_tpn-0.5;
    flag_out_alt(I) = 1;
    flag_out_alt(II) = 1;
    [~,med_peA,UCI_peA,BCI_peA,~,~,Edist_peA]=getCI(pol_estiCT,lev,flag_out_alt);
    [~,med_petA,~,~,~,~,~]=getCI(pol_esti_true,lev,flag_out_alt);
    % Do the graphs: 
    
    %% Get the Graphs
    dert_stop = 1;

    % Histogram
    g2 = figure(2);
    hold on
    histogram(Edist_pe)
    p1 = xline(med_pe,'r-');
    p2 = xline(med_pet,'b--');
    p3 = xline(UCI_pe,'k--');
    p4 = xline(BCI_pe,'k--');
    legend([p1,p2,p3],{'Estimated Median','True Median','90% CI'})
    title('Policy Effect ($\gamma$), All OLP','interpreter','latex')
    saveas(g2,['output/g2_', fapp, '.png']);
    
    % Cut on median overlap
    g1 = figure(1);
    hold on
    histogram(Edist_peA)
    p1 = xline(med_peA,'r-');
    p2 = xline(med_petA,'b--');
    p3 = xline(UCI_peA,'k--');
    p4 = xline(BCI_peA,'k--');
    legend([p1,p2,p3],{'Estimated Median','True Median','90% CI'})
    title('Policy Effect ($\gamma$) Median OLP','interpreter','latex')
    saveas(g1,['output/g1_', fapp, '.png']);
    
    % Cumulative gamma
    g3 = figure(3);
    hold on
    p1 = plot(1:par.nt,med_cg,'m-','linewidth',1.1);
    p2 = plot(1:par.nt,mean_cg,'m--','linewidth',1.1);
    p3 = plot(1:par.nt,med_cgt,'b--','linewidth',1.1);
    p4 = plot(1:par.nt,UCI_cg,'k--','linewidth',1.1);
    plot(1:par.nt,BCI_cg,'k--','linewidth',1.1)
    %ylim([-100,2])
    xlim([par.tp-1,med_tpn+1])
    yline(0,'-','Color',[0.5,0.5,0.5])
    xline(par.tp,'-','Color',[0.5,0.5,0.5])
    xline(med_tpn,'--','Color',[0.5,0.5,0.5])
    legend([p1,p2,p3,p4],{'Esti Median','Esti Mean','True Median','90% CI'})
    title('Policy Effect ($\gamma_{s}$)','interpreter','latex')
    xlabel('Stage')
    saveas(g3,['output/g3_', fapp, '.png']);

        ym = 215;
 
    % The Regions
    data.rate(:,1) = data.oriC;%  Dashed Blue    
    data.rate(:,2) = data.oriT;
    g4 = figure(4);
    hold on
    p1 = plot(1:par.nt,data.rate(:,1),'bo','linewidth',1.1);
    p2 = plot(1:par.nt,data.rate(:,2),'r^','linewidth',1.1);
    p3 = plot(1:par.nt,med_cd,'b-','linewidth',1.1);
    p4 = plot(1:par.nt,med_td,'r-','linewidth',1.1);
    p5 = plot(1:par.nt,UCI_cd,'b--','linewidth',1.1);
    p6 = plot(1:par.nt,BCI_cd,'b--','linewidth',1.1);
    p7 = plot(1:par.nt,UCI_td,'r--','linewidth',1.1);
    p8 = plot(1:par.nt,BCI_td,'r--','linewidth',1.1);
    xline(par.tp,'k-','linewidth',1.1);
    ylim([0,ym])
    xlim([10,80])
    legend([p3,p5],{'Median','90% CI'})
    title('Flow of Deaths','interpreter','latex')
    xlabel('Time')
    ylabel('XD')
    saveas(g4,['output/g4_', fapp, '.png']);
    end
    
 
    
    
  %% Save Results
  
  


  
    fname=['output/olpCT_', par.fapp, '.txt'];
    save(fname,'olpCT','-ascii','-double','-tabs');

    fname=['output/psi_estiCT_', par.fapp, '.txt'];
    save(fname,'psi_estiCT','-ascii','-double','-tabs');
    
    fname=['output/LS_CT_', par.fapp, '.txt'];
    save(fname,'LS_CT','-ascii','-double','-tabs');
    
    fname=['output/pol_estiCT_', par.fapp, '.txt'];
    save(fname,'pol_estiCT','-ascii','-double','-tabs');
    
    fname=['output/fp_norm', par.fapp, '.txt'];
    save(fname,'fp_norm','-ascii','-double','-tabs');
    
    fname=['output/tp_normCT_', par.fapp, '.txt'];
    save(fname,'tp_normCT','-ascii','-double','-tabs');
    
    fname=['output/eflagCT_', par.fapp, '.txt'];
    save(fname,'eflagCT','-ascii','-double','-tabs');
    
    fname=['output/olpTC_', par.fapp, '.txt'];
    save(fname,'olpTC','-ascii','-double','-tabs');

    fname=['output/psi_estiTC_', par.fapp, '.txt'];
    save(fname,'psi_estiTC','-ascii','-double','-tabs');
    
    fname=['output/LS_TC_', par.fapp, '.txt'];
    save(fname,'LS_TC','-ascii','-double','-tabs');
    
    fname=['output/pol_estiTC_', par.fapp, '.txt'];
    save(fname,'pol_estiTC','-ascii','-double','-tabs');
    
    fname=['output/tp_normTC_', par.fapp, '.txt'];
    save(fname,'tp_normTC','-ascii','-double','-tabs');
    
    fname=['output/eflagTC_', par.fapp, '.txt'];
    save(fname,'eflagTC','-ascii','-double','-tabs');
    

    
    
    %% Save Results


disp('End of Routine')
end
%-------------------------------------------------------------------------
function [datas]=smooth_ts(data,par,opt)
%{
Input:
    data.outcome: Raw Series
Output:
    data.soutcome: Smooth Outcome
%}
% Keep original data:
    datas.TO = data.T;
    datas.CO = data.C;
    datas.T = data.T;
    datas.C = data.C;

    if opt.smooth==1
        pos = par.nt;
    elseif opt.smooth==2
        pos = par.tu;
    end

    if opt.smooth~=0
       
       try
           if opt.type_sth==0     % Poynomial Smoother
               datas.T(1:pos) = csaps(par.time(1:pos),data.T(1:pos),par.sp,par.time(1:pos));
               datas.C(1:pos) = csaps(par.time(1:pos),data.C(1:pos),par.sp,par.time(1:pos));
           elseif opt.type_sth==1 %  Cheby Polynomials with Cheby nodes
               datas.T(1:pos) = coll(par.time(1:pos),data.T(1:pos),par,opt);
               datas.C(1:pos) = coll(par.time(1:pos),data.C(1:pos),par,opt);
           elseif opt.type_sth==2 % B-Splines 
               datas.T(1:pos) = spline_reg(par.time(1:pos),data.T(1:pos),par,opt);
               datas.C(1:pos) = spline_reg(par.time(1:pos),data.C(1:pos),par,opt);
           else                   % Moving average 
               datas.T(1:pos) = smoothdata(data.T(1:pos),'movmean',par.mov_w);
               datas.C(1:pos) = smoothdata(data.C(1:pos),'movmean',par.mov_w);
           end
       catch
           disp('Error While Smoothing');
           error('Break');
       end
    end
    end

  
%-------------------------------------------------------------------------
function [o_vals]=norm_func(data,par,opt)
close all
fapp = par.fapp;

%{
figure(100)
hold on
plot(data.CO,'bx')
plot(data.C,'b--')
plot(data.TO,'rx')
plot(data.T,'r--')
%}



% Keep full series
par.nt_long = length(data.C);
par.time_long = par.time;

% Trim to use prepolicy data for normalization
par.time = par.time(1:par.tu,1);
data.C = data.C(1:par.tu,1);

par.nt = length(data.C);

% Recompute some location parameters according to the trimming
data.T = data.T(1:par.tu);


par.tvec_c      = (1:par.nt)';
par.tvec_c_long = (1:par.nt_long)';



% Initial guess for mapping parameters
%{ 
phi(1): Magnitude Shifter
phi(2): Time Shifter
phi(3): Speed Shifter
%}
%% Mapping C to T
m_guess = [0.79,5.9 ,1.07];
%m_guess = [0.5,2 ,1];
% Activar cuando es a manubrio
%{
opt.CT = 1;
[cf.phiCT,x0] = make_x(data,opt,m_guess,par); 
func_df(x0,data.T,data.C,opt,par)
dert_stop = 1;
%}
opt.CT = 1;
[cf.phiCT,x0] = make_x(data,opt,m_guess,par); 
if opt.manual_guess==1
    if size(cf.phiCT,2)~=par.nm+1
        disp('ERROR: The number of manual normalization params, is not the correct one')
    end
end
if opt.sgues==1
    [x0,~,~] = fmincon(@(x)func_df(x,data.T,data.C,opt,par),x0,[],[],[],[],[0,-Inf,-Inf],[Inf,Inf,Inf],[],opt.sol_fmin);
end
[xCT,fvalCT,eflagCT] = fminunc(@func_df,x0,opt.sol,data.T,data.C,opt,par);

if eflagCT==0
    disp('WARNING! SOLVER C to T QUIT DUE TO MAX ITER')
end
cf.phiCT = xCT;
disp(['Coefficients in numerical mapping C to T: ', num2str(cf.phiCT)]);
if par.nm == 2
disp(['Coefficients in numerical mapping C to T: ', num2str([1/cf.phiCT(1),-cf.phiCT(2)/cf.phiCT(3),1/cf.phiCT(3)])]);
end

%func_df(xCT,data.T,data.C,opt,par)

aux1=func_stage(par.tu,cf.phiCT,0,par);
tp_normCT_1 = interp1(par.tvec_c_long,par.time_long,aux1); % Convert to time units
aux1=func_stage(par.tu,cf.phiCT,1,par);
tp_normCT_2 = interp1(par.tvec_c_long,par.time_long,aux1);
%% Mapping T to C
%m_guess = [1/0.79,-5.9/1.07 ,1/1.07];
m_guess = [1/xCT(1),-xCT(2)/xCT(3) ,1/xCT(3)];
opt.CT = 0;
[cf.phiTC,x0] = make_x(data,opt,m_guess,par); 
if opt.sgues ==1
   [x0,~,~] = fmincon(@(x)func_df(x,data.C,data.T,opt,par),x0,[],[],[],[],[0,-Inf,-Inf],[Inf,Inf,Inf],[],opt.sol_fmin);
end
[xTC,fvalTC,eflagTC] = fminunc(@func_df,x0,opt.sol,data.C,data.T,opt,par);
if eflagTC==0
    disp('WARNING! SOLVER T to C QUIT DUE TO MAX ITER')
end
cf.phiTC = xTC;
disp(['Coefficients in numerical mapping T to C: ', num2str(cf.phiTC)]);
if par.nm == 2
disp(['Coefficients in numerical mapping T to C: ', num2str([1/cf.phiTC(1),-cf.phiTC(2)/cf.phiTC(3),1/cf.phiTC(3)])]);
end

aux1=func_stage(par.tu,cf.phiTC,0,par);
tp_normTC_1 = interp1(par.tvec_c_long,par.time_long,aux1);
aux1=func_stage(par.tu,cf.phiTC,1,par);
tp_normTC_2 = interp1(par.tvec_c_long,par.time_long,aux1);


aux1=func_stage(1,cf.phiCT,0,par);
yr_init = interp1(par.tvec_c_long,par.time_long,aux1,'linear',NaN);
disp('Results from C to T:')
disp(['Year in ',par.Cname,': ', num2str(par.tp), ' Stage in ',par.Tname,': ', num2str(tp_normCT_1)]);
disp(['Year in ',par.Cname,': ', num2str(par.t0), ' Stage in ',par.Tname,': ', num2str(yr_init)]);
disp('Results from T to C:')
disp(['Year in ',par.Tname,': ', num2str(par.tp), ' Stage in ',par.Cname,': ', num2str(tp_normTC_1)]);
%disp(['Year in ',par.Tname,': ', num2str(par.t0), ' Stage in ',par.Tname,': ', num2str(yr_init)]);

%% Compute Relevant Output:
% Normalized series
svec = par.tvec_c_long;
tvec = func_stage(svec,cf.phiCT(:),0,par);    % corresponding time in ROECD (i.e., when will ROECD look like Germany)

C_norm = cf.phiCT(1).*interp1(tvec(1:par.tu),data.C,par.tvec_c_long,'linear',NaN);
CO_norm = cf.phiCT(1).*interp1(tvec,data.CO,par.tvec_c_long,'linear',NaN);

tvec2 = func_stage(svec,cf.phiCT(:),1,par);    
T_norm = 1/cf.phiCT(1).*interp1(tvec2(1:par.tu),data.T,par.tvec_c_long,'linear',NaN);
TO_norm = 1/cf.phiCT(1).*interp1(tvec2,data.TO,par.tvec_c_long,'linear',NaN);

tvec_TC = func_stage(svec,cf.phiTC(:),0,par);    % corresponding time in ROECD (i.e., when will ROECD look like Germany)
T_norm_TC = cf.phiTC(1).*interp1(tvec_TC(1:par.tu),data.T,par.tvec_c_long,'linear',NaN);
TO_norm_TC = cf.phiTC(1).*interp1(tvec_TC,data.TO,par.tvec_c_long,'linear',NaN);

tvec2_TC = func_stage(svec,cf.phiTC(:),1,par);   
C_norm_TC = (1/cf.phiTC(1)).*interp1(tvec2_TC(1:par.tu),data.C,par.tvec_c_long,'linear',NaN);
CO_norm_TC = (1/cf.phiTC(1)).*interp1(tvec2_TC,data.CO,par.tvec_c_long,'linear',NaN);

%
if opt.ilog==1
svec = par.tvec_c_long;
tvec = func_stage(svec,cf.phiCT(:),0,par);    % corresponding time in ROECD (i.e., when will ROECD look like Germany)

C_norm = exp(cf.phiCT(1).*log(interp1(tvec(1:par.tu),data.C,par.tvec_c_long,'linear',NaN)));
CO_norm = exp(cf.phiCT(1).*log(interp1(tvec,data.CO,par.tvec_c_long,'linear',NaN)));

tvec2 = func_stage(svec,cf.phiCT(:),1,par);    
T_norm = exp(1/cf.phiCT(1).*log(interp1(tvec2(1:par.tu),data.T,par.tvec_c_long,'linear',NaN)));
TO_norm = exp(1/cf.phiCT(1).*log(interp1(tvec2,data.TO,par.tvec_c_long,'linear',NaN)));

tvec_TC = func_stage(svec,cf.phiTC(:),0,par);    % corresponding time in ROECD (i.e., when will ROECD look like Germany)
T_norm_TC = exp(cf.phiTC(1).*log(interp1(tvec_TC(1:par.tu),data.T,par.tvec_c_long,'linear',NaN)));
TO_norm_TC = exp(cf.phiTC(1).*log(interp1(tvec_TC,data.TO,par.tvec_c_long,'linear',NaN)));

tvec2_TC = func_stage(svec,cf.phiTC(:),1,par);   
C_norm_TC = exp((1/cf.phiTC(1)).*log(interp1(tvec2_TC(1:par.tu),data.C,par.tvec_c_long,'linear',NaN)));
CO_norm_TC = exp((1/cf.phiTC(1)).*log(interp1(tvec2_TC,data.CO,par.tvec_c_long,'linear',NaN)));
end
%}
%

%{
figure(1)
hold on
plot(par.tvec_c_long,data.CO,'b-')
plot(par.tvec_c_long,CO_norm,'b--')
plot(par.tvec_c_long,data.TO,'r-')
%
figure(2)
hold on
plot(par.tvec_c_long,data.CO,'b-')
plot(par.tvec_c_long,TO_norm,'r--')
plot(par.tvec_c_long,data.TO,'r-')
%
figure(3)
hold on
plot(par.tvec_c_long,data.CO,'b-')
plot(par.tvec_c_long,CO_norm_TC,'b--')
plot(par.tvec_c_long,data.TO,'r-')
%
figure(4)
hold on
plot(par.tvec_c_long,data.CO,'b-')
plot(par.tvec_c_long,TO_norm_TC,'r--')
plot(par.tvec_c_long,data.TO,'r-')
%}

% Policy Effect CT
In = sum(isnan(CO_norm([par.tu-1,par.tu,par.tu+1])));
% Same as gamma
pdiff_CT = NaN(size(data.TO));
tp_normCT_3=func_stage(par.tu,cf.phiCT,0,par);
if tp_normCT_3>size(data.TO,1)
    % Extreme case
    disp('ERROR: Detected extreme case CT')
    tp_normCT_3 = size(data.TO,1);
end

if In==0

  %% Estimated Policy Effect
      if opt.ilog ==1
          auxCO_norm = exp(cf.phiCT(1).*log(data.CO));
      else
          auxCO_norm = cf.phiCT(1).*data.CO;
      end
     Tpol = interp1(par.tvec_c_long,data.TO,tp_normCT_3,'linear');
     ldr =  Tpol*(tp_normCT_3-floor(tp_normCT_3));% Last day red
     ldb =  auxCO_norm(par.tu)*(tp_normCT_3-floor(tp_normCT_3));% Last day blue

     for ij = par.tu:floor(tp_normCT_3)         
        unr=sum(data.TO(par.tu:ij)); % Area under Blue         
        unb=sum(CO_norm(par.tu:ij));
        pdiff_CT(ij) = (unr-unb)/unb;
     end
     unr=sum(data.TO(par.tu:floor(tp_normCT_3)))+ldr; % Area under Blue
     unb=sum(CO_norm(par.tu:floor(tp_normCT_3)))+ldb;
     gammarCT = (unr-unb)/unb;
     olp_CT = tp_normCT_1-par.tp;
     LS_CT = unr-unb;
     DLS_CT = LS_CT/olp_CT;
      % Normalize fvalCT
     
     aux_Cnorm = cf.phiCT(1).*interp1(tvec(1:par.tu),data.C(1:par.tu),svec(1:par.tu),'linear',NaN);     
     avT = mean(data.TO(1:par.tu)); % Average deaths in reference region
     I = not(isnan(aux_Cnorm));   
     numpoCT = sum(I);  % Number of matched points
     fvalCT =  fvalCT./(avT.*numpoCT); % Normalize it:

     flag_out = 0;
     % Fish for outliers
     if isnan(gammarCT)
        flag_out = 1; % wrong sign
     end
     if olp_CT<0
        flag_out = 1;  % Flip
     end      
     if LS_CT>=0
        flag_out = 1; % Flip
     end
     if eflagCT<=0
       flag_out = 1; % Non Convergence, or corner solution
     end  
else
    %I = not(isnan(ln_gdpn_rest_sc));

     gammarCT = NaN;
     olp_CT = NaN;
     LS_CT = NaN;
     DLS_CT = NaN;    
    
end

% Policy Effect CT
In = sum(isnan(TO_norm_TC([par.tu-1,par.tu,par.tu+1])));
% Same as gamma

pdiff_TC = NaN(size(data.TO));
tp_normTC_3=func_stage(par.tu,cf.phiTC,0,par);
if tp_normTC_3<1
    % Extreme case
    disp('ERROR: Detected extreme case TC')
    tp_normTC_3 = 1;
end
if In==0

    %% Estimated Policy Effect
      if opt.ilog==1
          auxTO_norm = exp(cf.phiTC(1).*log(data.TO));
      else
          auxTO_norm = cf.phiTC(1).*data.TO;
      end
     Cpol = interp1(par.tvec_c_long,data.CO,tp_normTC_3,'linear');
     ldr =  auxTO_norm(par.tu)*(ceil(tp_normTC_3)-tp_normTC_3);% Last day red
     ldb =  Cpol*(ceil(tp_normTC_3)-tp_normTC_3);% Last day blue

     for ij = ceil(tp_normTC_3):par.tu
        unr=sum(TO_norm_TC(ceil(tp_normTC_3):ij))+ldr; % Area under Blue
        unb=sum(data.CO(ceil(tp_normTC_3):ij))+ldb;
        pdiff_TC(ij) = (unr-unb)/unb;
     end
     unr=sum(TO_norm_TC(ceil(tp_normTC_3):par.tu))+ldr; % Area under Blue
     unb=sum(data.CO(ceil(tp_normTC_3):par.tu))+ldb;
     gammarTC = (unr-unb)/unb;
     olp_TC = par.tp-tp_normTC_1;
     LS_TC = unr-unb;
     DLS_TC = LS_TC/olp_TC;
     
         aux_Tnorm = cf.phiTC(1).*interp1(tvec_TC(1:par.tu),data.T(1:par.tu),svec(1:par.tu),'linear',NaN);     
     avC = mean(data.CO(1:par.tu)); % Average deaths in reference region
     I = not(isnan(aux_Tnorm));   
     numpoTC = sum(I);  % Number of matched points    
     fvalTC =  fvalTC./(avC.*numpoTC); % Normalize it:
     

else
    %I = not(isnan(ln_gdpn_rest_sc));

     gammarTC = NaN;
     olp_TC = NaN;
     LS_TC = NaN;
     DLS_TC = NaN;    
    
end

% Display Policy effects:
disp(['Policy Effect C to T: ',num2str(round(gammarCT,3))])
disp(['Policy Effect T to C: ',num2str(round(gammarTC,3))])
disp(['fvalCT  : ',num2str(fvalCT)])
disp(['fvalTC  : ',num2str(fvalTC)])
% Create Flag:
[~,pos11] = min(abs(par.tvec_c-par.tp));
[~,pos21] = min(abs(par.tvec_c-floor(tp_normCT_1)));
[~,pos12] = min(abs(par.tvec_c-ceil(tp_normCT_2)));
[~,pos22] = min(abs(par.tvec_c-par.tp));

[~,pos11TC] = min(abs(par.tvec_c-par.tp));
[~,pos21TC] = min(abs(par.tvec_c-ceil(tp_normTC_1)));
[~,pos12TC] = min(abs(par.tvec_c-floor(tp_normTC_2)));
[~,pos22TC] = min(abs(par.tvec_c-par.tp));

I = pos21-pos11<0;    
II = eflagCT==-1||eflagCT==-3;
III = isnan(tp_normCT_1);
o_vals.flagCT = 0;
if I ==1 || III==1
    o_vals.flagCT = 2;
end
if II ==1
    o_vals.flagCT = 1;
end


if I==1 || II==1 || III==1
    disp(['Mapping CT Not Working for region: ',par.Cname,' Flag: ',num2str(o_vals.flagCT)])
    o_vals.olpCT = -999;
    o_vals.pol_estiCT = -999;
    o_vals.LS_CT = -999;  
    o_vals.psi_estiCT = [-999,-999,-999];
    o_vals.tp_normCT  = -999;
    o_vals.fp_norm  = -999;
    o_vals.fvalCT =  -999;
    o_vals.flag_out = 1;
    o_vals.pdiff_CT = ones(size(pdiff_CT)).*(-999);
    
else     

    o_vals.olpCT = olp_CT;
    o_vals.pol_estiCT = gammarCT;  
    o_vals.LS_CT = LS_CT;  
    o_vals.psi_estiCT = cf.phiCT;
    o_vals.tp_normCT  = tp_normCT_1;
    o_vals.fp_norm  = yr_init;
    o_vals.fvalCT =  fvalCT;
    o_vals.flag_out = flag_out;
    o_vals.pdiff_CT = pdiff_CT;
    
    o_vals.TO = data.TO;
    o_vals.CO = data.CO;  
    o_vals.CO_norm = CO_norm;
    o_vals.time = par.time_long;
end


I = pos11TC-pos21TC<0;    
II = eflagTC==-1||eflagTC==-3;
III = isnan(tp_normTC_1);
o_vals.flagTC = 0;
if I ==1 || III==1
    o_vals.flagTC = 2;
end
if II ==1
    o_vals.flagTC = 1;
end

if I==1 || II==1 || III==1
    disp(['Mapping TC Not Working for region: ',par.Cname,' Flag: ',num2str(o_vals.flag)])
    o_vals.olpTC = -999;
    o_vals.pol_estiTC = -999;
     o_vals.LS_TC = -999;  
    o_vals.psi_estiTC = [-999,-999,-999];
    o_vals.tp_normTC  = -999;
    o_vals.pdiff_TC = ones(size(pdiff_TC)).*(-999);

    
else     


    o_vals.olpTC = olp_TC;%par.tp-tp_normTC_1;
    o_vals.pol_estiTC = gammarTC;
    o_vals.LS_TC = LS_TC;  
    o_vals.psi_estiTC = cf.phiTC;
    o_vals.tp_normTC  = tp_normTC_1;
    o_vals.fvalTC =  fvalTC;
    o_vals.pdiff_TC = pdiff_TC;

end



%% Generate and save Figures:
%
if par.vis_ind==1
if isnan(tp_normTC_1)
    tp_normTC_1 = par.t0;
    tp_normTC_3 = par.t0;
end

if isnan(tp_normCT_1)
    tp_normCT_1 = par.t1;
    tp_normCT_3 = par.t1;
end



% Generate title:
aosCT = sprintf('%.2f,' , [cf.phiCT,o_vals.flagCT]);
aosCT = aosCT(1:end-1);  % strip final comma

aosTC = sprintf('%.2f,' , [cf.phiTC,o_vals.flagTC]);
aosTC = aosTC(1:end-1);  % strip final comma

%% Mapping Function
f1=figure(1);
set(f1, 'Visible',par.visible);
hold on
plot(svec,tvec,'k-','linewidth',1.1)
plot(svec,svec,'k--','linewidth',1)
legend(par.Cname,'Location','Best')
xlabel('Time ')
ylabel('Stage')
title(aosCT);
print('-depsc',['output/f1_', fapp, '.eps']);
saveas(f1,['output/f1_', fapp,'.png']);
hold off

%% Before Normalization
ymin = min([min(data.TO),min(data.CO),min(CO_norm)]);
ymax = max([max(data.TO),max(data.CO),max(CO_norm)])*1.1;
ymax2 = max([max(data.TO),max(CO_norm)])*1.1;

f2=figure(2);
set(f2, 'Visible',par.visible);
hold on
plot(par.time_long(1:end),data.TO,'r^','linewidth',opt.lw)
plot(par.time_long,data.CO,'bo','linewidth',opt.lw)
plot(par.time_long(1:par.tu),data.T(1:par.tu),'r-','linewidth',opt.lw)
plot(par.time_long(1:par.tu),data.C(1:par.tu),'b-','linewidth',opt.lw)
xline(par.tp,'k-','linewidth',1.1)
legend(par.Tname,par.Cname,'Location','Best','interpreter','latex','FontSize',opt.fz1,'box','off')
xlim([par.t0,par.t1]);
ylim([ymin,ymax]);
set(gca, 'FontSize',opt.fz1);
ylabel(par.outcome_name ,'interpreter','latex','FontSize',opt.fz2)
xlabel(par.time_name,'interpreter','latex','FontSize',opt.fz2)
print('-depsc',['output/f2_', fapp, '.eps']);
saveas(f2,['output/f2_', fapp,'.png']);
hold off

%% After Normalization C to T 
f3=figure(3);
set(f3, 'Visible',par.visible);
hold on
plot(par.time_long(1:end),data.TO,'r^','linewidth',opt.lw)
plot(par.time_long,CO_norm,'bo','linewidth',opt.lw)
plot(par.time_long(1:par.tu),C_norm(1:par.tu),'bx--','linewidth',opt.lw)
plot(par.time_long(1:par.tu),data.T(1:par.tu),'r-','linewidth',opt.lw)
area([par.tp,tp_normCT_1],[9999 9999],-9999,'FaceColor','r','FaceAlpha',.1,'EdgeColor','None','linestyle','none');
xline(par.tp,'k-','linewidth',1.1)
xline(tp_normCT_1,'k--','linewidth',1.1)
legend(par.Tname, [par.Cname,' Normalized'],'Location','Best','interpreter','latex','FontSize',opt.fz1,'box','off')
xlim([par.t0,par.t1]);
ylim([ymin,ymax2]);
set(gca, 'FontSize',opt.fz1);
ylabel(par.outcome_name,'interpreter','latex','FontSize',opt.fz2)
xlabel('Stage','interpreter','latex','FontSize',opt.fz2)
print('-depsc',['output/f3_', fapp, '.eps']);
saveas(f3,['output/f3_', fapp, '.png']);
hold off

% Zoom 
if par.tp-tp_normCT_1>0
% A flip
    
f4=figure(4);
set(f4, 'Visible',par.visible);
hold on
plot(par.time_long(1:end),data.TO,'r^','MarkerSize',8,'linewidth',opt.lw)
plot(par.time_long,CO_norm,'bo','MarkerSize',8,'linewidth',opt.lw)
xline(par.tp,'k-','linewidth',1.1)
xline(tp_normCT_1,'k--','linewidth',1.1)
legend(par.Tname, [par.Cname,' Normalized'],'Location','Best','interpreter','latex','FontSize',opt.fz1,'box','off')
set(gca, 'FontSize',opt.fz1);
ylabel(par.outcome_name ,'interpreter','latex','FontSize',opt.fz2)
xlabel('Stage','interpreter','latex','FontSize',opt.fz2)
print('-depsc',['output/f4_', fapp, '.eps']);
saveas(f4,['output/f4_', fapp, '.png']);
hold off    
    
f5 = figure(5);
set(f5, 'Visible',par.visible);
hold on
plot(par.time_long,pdiff_CT,'ko','MarkerFaceColor','k','MarkerSize',10)
plot(par.time_long,zeros(par.nt_long,1),'Color', [0 0 0]+0.2)
xline(par.tp,'k-','linewidth',1.1)
xline(tp_normCT_1,'k--','linewidth',1.1)
legend('$\gamma$','Location','northwest','FontSize',opt.fz1,'interpreter','latex','box','off')
set(gca, 'FontSize',opt.fz1);
ylabel('$\gamma(s)$','interpreter','latex','FontSize',opt.fz2)
xlabel('Stage','interpreter','latex','FontSize',opt.fz2)
print('-depsc',['output/f5_', fapp, '.eps']);
saveas(f5,['output/f5_', fapp, '.png']);
hold off  
    

else
ymin = min([min(data.TO(par.tu-1:ceil(tp_normCT_3)+1)),min(CO_norm(par.tu-1:ceil(tp_normCT_3)+1))])*0.95;
ymax = max([max(data.TO(par.tu-1:ceil(tp_normCT_3)+1)),max(CO_norm(par.tu-1:ceil(tp_normCT_3)+1))])*1.05;
   
f4=figure(4);
set(f4, 'Visible',par.visible);
hold on
plot(par.time_long(1:end),data.TO,'r^','MarkerSize',8,'linewidth',opt.lw2)
plot(par.time_long,CO_norm,'bo','MarkerSize',8,'linewidth',opt.lw2)
area([par.tp,tp_normCT_1],[9999 9999],-9999,'FaceColor','r','FaceAlpha',.1,'EdgeColor','None','linestyle','none');
xline(par.tp,'k-','linewidth',1.1)
xline(tp_normCT_1,'k--','linewidth',1.1)
legend(par.Tname, [par.Cname,' Normalized'],'Location','Best','interpreter','latex','FontSize',opt.fz1,'box','off')
xlim([par.tp-1,tp_normCT_1+1]);
ylim([ymin,ymax])
set(gca, 'FontSize',opt.fz1);
ylabel(par.outcome_name ,'interpreter','latex','FontSize',opt.fz2)
xlabel('Stage','interpreter','latex','FontSize',opt.fz2)
print('-depsc',['output/f4_', fapp, '.eps']);
saveas(f4,['output/f4_', fapp, '.png']);
hold off

%% Policy Effect Zoom
alt_yr = [par.tp,tp_normCT_1];
sizealt=size(alt_yr,2);

if o_vals.pol_estiCT<0
ymin2 = min([min(pdiff_CT),o_vals.pol_estiCT])*1.05;
ymax2 = +eps;
    if isnan(ymin2)  
        ymin2 = -0.1; 
    end
else
ymin2 = -eps;
ymax2 = max([max(pdiff_CT),o_vals.pol_estiCT])*1.05;
    if isnan(ymax2)  
        ymax2 = 0.1; 
    end
end
    

f5 = figure(5);
set(f5, 'Visible',par.visible);
hold on
plot(alt_yr,repmat(o_vals.pol_estiCT,[1,sizealt]),'m-','linewidth',opt.lw+1)
plot(par.time_long,pdiff_CT,'ko','MarkerFaceColor','k','MarkerSize',10)
plot(tp_normCT_1,gammarCT,'ko','MarkerFaceColor','k','MarkerSize',10)
plot(par.time_long,zeros(par.nt_long,1),'Color', [0 0 0]+0.2)
area([par.tp,tp_normCT_1],[9999 9999],-9999,'FaceColor','r','FaceAlpha',.1,'EdgeColor','None','linestyle','none');
xline(par.tp,'k-','linewidth',1.1)
xline(tp_normCT_1,'k--','linewidth',1.1)
%xlim([par.minyear,par.maxyear])
xlim([par.tp-1,tp_normCT_1+1]);
ylim([ymin2,ymax2])
legend('$\gamma$','Location','northwest','FontSize',opt.fz1,'interpreter','latex','box','off')
set(gca, 'FontSize',opt.fz1);
ylabel('$\gamma(s)$','interpreter','latex','FontSize',opt.fz2)
xlabel('Stage','interpreter','latex','FontSize',opt.fz2)
print('-depsc',['output/f5_', fapp, '.eps']);
saveas(f5,['output/f5_', fapp, '.png']);
hold off

f9 = figure(9);
set(f9, 'Visible',par.visible);
subplot(2,1,1)
hold on
plot(par.time_long(1:end),data.TO,'r^','MarkerSize',8,'linewidth',opt.lw2)
plot(par.time_long,CO_norm,'bo','MarkerSize',8,'linewidth',opt.lw2)
area([par.tp,tp_normCT_1],[9999 9999],-9999,'FaceColor','r','FaceAlpha',.1,'EdgeColor','None','linestyle','none');
xline(par.tp,'k-','linewidth',1.1)
xline(tp_normCT_1,'k--','linewidth',1.1)
%legend(par.Tname, [par.Cname,' Normalized'],'Location','Best','interpreter','latex','FontSize',opt.fz1,'box','off')
xlim([par.tp-1,tp_normCT_1+1]);
ylim([ymin,ymax])
set(gca, 'FontSize',opt.fz1);
ylabel(par.outcome_name ,'interpreter','latex','FontSize',opt.fz2)
xlabel('Stage','interpreter','latex','FontSize',opt.fz2)
subplot(2,1,2)
hold on
plot(alt_yr,repmat(o_vals.pol_estiCT,[1,sizealt]),'m-','linewidth',opt.lw+1)
plot(par.time_long,pdiff_CT,'ko','MarkerFaceColor','k','MarkerSize',10)
plot(tp_normCT_1,gammarCT,'ko','MarkerFaceColor','k','MarkerSize',10)
plot(par.time_long,zeros(par.nt_long,1),'Color', [0 0 0]+0.2)
area([par.tp,tp_normCT_1],[9999 9999],-9999,'FaceColor','r','FaceAlpha',.1,'EdgeColor','None','linestyle','none');
xline(par.tp,'k-','linewidth',1.1)
xline(tp_normCT_1,'k--','linewidth',1.1)
%xlim([par.minyear,par.maxyear])
xlim([par.tp-1,tp_normCT_1+1]);
ylim([ymin2,ymax2])
legend('$\gamma$','Location','northwest','FontSize',opt.fz1,'interpreter','latex','box','off')
set(gca, 'FontSize',opt.fz1);
ylabel('$\gamma(s)$','interpreter','latex','FontSize',opt.fz2)
xlabel('Stage','interpreter','latex','FontSize',opt.fz2)
print('-depsc',['output/f9_', fapp, '.eps']);
saveas(f9,['output/f9_', fapp, '.png']);
hold off  


end

%% For mapping T to C

%% After Normalization T to C 
aux_val = floor(tp_normTC_3);
if tp_normTC_1-par.tp>0
    % its a flip
    aux_val = par.tu;
end
ymin = min([min(data.TO),min(data.CO),min(CO_norm)]);
ymax3 = max([max(TO_norm_TC),max(data.CO)])*1.1;

f6=figure(6);
set(f6, 'Visible',par.visible);
hold on
plot(par.time_long,data.CO,'bo','linewidth',opt.lw)
plot(par.time_long,TO_norm_TC,'r^','linewidth',opt.lw)
plot(par.time_long(1:aux_val),T_norm_TC(1:aux_val),'rx--','linewidth',opt.lw)
plot(par.time_long(1:aux_val),data.C(1:aux_val),'b-','linewidth',opt.lw)
area([tp_normTC_1,par.tp],[9999 9999],-9999,'FaceColor','r','FaceAlpha',.1,'EdgeColor','None','linestyle','none');
xline(par.tp,'k-','linewidth',1.1)
xline(tp_normTC_1,'k--','linewidth',1.1)
legend(par.Cname, [par.Tname,' Normalized'],'Location','Best','interpreter','latex','FontSize',opt.fz1,'box','off')
xlim([par.t0,par.t1]);
ylim([ymin,ymax3]);
set(gca, 'FontSize',opt.fz1);
ylabel(par.outcome_name,'interpreter','latex','FontSize',opt.fz2)
xlabel('Stage','interpreter','latex','FontSize',opt.fz2)
print('-depsc',['output/f6_', fapp, '.eps']);
saveas(f6,['output/f6_', fapp, '.png']);
hold off

% Zoom 
if tp_normTC_1-par.tp>0
% A flip
    
f7=figure(7);
set(f7, 'Visible',par.visible);
hold on
plot(par.time_long,data.CO,'bo','MarkerSize',8,'linewidth',opt.lw)
plot(par.time_long,TO_norm_TC,'r^','MarkerSize',8,'linewidth',opt.lw)
xline(par.tp,'k-','linewidth',1.1)
xline(tp_normTC_1,'k--','linewidth',1.1)
legend(par.Cname, [par.Tname,' Normalized'],'Location','Best','interpreter','latex','FontSize',opt.fz1,'box','off')
set(gca, 'FontSize',opt.fz1);
ylabel(par.outcome_name ,'interpreter','latex','FontSize',opt.fz2)
xlabel('Stage','interpreter','latex','FontSize',opt.fz2)
print('-depsc',['output/f7_', fapp, '.eps']);
saveas(f7,['output/f7_', fapp, '.png']);
hold off    
    
f8 = figure(8);
set(f8, 'Visible',par.visible);
hold on
plot(par.time_long,pdiff_TC,'ko','MarkerFaceColor','k','MarkerSize',10)
plot(par.time_long,zeros(par.nt_long,1),'Color', [0 0 0]+0.2)
xline(par.tp,'k-','linewidth',1.1)
xline(tp_normTC_1,'k--','linewidth',1.1)
legend('$\gamma$','Location','northwest','FontSize',opt.fz1,'interpreter','latex','box','off')
set(gca, 'FontSize',opt.fz1);
ylabel('$\gamma(s)$','interpreter','latex','FontSize',opt.fz2)
xlabel('Stage','interpreter','latex','FontSize',opt.fz2)
print('-depsc',['output/f8_', fapp, '.eps']);
saveas(f8,['output/f8_', fapp, '.png']);
hold off 



    

else
ymin = min([min(TO_norm_TC(floor(tp_normTC_3)-1:par.tu+1)),min(data.CO(floor(tp_normTC_3)-1:par.tu+1))])*0.95;
ymax = max([max(TO_norm_TC(floor(tp_normTC_3)-1:par.tu+1)),max(data.CO(floor(tp_normTC_3)-1:par.tu+1))])*1.05;
   

f7=figure(7);
set(f7, 'Visible',par.visible);
hold on
plot(par.time_long,data.CO,'bo','MarkerSize',8,'linewidth',opt.lw2)
plot(par.time_long,TO_norm_TC,'r^','MarkerSize',8,'linewidth',opt.lw2)
area([tp_normTC_1,par.tp],[9999 9999],-9999,'FaceColor','r','FaceAlpha',.1,'EdgeColor','None','linestyle','none');
xline(par.tp,'k-','linewidth',1.1)
xline(tp_normTC_1,'k--','linewidth',1.1)
legend(par.Cname, [par.Tname,' Normalized'],'Location','Best','interpreter','latex','FontSize',opt.fz1,'box','off')
xlim([tp_normTC_1-1,par.tp+1]);
ylim([ymin,ymax])
set(gca, 'FontSize',opt.fz1);
ylabel(par.outcome_name ,'interpreter','latex','FontSize',opt.fz2)
xlabel('Stage','interpreter','latex','FontSize',opt.fz2)
print('-depsc',['output/f7_', fapp, '.eps']);
saveas(f7,['output/f7_', fapp, '.png']);
hold off

%% Policy Effect Zoom
alt_yr = [tp_normTC_1,par.tp];
sizealt=size(alt_yr,2);

if o_vals.pol_estiCT<0
ymin2 = min([min(pdiff_TC),o_vals.pol_estiTC])*1.05;
ymax2 = +eps;
    if isnan(ymin2)  
        ymin2 = -0.1; 
    end
else
ymin2 = -eps;
ymax2 = max([max(pdiff_TC),o_vals.pol_estiTC])*1.05;
    if isnan(ymax2)  
        ymax2 = 0.1; 
    end
end

f8 = figure(8);
set(f8, 'Visible',par.visible);
hold on
plot(alt_yr,repmat(o_vals.pol_estiTC,[1,sizealt]),'m-','linewidth',opt.lw+1)
plot(par.time_long,pdiff_TC,'ko','MarkerFaceColor','k','MarkerSize',10)
plot(tp_normTC_1,(ldr-ldb)/ldb,'ko','MarkerFaceColor','k','MarkerSize',10)
plot(par.time_long,zeros(par.nt_long,1),'Color', [0 0 0]+0.2)
area([tp_normTC_1,par.tp],[9999 9999],-9999,'FaceColor','r','FaceAlpha',.1,'EdgeColor','None','linestyle','none');
xline(par.tp,'k-','linewidth',1.1)
xline(tp_normTC_1,'k--','linewidth',1.1)
%xlim([par.minyear,par.maxyear])
xlim([tp_normTC_1-1,par.tp+1]);
ylim([ymin2,ymax2])
legend('$\gamma$','Location','northwest','FontSize',opt.fz1,'interpreter','latex','box','off')
set(gca, 'FontSize',opt.fz1);
ylabel('$\gamma(s)$','interpreter','latex','FontSize',opt.fz2)
xlabel('Stage','interpreter','latex','FontSize',opt.fz2)
print('-depsc',['output/f8_', fapp, '.eps']);
saveas(f8,['output/f8_', fapp, '.png']);
hold off  

f10 = figure(10);
set(f10, 'Visible',par.visible);
subplot(2,1,1)
hold on
plot(par.time_long,data.CO,'bo','MarkerSize',8,'linewidth',opt.lw2)
plot(par.time_long,TO_norm_TC,'r^','MarkerSize',8,'linewidth',opt.lw2)
area([tp_normTC_1,par.tp],[9999 9999],-9999,'FaceColor','r','FaceAlpha',.1,'EdgeColor','None','linestyle','none');
xline(par.tp,'k-','linewidth',1.1)
xline(tp_normTC_1,'k--','linewidth',1.1)
legend(par.Cname, [par.Tname,' Normalized'],'Location','Best','interpreter','latex','FontSize',opt.fz1,'box','off')
xlim([tp_normTC_1-1,par.tp+1]);
ylim([ymin,ymax])
set(gca, 'FontSize',opt.fz1);
ylabel(par.outcome_name ,'interpreter','latex','FontSize',opt.fz2)
xlabel('Stage','interpreter','latex','FontSize',opt.fz2)
subplot(2,1,2)
hold on
plot(alt_yr,repmat(o_vals.pol_estiTC,[1,sizealt]),'m-','linewidth',opt.lw+1)
plot(par.time_long,pdiff_TC,'ko','MarkerFaceColor','k','MarkerSize',10)
plot(tp_normTC_1,(ldr-ldb)/ldb,'ko','MarkerFaceColor','k','MarkerSize',10)
plot(par.time_long,zeros(par.nt_long,1),'Color', [0 0 0]+0.2)
area([tp_normTC_1,par.tp],[9999 9999],-9999,'FaceColor','r','FaceAlpha',.1,'EdgeColor','None','linestyle','none');
xline(par.tp,'k-','linewidth',1.1)
xline(tp_normTC_1,'k--','linewidth',1.1)
%xlim([par.minyear,par.maxyear])
xlim([tp_normTC_1-1,par.tp+1]);
ylim([ymin2,ymax2])
legend('$\gamma$','Location','northwest','FontSize',opt.fz1,'interpreter','latex','box','off')
set(gca, 'FontSize',opt.fz1);
ylabel('$\gamma(s)$','interpreter','latex','FontSize',opt.fz2)
xlabel('Stage','interpreter','latex','FontSize',opt.fz2)
print('-depsc',['output/f10_', fapp, '.eps']);
saveas(f10,['output/f10_', fapp, '.png']);
hold off   

end
end % visible
end
%---------------------------------------------------------------------------
function [pphi,x] = make_x(data,opt,m,par)
% Creates Coefficient guesses
    if opt.manual_guess==1
        pphi = m;
        x = pphi;
    elseif opt.manual_guess == 2
         pg = pol_guess(par,data);
        pphi(1,1) = pg(1); % Phi0
        pphi(1,2) = pg(2); % Phi1
        pphi(1,3) = pg(3); % Phi2
        x = pphi;
    else
      
        % Detect whether series is decreasing or increasing in the basis
        % of the last two points and T
        I = data.T(par.tu)-data.T(1)<0;
        if I ==1 % Debresing
            Cval = data.C(par.tu);
            [~,Tpos] = min(abs(data.T-Cval));% find the closest value in C
            ddl = abs(par.tu-Tpos);
        else % Increasing
            Tval = data.T(par.tu);
            [~,Cpos] = min(abs(data.C-Tval));% find the closest value in C
            ddl = abs(par.tu-Cpos);
        end
        if opt.CT==1
        
        pphi(1,1) = max(data.T)/max(data.C);
        pphi(1,2) = ddl;
        pphi(1,3) = 1;
        else
       
        pphi(1,1) = max(data.C)/max(data.T);
        pphi(1,2) = -ddl;
        pphi(1,3) = 1;
        end
        x = pphi;
    end

    
end
% -------------------------------------------------------------------------
function [fv,pphi,nt] = func_df(x,dta_T,dta_C,opt,par)
% Distance function
[dist,~,~,~,pphi] = func_dist(x,dta_T,dta_C,opt,par);

    nt = length(dist);
    if (opt.wght==1)
        wght = sqrt((1:nt)');
        wght(nt) = 100.0*wght(nt);
    else
        wght = ones(nt,1);
    end
    wdist = wght.*dist;
  
    fv = 0.5*(wdist')*wdist;    


end
% -------------------------------------------------------------------------
function [dist,pred_C_sc,svec,tvec,pphi] = func_dist(x,dta_T,dta_C,opt,par)
% Transform to logs
if opt.ilog==1
    dta_C = log(dta_C);
    dta_T = log(dta_T);
end

% get the coefficients
pphi = x;

% stage vector corresponding to time before unification
svec = (1:par.tu)';                   % time=stage in T
tvec = func_stage(svec,pphi(:),0,par);    % corresponding time in C

% Predicted scaled values in RO_XX with basic trimming according to C:
% This is the same either from C to T or T to C
if opt.ilog==1
    if par.cut==0
    I = dta_C>log(eps);
    II = dta_T>log(eps);
    dta_T = dta_T(II);
    aux_svec = svec(II); % save for the graph manubrio
    else
    I = dta_C>log(par.cut);
    II = dta_T>log(par.cut);
    dta_T = dta_T(II);
    aux_svec = svec(II); % save for the graph manubrio
    end
else
   
    I = dta_C>par.cut;
    II = dta_T>par.cut;
    dta_T = dta_T(II);
    aux_svec = svec(II); % save for the graph manubrio
  
end
try
    pred_C_sc = interp1(tvec(I),dta_C(I),svec(II),'linear',NaN);
    I = isnan(pred_C_sc);
    if sum(I)>1
        dert_stop =1; 
    end
catch
    dert_stop = 1;
    disp('Interpolation Error')
    pred_C_sc = NaN(size(par.tvec_c(II)));
end

% add the scaling factor
pred_C_sc = pphi(1).*pred_C_sc;


II = isnan(pred_C_sc);
I = not(isnan(pred_C_sc));
if sum(II)>0
    dert_stop = 1;
end
dist = dta_T(I) - pred_C_sc(I);

if isempty(dist)
    dist = 10e10;
end
% Compute Number of period matched:
if size(dist,1)<par.min_match
    disp('Not Enough matched points')
    dist = 1e10;
end

% Activar cuando es a manubrio
%{

figure(200)
hold on
plot(aux_svec,pred_C_sc,'b--x')
plot(aux_svec,dta_T,'r--')
plot(svec,dta_C,'b-')

%}

end

% -------------------------------------------------------------------------
function y = func_stage(x,pphi,back,par)
if back ==0 % in years of C
    y = pphi(2);
    for i = 3:(par.nm+1)
        y = y + (x.^(i-2)).*pphi(i);
    end
else % in years of T
    % Could Compute the inverse, but too comversome
    y = x./pphi(3)- pphi(2);
end

end

function [IDX] = rnd_PS(set,n)

% Draws a random sample from the power set without replacement
% If the number of samples equals to the powerset, then remove the empty
% set
% Learn to remove madrid.
rng(1234)  % set seed
nel = size(set',1);
y = randsample(2^nel,n); 
IDX = de2bi(y);
if size(IDX,2)==nel+1
    % Remove the zeros
    I = IDX(:,end)==1;
    IDX = IDX(:,1:end-1);
    IDX = IDX(not(I),:);

end
vars = 1; 
%function(x) v[intToBits(x) > 0])


end
%--------------------------------------------------------------------------
function cf = spline_reg(x_data,y_data,par,opt)

x_max = max(x_data);
x_min = min(x_data);
if opt.cheb_nodes==1

    [knots] = cheb_roots(par.m);
    x = ((sec(pi/(2*par.m))*knots+1)./2)*(x_max-x_min)+x_min;
    knots = x;
elseif opt.cheb_nodes ==2
    % Evenly spaced 
    x = linspace(x_min,x_max,par.m);
    knots = augknt(x,4);
else
    knots = par.m; % Automatic by matlab
end
% 4 to have 2 continous derivatives
cf_func = spap2(knots,par.n,x_data,y_data);
cf = fnval(cf_func,x_data);
dert_stop = 1;

end
%--------------------------------------------------------------------------
function [theta_guess,fx] = guess_theta(x_data,y_data,par,opt)
% To find a guess for theta I do the following: 
% I fit a polinomial to the data, with this I get continous values. 
% Then I use thos values to solve for analytical thetas which are my guess:
% 

x_max = max(x_data);
x_min = min(x_data);
zi = cheb_roots(par.m);
% Step 5: Convert nods to x interval x_max,x_min
knots = (sec(pi/(2*par.m))*zi+1)/2*(x_max-x_min)+x_min;
    
fx = csaps(x_data,y_data,par.sp,knots);
theta_guess = NaN(par.n,1);
val_aux = NaN(par.m,1);
% Here I compute the analitical thetas accordinf to Nakajima:

for j = 1:par.n
    if j ==1
        theta_guess(j) = 1/par.m*(sum(fx));
    else
        for i = 1:par.m
            %val_aux(i) = fx(i).*chebyshevT(j-1,zi(i));
            val_aux(i) = fx(i).*Tox(j-1,zi(i));
        end
        I = isnan(val_aux);
        if sum(I)>0
            disp('some error on polinomial')
        end
        theta_guess(j) = (2/par.m)*sum(val_aux);
    end
end


dert_stop = 1;
%{
figure(101)
hold on
plot(x_data,y_data,'k-o','markerfacecolor','k','linewidth',0.5)
plot(knots,fx,'b-','linewidth',1.1)
dert_Stop = 1;
%}
end
%--------------------------------------------------------------
function val = Tox(or,x)
%{
Manual Chabyshev polinomials, cuz its faster:
%}
    if or ==0
        val = 1;
    elseif or ==1
        val = x;
    elseif or ==2
        val = 2*(x^2)-1;
    elseif or ==3
        val = 4*(x^3)-3*(x);
    elseif or ==4
        val = 8*(x^4)-8*(x^2)+1;
    elseif or ==5
        val = 16*(x^5)-20*(x^3)+5*(x);
    elseif or ==6
        val = 32*(x^6)-48*(x^4)+18*(x^2)-1;
    elseif or ==7      
        val = 64*(x^7)-112*(x^5)+56*(x^3)-7*x;
    elseif or ==8      
        val = 128*(x^8)-256*(x^6)+160*(x^4)-32*(x^2)+1;
    end
end
%--------------------------------------------------------------------------
function [roots_vec] = cheb_roots(m)
%{
Computes the roots of the chebychev polinomial of order m
%}
i_aux = (1:m)';
roots_vec = -cos((2.*i_aux-1).*pi./(2*m));

end

%--------------------------------------------------------------------------
function [fv,chebyval] =  dist_cheby(theta,x_data,y_data,par)
% Evaluate the Cheby polinomial


% Step 1: 
x_max = max(x_data);
x_min = min(x_data);
% Step 2: Pick the order of the Chebyshev Polinomial
n = par.n;
nobs = par.nobs;
% Step 3: Pick the number of Collocation points
m = par.m;  %par.nobs;% Number of collocation points same as order or more

% Compute solution thetas, I need data counterpart to do this: 
% So I just evaluate the Polinomial 
% Convert data points x to -1,1 
zi = (1/sec(pi/(2*m))).*(((2*(x_data-x_min))/(x_max-x_min))-1);

cheby_poli = NaN(par.n,par.nobs);

for i = 1:n
    for ii = 1:nobs
        %cheby_poli(i,ii) = theta(i)*chebyshevT(i-1,zi(ii));
        cheby_poli(i,ii) = theta(i)*Tox(i-1,zi(ii));
    end
end

%z = ((sec(pi/(2*par.m))*roots_vec+1)./2)*(z_max-z_min)+z_min;
%z = (1/sec(pi/(2*m))).*(((2*(x_data-x_min))/(x_max-x_min))-1);

I = isnan(cheby_poli);
if sum(I,'all')>0
    disp('Something wrong is happening')
end
chebyval = sum(cheby_poli,1)';
% Get the difference: 
I = not(isnan(chebyval));
dist = y_data(I) - chebyval(I);

if isempty(dist)
    disp('Why is dist empty')
    dist = 10e10;
end

 wdist = dist;
  
 fv = 0.5*(wdist')*wdist;   

end
%--------------------------------------------------------------------
function [cf,theta_sol,fval]=coll(x_data,y_data,par,opt)
%{
x, y = data input, does not allow for missing values
n = degree of the polinomial n_pol-1 = order of the polinomial

%}


% Guess theta (Cheby Coeffs)

[theta_x0,~] = guess_theta(x_data,y_data,par,opt);


dert_stop = 1;
[theta_sol,fval,ef] =  fminunc(@(theta) dist_cheby(theta,x_data,y_data,par),theta_x0,opt.optmin);
% Evauate fit
%profile on
[~,cf] = dist_cheby(theta_sol,x_data,y_data,par);
%profile viewer

dert_Stop = 1;



end
%---------------------------------------------------------------------
%--------------------------------------------------------------------------
function pg = pol_guess(par,data)
% This fits a polinomial prepolicy and uses the analitical coefficients as
% guesses
% Chooses the guess with positive phi(2)

    pDC = polyfit(par.time(1:par.tu),data.C(1:par.tu),2);
    pDT = polyfit(par.time(1:par.tu),data.T(1:par.tu),2);
    auxDC = polyval(pDC,par.time(1:par.tu));
    auxDT = polyval(pDT,par.time(1:par.tu));
    % Normalize analically
    pa.A = (pDT(1)*pDC(2))/(pDC(1)*pDT(2));
    pa.B = (pDT(1)*2)/pDT(2);
    pa.C = pDT(2)/pDT(3);

    pa.c1 = pDC(2)*pa.A-pDC(3)*pa.C;
    pa.c2 = pDC(2)*pa.B+2*pDC(1)*pa.A-pDC(2)*pa.C;
    pa.c3 = 2*pDC(1)*pa.B-pDC(1)*pa.C;

    rvals = roots([pa.c3,pa.c2,pa.c1]);

    pa.phi1 = rvals(2); % choose this because it gives a positive phi2
    pa.phi2 = pa.A+pa.B*pa.phi1;
    pa.phi0 = (pa.phi2^2)*pDC(1)/pDT(1);

    if pa.phi2<0
        pa.phi1 = rvals(1); % choose this because it gives a positive phi2
        pa.phi2 = pa.A+pa.B*pa.phi1;
        pa.phi0 = (pa.phi2^2)*pDC(1)/pDT(1);
    end

    if pa.phi2<0
        pa.phi1 = 1; % choose this because it gives a positive phi2
        pa.phi2 = 2;
        pa.phi0 = 0.5;
        disp('Polinomial guess failed: using default')
    end
    pg(1) = 1/pa.phi0;
    pg(2) = -pa.phi1/pa.phi2;
    pg(3) = 1/pa.phi2;

    dert_stop = 1;
    %{
    tnorm = pg(2)+pg(3)*par.time(1:par.tu);
    figure(321)
    hold on
    plot(par.time(1:par.tu),auxDC,'b-')
    plot(par.time(1:par.tu),auxDT,'ro-')
    plot(tnorm,auxDC.*pg(1),'bx-')
    %}

    %{
    figure(301)
    hold on
    plot(par.time,data.C,'b--')
    plot(par.time,data.T,'r--')
    plot(par.time(1:par.tu),auxDC,'bx-')
    plot(par.time(1:par.tu),auxDT,'rx-')
    ylim([-50,250])
    xlim([0,80])
    %}

end
%
function [S_mean,S_median,UCI,BCI,ntest,nobs,S_Edist]=getCI(ser,lev,flag)
%{
This function computes the confidence bands for a given empirical
distribution using the percentile method.

This function also performs the loops
nc: number of Columns, 
nb: number of Boorstrap iterations
ser: series, should be NBxNC
flag: if 1 it is an outlier and should be taken out
nobs: number of healthy iterations-NaN's


%}
par.nc = size(ser,2);
%par.nb = nb;
S_Edist = [];              % Empirical distribution for histogram
S_median = NaN(par.nc,1);
S_mean   = NaN(par.nc,1);
UCI = NaN(par.nc,1);
BCI = NaN(par.nc,1);
CI_BU = NaN(par.nc,2);     % (B)Bottom (U)Upper
ntest = NaN(par.nc,1);
nobs = NaN(par.nc,1);

I = logical(flag);
ser = ser(not(I),:);

for rc = 1:par.nc
    sort_gam = sort(ser(:,rc));
    % exclude the NaNs, if only NaN's then skip
    I = isnan(sort_gam);
    if sum(I)==size(sort_gam,1)
        continue
    end
    sort_gam = sort_gam(not(I));
    S_Edist = sort_gam;
    nb_t = size(sort_gam,1);
    nobs(rc,1) = nb_t;   % Number of observations, this excludes the nans
    
    cut = round(nobs(rc,1)*lev/2); % cut for computing confidence bounds
    if cut ==0
        CI_BU(rc,:) = [-999,-999];
        UCI(rc,1) = -999;
        BCI(rc,1) = -999;
        S_median(rc,1) = -999;
        S_mean(rc,1) = -999;
    else
        CI_BU(rc,:) = [sort_gam(cut),sort_gam(nb_t-cut+1)];
        UCI(rc,1) = CI_BU(rc,2);
        BCI(rc,1) = CI_BU(rc,1);
        S_median(rc,1) = median(sort_gam);
        S_mean(rc,1) = mean(sort_gam);
        try
           % JB test with H0: Normal
            %[~,ntest(rc,1)] = jbtest(sort_gam);
            [~,p_jb] = jbtest(sort_gam);
            ntest(rc,1) = p_jb;%nb_t; 
        catch
        end
    end

end


end