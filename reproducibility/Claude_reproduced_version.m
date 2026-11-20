%% ah_onehit.m
% One-hit ethanol simulation of the flux-centric, coarse-grained hepatocyte
% model (Shoraka & Vadigepalli). Written from the manuscript equations only.
% Requires MATLAB R2016b+ (local functions in scripts, implicit expansion).
% written with OPUS 5.5 HIGH

clear; clc;
P = defaultParams();

% Representative individuals, Sec. 3.7 "Heterogeneous failure scenario" (Fig. 4)
%          label                    kd     y1    y2
cases = {'Healthy',               0.010, 0.20, 0.25
         'Diseased (ALD)',        0.005, 0.50, 0.60
         'Failed, compensated',   0.003, 0.75, 0.50
         'Failed, decompensated', 0.001, 0.75, 0.60};

res = cell(size(cases,1),1);
fprintf('%-24s %9s %9s %9s %11s %5s %5s  %s\n', ...
    'case','netUrea','netMET','netFABs','dGSH_T','OxF','MetD','regime');
for c = 1:size(cases,1)
    Pc = P;  Pc.kd = cases{c,2};  Pc.y1 = cases{c,3};  Pc.y2 = cases{c,4};
    tic; out = simulateOneHit(Pc); el = toc;
    res{c} = out;
    fprintf('%-24s %+9.3f %+9.3f %+9.3f %+11.2f %5d %5d  %s  (%.0f s)\n', ...
        cases{c,1}, out.netDev(3), out.netDev(2), out.netDev(4), out.GSHT, ...
        out.oxFail, out.metDecomp, out.regime, el);
end
plotResults(res, cases, P);

%% ======================================================================
function P = defaultParams()
% Modules: 1 GSH, 2 MET, 3 Urea, 4 FABs
P.modNames = {'GSH','MET','Urea','FABs'};
P.v0 = [2.0737 0.4002 0.6193 0.9068];        % mM/min, Sec. 3.6 (sums to N = 4)

% Adjacency A(i,k): resources of module k usable by module i (Sec. 3.4)
P.A = [0 7 4 2
       5 0 2 1
       6 5 0 4
       2 1 3 0];

% Metabolites: 1 GSH_mito, 2 SAMe_cyto, 3 ammonia_mito, 4 acetyl-CoA_mito, 5 malonyl-CoA_mito
P.metNames = {'[GSH]_{mito}','[SAMe]_{cyto}','[ammonia]_{mito}', ...
              '[a-CoA]_{mito}','[malonyl-CoA]_{mito}'};
P.class = [2 2 1 1 1];          % 1 = critically high (consumed), 2 = critically low (produced)
P.iGSH = 1;  P.iAcCoA = 4;

% Coupling matrix L (modules x metabolites). Supplement not available:
% unit stoichiometry per module-metabolite association (Sec. 3.4 list).
P.L = [1 0 0 0 0      % GSH  produces GSH
       0 1 0 0 0      % MET  produces SAMe
       0 0 1 1 0      % Urea consumes ammonia, acetyl-CoA
       0 0 0 1 1];    % FABs consumes acetyl-CoA, malonyl-CoA (Eq. 30)

P.m0   = [10 0.08 0 0 0];        % initial loads (mM)
P.mbar = [10 0.08 Inf Inf Inf];  % admissible upper bound (class 2) = set point
P.crit = [7.5 0.04 1 2 0.3];     % critical concentrations (Sec. 3.6)

% Allocation / kinetics (Table 1, Sec. 3.6)
P.Dt    = 60;      % characteristic clearance time (min)
P.beta  = 7;
P.theta = 5;
P.r     = 0.01;    % 1/min
P.eps   = 1e-12;
P.tau   = 1e4;     % mM^-2 min^-1
P.y3    = 1;       % acetate -> acetyl-CoA rate (mM/min)
P.kEtOH = 0.05;    % zero-order ethanol clearance (mM/min)

% Patient-specific (overwritten per case)
P.kd = 0.01;  P.y1 = 0.2;  P.y2 = 0.25;

% Perturbation and numerics (Sec. 3.7)
P.dose  = 17.4;           % mM (BAC 0.08 g/dL)
P.tPert = 200;            % min
P.Tend  = 200 + 4520;     % ~72 h after exposure
P.dt    = 0.005;          % Forward Euler step (min)
P.saveEvery = 1;          % min, recording resolution

% Failure criteria (Sec. 2.3)
P.deltaT = 100;           % terminal window (min)
P.deltaI = 0.1;           % decompensation threshold
P.epsOx  = 1e-9;          % "[H2O2] > 0" (Sec. 2.6); use 1 for the cytotoxic line in Fig. 4
end

%% ======================================================================
function out = simulateOneHit(P)
nS = numel(P.v0);  nM = numel(P.m0);
v0 = P.v0;  L = P.L;  dt = P.dt;

% ---- Precomputed steady-state structure --------------------------------
T     = v0 * L;                       % Eq. 3/9 : baseline throughput per metabolite
C     = (L .* v0.') ./ T;             % Eq. 4/10: capacity fractions
critS = P.crit .* L ./ sum(L,1);      % module share of critical concentration
Prow  = P.A ./ sum(P.A,2);            % Eq. 21  : routing p(i,k)
isC1  = (P.class == 1);

[pI,pJ] = find(L > 0);  pI = pI.';  pJ = pJ.';     % module-metabolite pairs
lin = sub2ind([nS nM], pI, pJ);
Cp = C(lin);  Sp = critS(lin);  Tp = T(pJ);
mbP = P.mbar(pJ);  crP = P.crit(pJ);
k1 = find(isC1(pJ));  k2 = find(~isC1(pJ));
nP = numel(lin);

% ---- State --------------------------------------------------------------
v = v0;  m = P.m0;  H = 0;  a = 0;  ac = 0;  dosed = false;

nSteps = round(P.Tend/dt);
stride = max(1, round(P.saveEvery/dt));
nRec   = floor(nSteps/stride) + 1;
Rt = zeros(nRec,1); Rv = zeros(nRec,nS); Rm = zeros(nRec,nM);
RH = zeros(nRec,1); Ra = zeros(nRec,1);
rec = 1;  Rv(1,:) = v;  Rm(1,:) = m;

netDev = zeros(1,nS);  termDev = zeros(1,nS);  Hmax = 0;
preSum = 0;  preN = 0;  postSum = 0;
tTerm = P.Tend - P.deltaT;

for n = 1:nSteps
    t = (n-1)*dt;

    % ---- 1. Demand, availability, burden per pair (Eqs. 5-7, 11-13) ----
    mp = m(pJ);  Cm = Cp .* mp;
    dP = zeros(1,nP);  aP = dP;  zP = dP;

    % Class 1: critically high
    dP(k1) = Cp(k1) .* (mp(k1)/P.Dt + Tp(k1));
    aP(k1) = max(0, 1 - Cm(k1)./Sp(k1));
    zP(k1) = Cm(k1) ./ Sp(k1);

    % Class 2: critically low
    for q = k2
        dP(q) = Cp(q) * ((mbP(q) - mp(q))/P.Dt + Tp(q));
        if Cm(q) > mbP(q)
            aP(q) = 1;  zP(q) = 0;
        elseif Cm(q) > Sp(q)
            aP(q) = (Cm(q) - Sp(q)) / (mbP(q) - Sp(q));
            zP(q) = (mbP(q) - Cm(q)) / (mbP(q) - crP(q));
        else
            aP(q) = 0;
            zP(q) = (mbP(q) - Cm(q)) / (mbP(q) - crP(q));
        end
        if mp(q) <= crP(q)
            zP(q) = min(P.theta, Sp(q)/max(Cm(q), realmin));
        end
    end

    % ---- 2. Aggregate to modules (Eqs. 14-17) -------------------------
    Dm = -Inf(nS,nM); Dm(lin) = dP; dS = max(Dm,[],2).';
    Am =  Inf(nS,nM); Am(lin) = aP; aS = min(Am,[],2).';
    Zm = -Inf(nS,nM); Zm(lin) = zP; zS = max(Zm,[],2).';

    % ---- 3. Resource reallocation (Eqs. 16, 18-23) --------------------
    s   = aS .* v;                                  % shareable pool
    w   = exp(P.beta*(zS - max(zS)));  w = w/sum(w); % softmax priority
    Rq  = w .* max(0, dS - v);                      % request intensity
    RP  = Rq.' .* Prow;                             % R_i * p_ik
    X   = RP ./ (sum(RP,1) + P.eps) .* s;           % X(i,k) = x_{i<-k}
    dv  = P.r * (sum(X,2).' - sum(X,1));            % conserves sum(v)

    % ---- 4. Slow metabolite loads (Eqs. 2, 8, 26, 28-30) --------------
    flow = v * L;
    vm = zeros(1,nM);
    vm(isC1)  = -flow(isC1)  + T(isC1)  - P.kd;     % k_p = baseline throughput
    vm(~isC1) =  flow(~isC1) - T(~isC1) + P.kd;     % k_c = baseline throughput
    mN = m;
    g = isC1 & (m > 0 | vm >= 0);
    mN(g) = m(g) + vm(g)*dt;
    g = ~isC1 & ((m > 0 & m <= P.mbar) | (m <= 0 & vm >= 0));
    mN(g) = m(g) + vm(g)*dt;
    g = ~isC1 & (m > P.mbar);
    mN(g) = m(g) - T(g)*dt;
    mN = max(mN, 0);

    % ---- 5. Ethanol chemistry (Sec. 3.4, Eqs. 25-28) ------------------
    if ~dosed && t >= P.tPert
        a = a + P.dose;  dosed = true;
    end
    burned = min(a, P.kEtOH*dt);  a = a - burned;     % zero-order clearance
    tr = min(ac, P.y3*dt);                            % acetate -> acetyl-CoA
    ac = ac - tr + P.y2*burned;
    mN(P.iAcCoA) = mN(P.iAcCoA) + tr;
    H = H + P.y1*burned;                              % H2O2 production
    G = mN(P.iGSH);
    xi = min([P.tau*G^2*H*dt, H, G/2]);               % 2GSH + H2O2, stop-gapped
    H = H - xi;
    mN(P.iGSH) = G - 2*xi;

    % ---- 6. Advance --------------------------------------------------
    m = mN;  v = v + dv*dt;
    tn = n*dt;

    netDev = netDev + (v - v0)*dt;
    if tn > tTerm, termDev = termDev + abs(v - v0)*dt; end
    Hmax = max(Hmax, H);
    if tn <= P.tPert
        preSum = preSum + m(P.iGSH);  preN = preN + 1;
    else
        postSum = postSum + m(P.iGSH)*dt;
    end

    if mod(n, stride) == 0
        rec = rec + 1;
        Rt(rec) = tn; Rv(rec,:) = v; Rm(rec,:) = m; RH(rec) = H; Ra(rec) = a;
    end
end

% ---- Outputs and classification (Sec. 2.3) -----------------------------
out.t = Rt(1:rec);  out.v = Rv(1:rec,:);  out.m = Rm(1:rec,:);
out.H2O2 = RH(1:rec);  out.EtOH = Ra(1:rec);
out.netDev  = netDev;                     % integral of (v - v0)
out.termDev = termDev;                    % Eq. 1 per module
out.GSHT    = postSum - (preSum/preN)*(P.Tend - P.tPert);  % integrated GSH deviation
out.Hmax    = Hmax;
out.oxFail    = Hmax > P.epsOx;
out.metDecomp = any(termDev >= P.deltaI);
labels = {'none (O^c & M^c)','metabolic only (O^c & M)'; ...
          'oxidative only (O & M^c)','joint (O & M)'};
out.regime = labels{out.oxFail+1, out.metDecomp+1};
end

%% ======================================================================
function plotResults(res, cases, P)
nC = numel(res);
figure('Color','w','Position',[80 80 320*nC 720]);
for c = 1:nC
    R = res{c};
    subplot(3,nC,c);
    plot(R.t, R.H2O2, 'Color',[0.93 0.69 0.13]); hold on;
    yline(1,'r','Cytotoxic');
    title(cases{c,1}); if c==1, ylabel('[H_2O_2]_{mito} (mM)'); end

    subplot(3,nC,nC+c);
    burden = abs(R.m - P.m0) ./ P.crit;          % Eq. 24
    plot(R.t, burden);
    if c==1, ylabel('Metabolic burden'); end
    if c==nC, legend(P.metNames,'Location','eastoutside'); end

    subplot(3,nC,2*nC+c);
    plot(R.t, R.v ./ P.v0);
    xlabel('Time (min)'); if c==1, ylabel('Normalized velocity'); end
    if c==nC, legend(P.modNames,'Location','eastoutside'); end
end
end