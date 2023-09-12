function [outstr, randstr] = granger_causality_Seth2015(X, fs, source_labels, nsamples, nrand, verbose, plotting)
% X = (m x n) matrix of (sources x time)
% fs = sampling rate
% plotting = true;
% verbose = true;
% Adapted from :
% % Barnett L, Seth AK. The MVGC multivariate Granger causality toolbox: 
% % a new approach to Granger-causal inference. J Neurosci Methods. 2014 Feb 15;
% % 223:50-68. doi: 10.1016/j.jneumeth.2013.10.018. Epub 2013 Nov 5. PMID: 24200508.
% https://www.mathworks.com/matlabcentral/fileexchange/78727-the-multivariate-granger-causality-mvgc-toolbox
Xog = X;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
regmode   = 'LWR';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 20;     % maximum model order for model order estimation

acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)

tstat     = 'F';    % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDRD'; % multiple hypothesis test correction (see routine 'significance')

% fs        = 4;    % sample rate (Hz)
fres      = [];     % frequency resolution (empty for automatic calculation)

seed      = 0;      % random seed (0 for unseeded)
rng_seed(seed);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nvars = size(Xog,1); % number of variables

if ~isempty(nsamples) % break up X into trials of length nsamples
    n_trials = floor(size(Xog,2)/nsamples);
    X = NaN(nvars, nsamples, n_trials);
    for i = 1:n_trials
        X(:,:,i) = Xog(:, nsamples*(i-1) + 1: nsamples*i);
    end
else
    X = Xog;
end

%%



if verbose; ptic('\n*** tsdata_to_infocrit\n'); end
[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode, verbose);
if verbose; ptoc('*** tsdata_to_infocrit took '); end

% Plot information criteria.
if plotting==true
figure(1); clf;
plot_tsdata([AIC BIC]',{'AIC','BIC'},1/fs);
title('Model order estimation');
end
% amo = size(AT,3); % actual model order

if verbose
    fprintf('\nbest model order (AIC) = %d\n',moAIC);
    fprintf('best model order (BIC) = %d\n',moBIC);
%     fprintf('actual model order     = %d\n',amo);
end
% Select model order.

% if strcmpi(morder,'actual')
%     morder = amo;
%     if verbose; fprintf('\nusing actual model order = %d\n',morder); end
% elseif
if strcmpi(morder,'AIC')
    morder = moAIC;
    if verbose; fprintf('\nusing AIC best model order = %d\n',morder); end
elseif strcmpi(morder,'BIC')
    morder = moBIC;
    if verbose; fprintf('\nusing BIC best model order = %d\n',morder); end
else
    if verbose; fprintf('\nusing specified model order = %d\n',morder); end
end

% VAR model estimation (<mvgc_schema.html#3 |A2|>)

% Estimate VAR model of selected order from data.

if verbose; ptic('\n*** tsdata_to_var... '); end
[A,SIG] = tsdata_to_var(X,morder,regmode);
assert(~isbad(A),'VAR estimation failed - bailing out');
if verbose; ptoc; end

% Report information on the estimated VAR, and check for errors.
%
% _IMPORTANT:_ We check the VAR model for stability and symmetric
% positive-definite residuals covariance matrix. _THIS CHECK SHOULD ALWAYS BE
% PERFORMED!_ - subsequent routines may fail if there are errors here. If there
% are problems with the data (e.g. non-stationarity, colinearity, etc.) there's
% also a good chance they'll show up at this point - and the diagnostics may
% supply useful information as to what went wrong.

info = var_info(A,SIG,verbose);
assert(~info.error,'VAR error(s) found - bailing out');

%% Granger causality calculation: time domain  (<mvgc_schema.html#3 |A13|>)

% Calculate time-domain pairwise-conditional causalities from VAR model parameters
% by state-space method [4]. The VAR model is transformed into an equivalent state-
% space model for computation. Also return p-values for specified test (F-test or
% likelihood-ratio test; this is optional - if p-values are not required, then it
% is not necessary to supply time series |X|, regression mode |regmode|, or test
% specification |tstat|).

if verbose; ptic('*** var_to_pwcgc... '); end
[F,pval] = var_to_pwcgc(A,SIG,X,regmode,tstat);
if verbose; ptoc; end

% Check for failed GC calculation

assert(~isbad(F,false),'GC calculation failed - bailing out');

% Significance-test p-values, correcting for multiple hypotheses.

sig = significance(pval,alpha,mhtc);

% Plot time-domain causal graph, p-values and significance.

if plotting==true
figure(2); clf;
sgtitlex('Pairwise-conditional Granger causality - time domain');
subplot(1,3,1);
plot_pw(F);
set(gca, 'CLim', [0 nanmax(F(:))], 'XTickLabel', source_labels, 'YTickLabel', source_labels)

title('Pairwise-conditional GC');
subplot(1,3,2);
plot_pw(pval);

set(gca, 'CLim', [-.1 1], 'XTickLabel', source_labels, 'YTickLabel', source_labels)
title(['p-values (' tstat '-test)']);
subplot(1,3,3);
plot_pw(sig);
set(gca, 'CLim', [-.1 .5], 'XTickLabel', source_labels, 'YTickLabel', source_labels)
title(['Significant at \alpha = ' num2str(alpha)]);
colormap(viridis)
end

if false
    %% Granger causality calculation: frequency domain  (<mvgc_schema.html#3 |A14|>)
    
    % If not specified, we set the frequency resolution to something sensible. Warn if
    % resolution is very large, as this may lead to excessively long computation times,
    % and/or out-of-memory issues.
    
    if isempty(fres)
        fres = 2^nextpow2(info.acdec); % based on autocorrelation decay; alternatively, you could try fres = 2^nextpow2(nobs);
        if verbose; fprintf('\nfrequency resolution auto-calculated as %d (increments ~ %.2gHz)\n',fres,fs/2/fres); end
    end
    if fres > 20000 % adjust to taste
        if verbose; fprintf(2,'\nWARNING: large frequency resolution = %d - may cause computation time/memory usage problems\nAre you sure you wish to continue [y/n]? ',fres); end
        istr = input(' ','s'); if isempty(istr) || ~strcmpi(istr,'y'); fprintf(2,'Aborting...\n'); return; end
    end
    
    % Calculate spectral pairwise-conditional causalities at given frequency
    % resolution by state-space method.
    
    if verbose; ptic('\n*** var_to_spwcgc... '); end
    f = var_to_spwcgc(A,SIG,fres);
    assert(~isbad(f,false),'spectral GC calculation failed - bailing out');
    if verbose; ptoc; end
    
    % Plot spectral causal graph.
    
    if plotting==true
        figure(3); clf;
        sgtitlex('Pairwise-conditional Granger causality - frequency domain');
        plot_spw(f,fs);
    end
    %% Granger causality calculation: frequency domain -> time-domain  (<mvgc_schema.html#3 |A15|>)
    
    % Check that spectral causalities average (integrate) to time-domain
    % causalities, as they should according to theory.
    
    if verbose; fprintf('\nfrequency-domain GC integration check... '); end
    Fint = smvgc_to_mvgc(f); % integrate spectral MVGCs
    amax = maxabs(F+Fint)/2;
    if amax < 1e-5; amax = 1; end % in case all GCs very small
    mre = maxabs(F-Fint)/amax;
    if mre < 1e-5
        fprintf('OK (maximum relative error ~ %.0e)\n',mre);
    else
        fprintf(2,'WARNING: high maximum relative error ~ %.0e\n',mre);
    end
end


%%

outstr = [];
outstr.F    = F;
outstr.pval = pval;
outstr.sig  = sig;
outstr.nrand  = nrand;
outstr.prob  = NaN(nvars);

randstr      = [];

if nrand>0    
    F_rand = NaN(nrand, nvars, nvars);
    pval_rand = NaN(nrand, nvars, nvars);
    sig_rand = NaN(nrand, nvars, nvars);
    cshifts = randi(round([size(Xog,2)*.1, size(Xog,2)*.9]), [nvars, nrand]);
    % Xog = X;
    for randLoop = 1:nrand
        if nvars==2 && ndims(Xog)==3 %alreadys in source x smaple x trial form
            X_rand = Xog;
            ntrials = size(Xog,3);
            cshifts = randi(round([size(Xog,2)*.1, size(Xog,2)*.9]), [2, ntrials, nrand]);
            v = mod(randLoop,2)+1; % alternated shuffled matrix
            xtemp = squeeze(Xog(v,:,:));
            for i = 1:ntrials
                X_rand(v,:,i) = circshift(xtemp(:,i), cshifts(v, i, randLoop));
            end
        else
            
        X2 = Xog;
        if nvars==2 % only need to shift one
            X2(1,:) = circshift(Xog(1,:), cshifts(1, randLoop));
        else
            for i = 1:nvars
                X2(i,:) = circshift(Xog(i,:), cshifts(i, randLoop));
            end
        end
        if ~isempty(nsamples) % break up X into trials of length nsamples
            X_rand = NaN(nvars, nsamples, n_trials);
            for i = 1:n_trials
                % this method can sometimes shift things by the same
                % ammount, thus no overall shift
                X_rand(:,:,i) = X2(:, nsamples*(i-1) + 1: nsamples*i);
            end
        else
            X_rand = X2;
        end
        end
        % granger on rand data
%         [AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X, momax, icregmode);
        [A_rand,SIG_rand] = tsdata_to_var(X_rand, morder, regmode);
        [F_rand(randLoop, :, :), pval_rand(randLoop, :, :)] = var_to_pwcgc(A_rand, SIG_rand, X_rand, regmode, tstat);
        sig_rand(randLoop, :, :) = significance(pval_rand(randLoop, :, :), alpha, mhtc);
    end
    randstr.F    = F_rand;
    randstr.pval = pval_rand;
    randstr.sig  = sig_rand;
    %%
%     outstr.prob  = NaN(nvars);
    if plotting; figure(109); clf; end
    for i = 1:nvars
        for j = 1:nvars
            if i~=j
                outstr.prob(i,j) = sum(F(i,j)<F_rand(:,i,j))/nrand;
                if plotting
                    figure(109);
                    d = F_rand(:,i,j);
                    subplot(nvars, nvars, sub2ind([nvars nvars], j, i)); cla
%                     histogram(d, linspace(0, quantile(d, .95), 50), 'Normalization', 'probability')
                    [c, bins] = histcounts(d, 100, 'Normalization', 'probability');
                    dord = sort(d,'descend');
%                     dord = dord./max(dord);
                    hold on;
                    bins = bins(1:end-1) + abs(mean(diff(bins)));
                    plot(bins, c, 'k-')
                    plot([F(i,j) F(i,j)], [0 1.2*max(c)], 'r-')
%                     plot([outstr.prob(i,j) outstr.prob(i,j)], [0, F(i,j)], 'r-', 'LineWidth', 2)
%                     plot(linspace(0,1,nrand), dord, 'k-')
%                     xlim([-.2 1.2])
%                     ylim([-.2*max(d) 1.5*max(d)])
                    title(sprintf('rand prob = %2.5f', outstr.prob(i,j)))
                end
            end
        end
    end
end
% figure;

outstr.setup.regmode = regmode;
outstr.setup.icregmode = icregmode;
outstr.setup.morder  = morder;
outstr.setup.momax  = momax;
outstr.setup.acmaxlags  = acmaxlags;
outstr.setup.tstat = tstat;
outstr.setup.alpha = alpha;
outstr.setup.mhtc  = mhtc;
outstr.setup.seed = seed;

% outstr.setup.fres  = fres;
