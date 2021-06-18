function power = hepsret_simulate_detect_prob(effect_diff, std_eff, ...
    ntrials_cond, nsubj, nperms, alpha)
% simulation for assess how likely it is to detect an effect

hperm = NaN(1,nperms);

for np = 1:nperms
    msc1 = NaN(1, nsubj); msc2=msc1;
    for ns = 1:nsubj
        simc1 = -effect_diff/2 + std_eff(ns)*randn(ntrials_cond(ns,1),1);
        simc2 = +effect_diff/2 + std_eff(ns)*randn(ntrials_cond(ns,2),1);
        msc1(ns) = mean(simc1,1);
        msc2(ns) = mean(simc2,1);
    end
    hperm(np) = ttest(msc2, msc1, 'alpha', alpha);
end

power = sum(hperm)/numel(hperm);
fprintf('probability of detecting effect (power) = %g\n', power)

end
