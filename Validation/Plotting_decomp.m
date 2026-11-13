%% slice images
decay_vec=decay_min:(decay_max-decay_min)/m:decay_max;

figure()
t = tiledlayout(3,4);
title(t,"Net H2O2 directly via ethanol per decay value (k_d), one-hit")

for i = 1:length(decay_vec)
    ax = nexttile;
    
    idx = find(abs(decay_list - decay_vec(i)) < 1e-12);
    
    slice = reshape(total_H2O2_case1(idx), 11, 11);
    corresponding_y1 = reshape(y1_list(idx), 11, 11);
    corresponding_y2 = reshape(y2_list(idx), 11, 11);
    
    imagesc( corresponding_y2(1,:),corresponding_y1(:,1), (slice))
    set(ax, 'YDir', 'normal'); % Set to normal (bottom-to-top)
    title(ax, sprintf("k_d = %g", decay_vec(i)))
    xlabel(ax, "y_2")
    ylabel(ax, "y_1")
    colorbar()
end

figure()
t = tiledlayout(3,4);
title(t,"Net H2O2 directly via ethanol per decay value (k_d), two-hit")

for i = 1:length(decay_vec)
    ax = nexttile;
    
    idx = find(abs(decay_list - decay_vec(i)) < 1e-12);
    
    slice = reshape(total_H2O2_case2(idx), 11, 11);
    corresponding_y1 = reshape(y1_list(idx), 11, 11);
    corresponding_y2 = reshape(y2_list(idx), 11, 11);
    
    imagesc( corresponding_y2(1,:),corresponding_y1(:,1), (slice))
    set(ax, 'YDir', 'normal'); % Set to normal (bottom-to-top)
    title(ax, sprintf("k_d = %g", decay_vec(i)))
    xlabel(ax, "y_2")
    ylabel(ax, "y_1")
    colorbar()
end


figure()
t = tiledlayout(3,4);
title(t,"Net oxidative stress per decay value (k_d), one-hit")

for i = 1:length(decay_vec)
    ax = nexttile;
    
    idx = find(abs(decay_list - decay_vec(i)) < 1e-12);
    
    slice = reshape(total_oxidative_stress_case1(idx), 11, 11);
    corresponding_y1 = reshape(y1_list(idx), 11, 11);
    corresponding_y2 = reshape(y2_list(idx), 11, 11);
    
    imagesc( corresponding_y2(1,:),corresponding_y1(:,1), (slice))
    set(ax, 'YDir', 'normal'); % Set to normal (bottom-to-top)
    title(ax, sprintf("k_d = %g", decay_vec(i)))
    xlabel(ax, "y_2")
    ylabel(ax, "y_1")
    colorbar()
end


figure()
t = tiledlayout(3,4);
title(t,"Net oxidative stress per decay value (k_d), two-hit")

for i = 1:length(decay_vec)
    ax = nexttile;
    
    idx = find(abs(decay_list - decay_vec(i)) < 1e-12);
    
    slice = reshape(total_oxidative_stress_case2(idx), 11, 11);
    corresponding_y1 = reshape(y1_list(idx), 11, 11);
    corresponding_y2 = reshape(y2_list(idx), 11, 11);
    
    imagesc( corresponding_y2(1,:),corresponding_y1(:,1), (slice))
    set(ax, 'YDir', 'normal'); % Set to normal (bottom-to-top)
    title(ax, sprintf("k_d = %g", decay_vec(i)))
    xlabel(ax, "y_2")
    ylabel(ax, "y_1")
    colorbar()
end

