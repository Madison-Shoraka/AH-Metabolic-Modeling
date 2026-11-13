% 1. Plot the main matrix A
[initial_decomp,idx]=sort(met_decomp_0treatment(param_unstab));
treatment1=max(metabolic_decomp_1treatment')';
treatment1=treatment1(param_unstab);
treatment1=treatment1(idx);
treatment2=max(metabolic_decomp_2treatment')';
treatment2=treatment2(param_unstab);
treatment2=treatment2(idx);
treatment3=max(metabolic_decomp_3treatment')';
treatment3=treatment3(param_unstab);
treatment3=treatment3(idx);
A=[initial_decomp,treatment1,treatment2,treatment3];
imagesc(A);
colorbar;
hold on;
% 2. Find the row and column coordinates where B is 1
B=A>(initial_decomp);

[rows, cols] = find(B == 1);
% 3. Loop through and draw a red box around each tile
for i = 1:length(rows)
    r = rows(i);
    c = cols(i);
    % In imagesc, X is columns and Y is rows.
    % A single tile centered at (c, r) spans from -0.5 to +0.5.
    rectangle('Position', [c - 0.5, r - 0.5, 1, 1], 'EdgeColor', 'r', 'LineWidth', 2);
end
hold off;
disp("Percent of each treatment that is successful (better than none)")
sum(1-B)./length_no_stable


%%

y1=y1_list(index_failure_case1);
y1=y1(param_unstab);
y2=y2_list(index_failure_case1);
y2=y2(param_unstab);
kd=decay_list(index_failure_case1);
kd=kd(param_unstab);

X=[y1,y2,kd,met_decomp_0treatment(param_unstab),met_decomp_1treatment(param_unstab)<met_decomp_0treatment(param_unstab),met_decomp_2treatment(param_unstab)<met_decomp_0treatment(param_unstab),met_decomp_3treatment(param_unstab)<met_decomp_0treatment(param_unstab)];
[coeff,score,latent,tsquared,explained] = pca(X);
scatter(score(:,1),score(:,2),[],met_decomp_0treatment(param_unstab),'filled')
axis equal
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
zlabel('3rd Principal Component')
