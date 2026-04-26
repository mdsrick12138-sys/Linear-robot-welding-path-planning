clc; clear; close all;

%% 参数设置
input_file = 'D:\BY\ZXGH\point_cloud_00004sss.ply';
output_dir = 'D:\BY\ZXGH\clusters_trajectories\';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

z_offsets = [1.5, 1.5, 1.5];
num_slices = 100;
minDistance = 2.0;

region_params = [
    1/3, 0.3, 2;
    2/3, 0.34, 2;
    
];

numPoints = 100; % 重采样点数

%% 读点云
ptCloud_full = pcread(input_file);
points_full = ptCloud_full.Location;
z_vals_full = points_full(:,3);
max_z = max(z_vals_full);
min_z = min(z_vals_full);
z_range = max_z - min_z;

%% 顶部轨迹提取
z_threshold = max_z - 0.8;
top_points = points_full(z_vals_full > z_threshold, :);
x_edges = linspace(min(top_points(:,1)), max(top_points(:,1)), num_slices+1);
trajectory_simple = [];
for i = 1:num_slices
    in_slice = top_points(:,1) >= x_edges(i) & top_points(:,1) < x_edges(i+1);
    slice_points = top_points(in_slice, :);
    if ~isempty(slice_points)
        center_pt = mean(slice_points, 1);
        trajectory_simple = [trajectory_simple; center_pt];
    end
end

original_top_trajectory = trajectory_simple;
trajectory_simple(:,3) = trajectory_simple(:,3) - 3.0;

allTrajectories = {};
allTrajectoryPoints = {};
trajectoryCount = 1;
allTrajectories{trajectoryCount} = trajectory_simple;
allTrajectoryPoints{trajectoryCount} = top_points;
writematrix(trajectory_simple, fullfile(output_dir, 'trajectory_top_simple.csv'));

%% 区域聚类轨迹提取
for regionIdx = 1:size(region_params,1)
    region_pos = region_params(regionIdx, 1);
    thickness = region_params(regionIdx, 2);
    expected_trajectories = region_params(regionIdx, 3);

    z_target = max_z - z_range * region_pos;
    region_mask = (z_vals_full > (z_target - thickness)) & (z_vals_full < (z_target + thickness));
    region_points = points_full(region_mask, :);
    ptCloud_region = pointCloud(region_points);

    [labels, numClusters] = pcsegdist(ptCloud_region, minDistance);
    found = 0;

    for clusterIdx = 1:numClusters
        idx = (labels == clusterIdx);
        cluster_points = region_points(idx, :);
        if isempty(cluster_points), continue; end

        z_vals = cluster_points(:,3);
        max_z_cluster = max(z_vals);
        z_thresh = max_z_cluster - z_offsets(min(regionIdx, numel(z_offsets)));
        high_points = cluster_points(z_vals > z_thresh, :);
        if size(high_points,1) < 10, continue; end

        x_edges = linspace(min(high_points(:,1)), max(high_points(:,1)), num_slices+1);
        trajectory = [];
        for i = 1:num_slices
            in_slice = high_points(:,1) >= x_edges(i) & high_points(:,1) < x_edges(i+1);
            slice_points = high_points(in_slice, :);
            if ~isempty(slice_points)
                center_pt = mean(slice_points,1);
                trajectory = [trajectory; center_pt];
            end
        end
        if size(trajectory,1) < 2, continue; end

        trajectory(:,3) = trajectory(:,3) - 3.0;

        trajectoryCount = trajectoryCount + 1;
        allTrajectories{trajectoryCount} = trajectory;
        allTrajectoryPoints{trajectoryCount} = high_points;

        writematrix(trajectory, fullfile(output_dir, ...
            sprintf('trajectory_r%d_c%d.csv', regionIdx, clusterIdx)));

        found = found + 1;
        if found >= expected_trajectories
            break;
        end
    end
end

%% 插入顶部轨迹副本至 2/3 高度
middle_traj = original_top_trajectory;
target_z_middle = max_z - z_range * (2/3);
mean_z_orig = mean(middle_traj(:,3));
z_shift = target_z_middle - mean_z_orig;
middle_traj(:,3) = middle_traj(:,3) + z_shift;
middle_traj(:,3) = middle_traj(:,3) - 3.0;

trajectoryCount = trajectoryCount + 1;
allTrajectories{trajectoryCount} = middle_traj;
allTrajectoryPoints{trajectoryCount} = top_points;
writematrix(middle_traj, fullfile(output_dir, 'trajectory_middle_inserted.csv'));

%% 重采样轨迹
for i = 1:length(allTrajectories)
    allTrajectories{i} = resampleTrajectory(allTrajectories{i}, numPoints);
end

%% 图 1：原始点云
figure;
pcshow(ptCloud_full);
title('图 1：原始点云');
xlabel('X'); ylabel('Y'); zlabel('Z'); view(3); axis equal; grid on;

%% 图 2：仅拟合轨迹
figure; hold on;
colors_map = lines(length(allTrajectories));
for i = 1:length(allTrajectories)
    plot3(allTrajectories{i}(:,1), allTrajectories{i}(:,2), allTrajectories{i}(:,3), ...
        '-', 'Color', colors_map(i,:), 'LineWidth', 1.5);
end
xlabel('X'); ylabel('Y'); zlabel('Z');
title('图 2：仅拟合轨迹');
grid on; axis equal; view(3);

%% 图 2（旋转 180°）：仅拟合轨迹（绕X轴）
figure; hold on;

% 创建绕X轴旋转180度矩阵
theta = deg2rad(180);  % 180°
R_x = [1, 0,         0;
       0, cos(theta), -sin(theta);
       0, sin(theta),  cos(theta)];

% 计算轨迹中心（可选，让旋转围绕整体中心）
all_points_cat = cell2mat(allTrajectories'); 
centroid_traj = mean(all_points_cat, 1);

% 绘制旋转后的轨迹
for i = 1:length(allTrajectories)
    traj_rot = (allTrajectories{i} - centroid_traj) * R_x + centroid_traj;
    plot3(traj_rot(:,1), traj_rot(:,2), traj_rot(:,3), '-', 'LineWidth', 1.8);
end

xlabel('X'); ylabel('Y'); zlabel('Z');
title('图 2（旋转 180°）：仅拟合轨迹');
grid on; axis equal; view(3);


%% 图 3：拟合轨迹 + 拟合点云
figure; hold on;
for i = 1:length(allTrajectoryPoints)
    scatter3(allTrajectoryPoints{i}(:,1), allTrajectoryPoints{i}(:,2), allTrajectoryPoints{i}(:,3), ...
        3, colors_map(i,:), 'filled');
end
for i = 1:length(allTrajectories)
    plot3(allTrajectories{i}(:,1), allTrajectories{i}(:,2), allTrajectories{i}(:,3), ...
        '-', 'Color', colors_map(i,:), 'LineWidth', 2);
end
xlabel('X'); ylabel('Y'); zlabel('Z');
title('图 3：拟合轨迹 + 拟合点云');
grid on; axis equal; view(3);


%% 图 4：拟合轨迹 + 原始点云
figure; hold on;
pcshow(ptCloud_full, 'MarkerSize', 10);
for i = 1:length(allTrajectories)
    plot3(allTrajectories{i}(:,1), allTrajectories{i}(:,2), allTrajectories{i}(:,3), ...
        '-', 'Color', colors_map(i,:), 'LineWidth', 2);
end
xlabel('X'); ylabel('Y'); zlabel('Z');
title('图 4：拟合轨迹 + 原始点云');
grid on; axis equal; view(3);

%% 构建完整闭环焊接轨迹
start_top = [-145, -15, 430];
fullTrajectory = allTrajectories{1};
layer_indices = ones(size(fullTrajectory,1),1);
for i = 2:length(allTrajectories)
    traj = allTrajectories{i};
    if norm(fullTrajectory(end,:) - traj(end,:)) < norm(fullTrajectory(end,:) - traj(1,:))
        traj = flipud(traj);
    end
    fullTrajectory = [fullTrajectory; traj];
    layer_indices = [layer_indices; i*ones(size(traj,1),1)];
end

num_interp = 50;
segment_start = [linspace(start_top(1), fullTrajectory(1,1), num_interp)', ...
                 linspace(start_top(2), fullTrajectory(1,2), num_interp)', ...
                 linspace(start_top(3), fullTrajectory(1,3), num_interp)'];
segment_end = [linspace(fullTrajectory(end,1), start_top(1), num_interp)', ...
               linspace(fullTrajectory(end,2), start_top(2), num_interp)', ...
               linspace(fullTrajectory(end,3), start_top(3), num_interp)'];

fullTrajectory_closed = [segment_start; fullTrajectory; segment_end];

numPointsClosed = 2000;
fullTrajectory_closed = resampleTrajectory(fullTrajectory_closed, numPointsClosed);
layer_indices_closed = round(linspace(0, length(allTrajectories), numPointsClosed))';

writematrix(fullTrajectory_closed, fullfile(output_dir, 'full_welding_trajectory_closed.csv'));
%% 轨迹反转功能：反转闭环轨迹点顺序并保存
fullTrajectory_closed_reversed = flipud(fullTrajectory_closed);
layer_indices_closed_reversed = flipud(layer_indices_closed);  % 可选保留原layer顺序或重新编号

% 保存为新的CSV文件
writematrix(fullTrajectory_closed_reversed, fullfile(output_dir, 'full_welding_trajectory_closed_reversed.csv'));

% 保存为新的TXT文件（带层号）
txt_filename_reversed = fullfile(output_dir, 'full_welding_trajectory_closed_reversed_with_layers.txt');
fileID = fopen(txt_filename_reversed, 'w');
fprintf(fileID, 'X\tY\tZ\tLayer\n');
for i = 1:size(fullTrajectory_closed_reversed,1)
    fprintf(fileID, '%.6f\t%.6f\t%.6f\t%d\n', ...
        fullTrajectory_closed_reversed(i,1), ...
        fullTrajectory_closed_reversed(i,2), ...
        fullTrajectory_closed_reversed(i,3), ...
        layer_indices_closed_reversed(i));
end
fclose(fileID);
fprintf('反转后的闭环焊接轨迹保存至：\n%s\n', txt_filename_reversed);

txt_filename = fullfile(output_dir, 'full_welding_trajectory_closed_with_layers.txt');
fileID = fopen(txt_filename, 'w');
fprintf(fileID, 'X\tY\tZ\tLayer\n');
for i = 1:size(fullTrajectory_closed,1)
    fprintf(fileID, '%.6f\t%.6f\t%.6f\t%d\n', ...
        fullTrajectory_closed(i,1), fullTrajectory_closed(i,2), fullTrajectory_closed(i,3), ...
        layer_indices_closed(i));
end
fclose(fileID);
fprintf('闭环焊接轨迹保存至：\n%s\n', txt_filename);

%% 构建包含轨迹的点云
trajectory_cloud = fullTrajectory_closed;
colors = repmat(uint8([255 0 0]), size(trajectory_cloud, 1), 1);
color_full = repmat(uint8([100 100 100]), size(points_full, 1), 1);

combined_points = [points_full; trajectory_cloud];
combined_colors = [color_full; colors];

ptCloud_combined = pointCloud(combined_points, 'Color', combined_colors);
pcwrite(ptCloud_combined, fullfile(output_dir, 'point_cloud_with_trajectory.ply'));

txt_filename_traj = fullfile(output_dir, 'combined_point_cloud_with_trajectory.txt');
fileID = fopen(txt_filename_traj, 'w');
fprintf(fileID, 'X\tY\tZ\n');
for i = 1:size(combined_points,1)
    fprintf(fileID, '%.6f\t%.6f\t%.6f\n', ...
        combined_points(i,1), combined_points(i,2), combined_points(i,3));
end
fclose(fileID);
fprintf('TXT格式点云（包含轨迹）已保存至：\n%s\n', txt_filename_traj);

%% 图 5：完整闭环焊接轨迹
figure; hold on;

plot3(fullTrajectory_closed(:,1), fullTrajectory_closed(:,2), fullTrajectory_closed(:,3), ...
    '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);

skip = 6;
idx = 1:skip:size(fullTrajectory_closed,1);
plot3(fullTrajectory_closed(idx,1), fullTrajectory_closed(idx,2), fullTrajectory_closed(idx,3), ...
    '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 15, 'LineWidth', 1.5);

for i = 1:min(length(allTrajectories), 6)
    traj = allTrajectories{i};
    plot3(traj(:,1), traj(:,2), traj(:,3), 'o', 'MarkerSize', 4, ...
        'Color', [rand rand rand], 'MarkerFaceColor', [rand rand rand]);
end

plot3(start_top(1), start_top(2), start_top(3), 'ms', 'MarkerSize', 10, 'MarkerFaceColor', 'm');
text(start_top(1), start_top(2), start_top(3)+1, 'StartTop', 'FontSize', 12, 'Color', 'm');

plot3(fullTrajectory(1,1), fullTrajectory(1,2), fullTrajectory(1,3), 'go', ...
    'MarkerSize', 10, 'MarkerFaceColor', 'g');
text(fullTrajectory(1,1), fullTrajectory(1,2), fullTrajectory(1,3), 'Start', 'FontSize', 12, 'Color', 'g');

plot3(fullTrajectory(end,1), fullTrajectory(end,2), fullTrajectory(end,3), 'bo', ...
    'MarkerSize', 10, 'MarkerFaceColor', 'b');
text(fullTrajectory(end,1), fullTrajectory(end,2), fullTrajectory(end,3), 'End', 'FontSize', 12, 'Color', 'b');

quiver3(start_top(1), start_top(2), start_top(3), ...
       fullTrajectory(1,1)-start_top(1), ...
       fullTrajectory(1,2)-start_top(2), ...
       fullTrajectory(1,3)-start_top(3), ...
       0, 'Color', 'k', 'LineWidth', 2, 'MaxHeadSize', 2);

xlabel('X'); ylabel('Y'); zlabel('Z');
grid on; axis equal;
title('图 5：完整闭环焊接轨迹（起点→焊接→回到起点）');
view(3);

%% 轨迹反转前后对比可视化
figure; hold on; grid on; axis equal;
title('轨迹反转前后对比 (相同路径，方向相反)');
xlabel('X'); ylabel('Y'); zlabel('Z');
view(3);

% 原始轨迹（反转前）
plot3(fullTrajectory_closed(:,1), fullTrajectory_closed(:,2), fullTrajectory_closed(:,3), ...
    '--b', 'LineWidth', 2, 'DisplayName', '反转前轨迹 (起点→终点)');

% 反转前轨迹方向指示
quiver3(fullTrajectory_closed(1,1), fullTrajectory_closed(1,2), fullTrajectory_closed(1,3), ...
    fullTrajectory_closed(10,1)-fullTrajectory_closed(1,1), ...
    fullTrajectory_closed(10,2)-fullTrajectory_closed(1,2), ...
    fullTrajectory_closed(10,3)-fullTrajectory_closed(1,3), ...
    0, 'b', 'LineWidth', 1.5, 'MaxHeadSize', 2, 'DisplayName', '反转前方向');

% 反转后轨迹
plot3(fullTrajectory_closed_reversed(:,1), fullTrajectory_closed_reversed(:,2), fullTrajectory_closed_reversed(:,3), ...
    '-r', 'LineWidth', 1.5, 'DisplayName', '反转后轨迹 (起点→终点)');

% 反转后方向指示
quiver3(fullTrajectory_closed_reversed(1,1), fullTrajectory_closed_reversed(1,2), fullTrajectory_closed_reversed(1,3), ...
    fullTrajectory_closed_reversed(10,1)-fullTrajectory_closed_reversed(1,1), ...
    fullTrajectory_closed_reversed(10,2)-fullTrajectory_closed_reversed(1,2), ...
    fullTrajectory_closed_reversed(10,3)-fullTrajectory_closed_reversed(1,3), ...
    0, 'r', 'LineWidth', 1.5, 'MaxHeadSize', 2, 'DisplayName', '反转后方向');

% 关键点标注
plot3(start_top(1), start_top(2), start_top(3), 'pm', 'MarkerSize', 12, 'MarkerFaceColor', 'm', 'DisplayName', '全局起点');
text(start_top(1), start_top(2), start_top(3)+5, '全局起点', 'FontSize', 10, 'Color', 'm');

% 反转前起终点
plot3(fullTrajectory_closed(1,1), fullTrajectory_closed(1,2), fullTrajectory_closed(1,3), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
text(fullTrajectory_closed(1,1), fullTrajectory_closed(1,2), fullTrajectory_closed(1,3)-5, '反转前起点', 'FontSize', 10, 'Color', 'g');

plot3(fullTrajectory_closed(end,1), fullTrajectory_closed(end,2), fullTrajectory_closed(end,3), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
text(fullTrajectory_closed(end,1), fullTrajectory_closed(end,2), fullTrajectory_closed(end,3)-5, '反转前终点', 'FontSize', 10, 'Color', 'b');

% 反转后起终点
plot3(fullTrajectory_closed_reversed(1,1), fullTrajectory_closed_reversed(1,2), fullTrajectory_closed_reversed(1,3), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(fullTrajectory_closed_reversed(1,1), fullTrajectory_closed_reversed(1,2), fullTrajectory_closed_reversed(1,3)+5, '反转后起点', 'FontSize', 10, 'Color', 'r');

plot3(fullTrajectory_closed_reversed(end,1), fullTrajectory_closed_reversed(end,2), fullTrajectory_closed_reversed(end,3), 'co', 'MarkerSize', 8, 'MarkerFaceColor', 'c');
text(fullTrajectory_closed_reversed(end,1), fullTrajectory_closed_reversed(end,2), fullTrajectory_closed_reversed(end,3)+5, '反转后终点', 'FontSize', 10, 'Color', 'c');

legend('Location', 'best');
view(-30, 30);

%% 旋转功能：绕X轴旋转180度 (保持形状不变)
rotation_angle = 180; % 旋转角度（度）
rotation_axis = 'x';  % 旋转轴 ('x', 'y', 'z')

% 计算点云中心
centroid = mean(points_full, 1);

% 创建旋转矩阵
theta = deg2rad(rotation_angle);
switch lower(rotation_axis)
    case 'x'
        R = [1, 0, 0;
             0, cos(theta), -sin(theta);
             0, sin(theta), cos(theta)];
    case 'y'
        R = [cos(theta), 0, sin(theta);
             0, 1, 0;
             -sin(theta), 0, cos(theta)];
    case 'z'
        R = [cos(theta), -sin(theta), 0;
             sin(theta), cos(theta), 0;
             0, 0, 1];
end

% 对原始点云应用旋转
points_rotated = (points_full - centroid) * R + centroid;
ptCloud_rotated = pointCloud(points_rotated);

% 对轨迹点应用同样的旋转
fullTrajectory_closed_rotated = (fullTrajectory_closed - centroid) * R + centroid;

% 保存旋转后的轨迹
writematrix(fullTrajectory_closed_rotated, fullfile(output_dir, 'full_welding_trajectory_rotated.csv'));

% 保存带层号的旋转轨迹
txt_filename_rotated = fullfile(output_dir, 'full_welding_trajectory_rotated_with_layers.txt');
fileID = fopen(txt_filename_rotated, 'w');
fprintf(fileID, 'X\tY\tZ\tLayer\n');
for i = 1:size(fullTrajectory_closed_rotated,1)
    fprintf(fileID, '%.6f\t%.6f\t%.6f\t%d\n', ...
        fullTrajectory_closed_rotated(i,1), ...
        fullTrajectory_closed_rotated(i,2), ...
        fullTrajectory_closed_rotated(i,3), ...
        layer_indices_closed(i));
end
fclose(fileID);

%% 图 6：旋转后的点云
figure;
pcshow(ptCloud_rotated);
title(sprintf('图 6：旋转后的点云 (绕%s轴旋转%d度)', upper(rotation_axis), rotation_angle));
xlabel('X'); ylabel('Y'); zlabel('Z'); view(3); axis equal; grid on;

%% 图 7：旋转后的轨迹
%% ...（前面代码保持不变）...
%% 图 7：旋转后的轨迹 + 拟合点云
%% 图 7：旋转后的轨迹 + 拟合点云
figure; hold on;

% 设置颜色映射
colors_map = lines(length(allTrajectories));

% 保存文件路径
rotated_traj_dir = 'C:\Users\bighero\Desktop\Si3N4_deep\ply_file\poukou\clusters_trajectories\';

% 如果目录不存在，创建它
if ~exist(rotated_traj_dir, 'dir')
    mkdir(rotated_traj_dir);
end

% 初始化合并数组
all_rotated_traj_concat = [];
all_rotated_points_concat = [];

for i = 1:length(allTrajectories)
    % 应用旋转到轨迹点
    traj_rotated = (allTrajectories{i} - centroid) * R + centroid;
    
    % 应用旋转到对应点云
    points_rotated = (allTrajectoryPoints{i} - centroid) * R + centroid;

    % 可视化点云
    scatter3(points_rotated(:,1), points_rotated(:,2), points_rotated(:,3), ...
        3, colors_map(i,:), 'filled', 'MarkerEdgeAlpha', 0.3, 'MarkerFaceAlpha', 0.3);
    
    % 可视化轨迹
    plot3(traj_rotated(:,1), traj_rotated(:,2), traj_rotated(:,3), ...
        '-', 'Color', colors_map(i,:), 'LineWidth', 2.5);

    % 累积轨迹和点云
    all_rotated_traj_concat = [all_rotated_traj_concat; traj_rotated];
    all_rotated_points_concat = [all_rotated_points_concat; points_rotated];
end

% --- 保存所有轨迹到一个 CSV 文件 ---
csv_all_traj = fullfile(rotated_traj_dir, 'all_rotated_trajectories.csv');
writematrix(all_rotated_traj_concat, csv_all_traj);
fprintf('✅ 已保存所有旋转轨迹到: %s\n', csv_all_traj);

% --- 保存所有点云为一个 PLY 文件 ---
ptCloud_all = pointCloud(all_rotated_points_concat);
ply_all_points = fullfile(rotated_traj_dir, 'all_rotated_pointcloud.ply');
pcwrite(ptCloud_all, ply_all_points);
fprintf('✅ 已保存所有拟合点云为: %s\n', ply_all_points);
% --- 将拟合点云在Z轴方向下移3mm ---
offset_z = +3;  % 向下为负方向
translated_points = all_rotated_points_concat;
translated_points(:,3) = translated_points(:,3) + offset_z;

% 创建平移后的点云对象
ptCloud_translated = pointCloud(translated_points);

% 保存仅包含平移后的拟合点云的PLY文件
ply_translated_points = fullfile(rotated_traj_dir, 'translated_only_pointcloud.ply');
pcwrite(ptCloud_translated, ply_translated_points);
fprintf('✅ 已保存Z轴下移3mm后的点云为: %s\n', ply_translated_points);

% --- 图形增强 ---
xlabel('X'); ylabel('Y'); zlabel('Z');
title(sprintf('图 7：旋转后的轨迹和拟合点云 (绕%s轴旋转%d度)', upper(rotation_axis), rotation_angle));
grid on; axis equal; view(3);

% 添加图例
legend_items = arrayfun(@(i) sprintf('轨迹 %d', i), 1:length(allTrajectories), 'UniformOutput', false);
h = legend(legend_items, 'Location', 'best');
set(h, 'FontSize', 10);

% 增强可视化效果
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
lighting gouraud;
material shiny;
camlight headlight;

%% ...（后续代码保持不变）...

%% ...（后续代码保持不变）...


%% 图 8：旋转前后的轨迹对比
figure; hold on; grid on; axis equal;
title(sprintf('旋转前后对比 (绕%s轴旋转%d度)', upper(rotation_axis), rotation_angle));
xlabel('X'); ylabel('Y'); zlabel('Z');

% 原始轨迹
plot3(fullTrajectory_closed(:,1), fullTrajectory_closed(:,2), fullTrajectory_closed(:,3), ...
    'b-', 'LineWidth', 2, 'DisplayName', '原始轨迹');

% 旋转后轨迹
plot3(fullTrajectory_closed_rotated(:,1), fullTrajectory_closed_rotated(:,2), fullTrajectory_closed_rotated(:,3), ...
    'r-', 'LineWidth', 2, 'DisplayName', '旋转后轨迹');

% 原始关键点
plot3(start_top(1), start_top(2), start_top(3), 'mo', 'MarkerSize', 10, ...
    'MarkerFaceColor', 'm', 'DisplayName', '原始起点');
plot3(fullTrajectory_closed(1,1), fullTrajectory_closed(1,2), fullTrajectory_closed(1,3), 'go', ...
    'MarkerSize', 8, 'MarkerFaceColor', 'g', 'DisplayName', '轨迹起点');
plot3(fullTrajectory_closed(end,1), fullTrajectory_closed(end,2), fullTrajectory_closed(end,3), 'bo', ...
    'MarkerSize', 8, 'MarkerFaceColor', 'b', 'DisplayName', '轨迹终点');

% 旋转后关键点
rotated_start = (start_top - centroid) * R + centroid;
rotated_traj_start = (fullTrajectory_closed(1,:) - centroid) * R + centroid;
rotated_traj_end = (fullTrajectory_closed(end,:) - centroid) * R + centroid;

plot3(rotated_start(1), rotated_start(2), rotated_start(3), 'm^', 'MarkerSize', 10, ...
    'MarkerFaceColor', 'm', 'DisplayName', '旋转后起点');
plot3(rotated_traj_start(1), rotated_traj_start(2), rotated_traj_start(3), 'g^', ...
    'MarkerSize', 8, 'MarkerFaceColor', 'g', 'DisplayName', '旋转后轨迹起点');
plot3(rotated_traj_end(1), rotated_traj_end(2), rotated_traj_end(3), 'b^', ...
    'MarkerSize', 8, 'MarkerFaceColor', 'b', 'DisplayName', '旋转后轨迹终点');

legend('Location', 'best');
view(3);

%% 镜像翻转功能（可选）
flip_axis = 'z'; % 可选 'x', 'y', 'z'
switch flip_axis
    case 'z'
        z_mirror = (max_z + min_z) / 2;
        fullTrajectory_closed_flipped = fullTrajectory_closed;
        fullTrajectory_closed_flipped(:,3) = 2 * z_mirror - fullTrajectory_closed(:,3);
    case 'y'
        y_vals = fullTrajectory_closed(:,2);
        y_mirror = (max(y_vals) + min(y_vals)) / 2;
        fullTrajectory_closed_flipped = fullTrajectory_closed;
        fullTrajectory_closed_flipped(:,2) = 2 * y_mirror - fullTrajectory_closed(:,2);
    case 'x'
        x_vals = fullTrajectory_closed(:,1);
        x_mirror = (max(x_vals) + min(x_vals)) / 2;
        fullTrajectory_closed_flipped = fullTrajectory_closed;
        fullTrajectory_closed_flipped(:,1) = 2 * x_mirror - fullTrajectory_closed(:,1);
    otherwise
        error('Unsupported flip axis');
end

% 重采样（可选，保持数量一致）
fullTrajectory_closed_flipped = resampleTrajectory(fullTrajectory_closed_flipped, numPointsClosed);

% 保存 flipped CSV
writematrix(fullTrajectory_closed_flipped, fullfile(output_dir, 'full_welding_trajectory_closed_flipped.csv'));

% 保存 flipped TXT 带层号
txt_filename_flipped = fullfile(output_dir, 'full_welding_trajectory_closed_flipped_with_layers.txt');
fileID = fopen(txt_filename_flipped, 'w');
fprintf(fileID, 'X\tY\tZ\tLayer\n');
for i = 1:size(fullTrajectory_closed_flipped,1)
    fprintf(fileID, '%.6f\t%.6f\t%.6f\t%d\n', ...
        fullTrajectory_closed_flipped(i,1), ...
        fullTrajectory_closed_flipped(i,2), ...
        fullTrajectory_closed_flipped(i,3), ...
        layer_indices_closed(i)); % 保留原layer编号
end
fclose(fileID);

fprintf('翻转后的闭环轨迹已保存至：\n%s\n', txt_filename_flipped);

%% 图 9：旋转后的完整闭环焊接轨迹（新图）
figure; hold on; grid on; axis equal;
title(sprintf('图 9：旋转后的完整闭环焊接轨迹 (绕%s轴旋转%d度)', upper(rotation_axis), rotation_angle));
xlabel('X'); ylabel('Y'); zlabel('Z');
view(3);

% 绘制旋转后的完整闭环轨迹
plot3(fullTrajectory_closed_rotated(:,1), fullTrajectory_closed_rotated(:,2), fullTrajectory_closed_rotated(:,3), ...
    'b-', 'LineWidth', 2, 'DisplayName', '旋转后轨迹');

% 绘制轨迹方向指示箭头
arrow_start_idx = 1;
arrow_end_idx = min(50, size(fullTrajectory_closed_rotated,1));
quiver3(fullTrajectory_closed_rotated(arrow_start_idx,1), fullTrajectory_closed_rotated(arrow_start_idx,2), fullTrajectory_closed_rotated(arrow_start_idx,3), ...
       fullTrajectory_closed_rotated(arrow_end_idx,1)-fullTrajectory_closed_rotated(arrow_start_idx,1), ...
       fullTrajectory_closed_rotated(arrow_end_idx,2)-fullTrajectory_closed_rotated(arrow_start_idx,2), ...
       fullTrajectory_closed_rotated(arrow_end_idx,3)-fullTrajectory_closed_rotated(arrow_start_idx,3), ...
       0, 'r', 'LineWidth', 1.5, 'MaxHeadSize', 1.5, 'DisplayName', '轨迹方向');

% 标记关键点
rotated_global_start = (start_top - centroid) * R + centroid;
rotated_traj_start = (fullTrajectory(1,:) - centroid) * R + centroid;
rotated_traj_end = (fullTrajectory(end,:) - centroid) * R + centroid;

% 全局起点
plot3(rotated_global_start(1), rotated_global_start(2), rotated_global_start(3), ...
    'pm', 'MarkerSize', 12, 'MarkerFaceColor', 'm', 'DisplayName', '全局起点');

% 轨迹起点（第一个子轨迹起点）
plot3(rotated_traj_start(1), rotated_traj_start(2), rotated_traj_start(3), ...
    'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', '轨迹起点');

% 轨迹终点（最后一个子轨迹终点）
plot3(rotated_traj_end(1), rotated_traj_end(2), rotated_traj_end(3), ...
    'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'DisplayName', '轨迹终点');

% 添加文本标签
text(rotated_global_start(1), rotated_global_start(2), rotated_global_start(3)+5, ...
    '全局起点', 'FontSize', 10, 'Color', 'm', 'FontWeight', 'bold');
text(rotated_traj_start(1), rotated_traj_start(2), rotated_traj_start(3)+5, ...
    '轨迹起点', 'FontSize', 10, 'Color', 'g', 'FontWeight', 'bold');
text(rotated_traj_end(1), rotated_traj_end(2), rotated_traj_end(3)+5, ...
    '轨迹终点', 'FontSize', 10, 'Color', 'b', 'FontWeight', 'bold');

% 添加图例
legend('Location', 'best');

% 添加网格和美观设置
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
box on;

%% 图 10：最终焊接轨迹在点云上的显示
figure(10); clf;
set(gcf, 'Position', [100, 100, 1200, 800]);

pcshow(ptCloud_rotated, 'MarkerSize', 10, 'Background', 'w'); 
hold on;

% 连续轨迹线
plot3(fullTrajectory_closed_rotated(:,1), ...
      fullTrajectory_closed_rotated(:,2), ...
      fullTrajectory_closed_rotated(:,3), ...
      'r-', 'LineWidth', 2);

% 抽点显示 marker
skip = 10;
idx = 1:skip:size(fullTrajectory_closed_rotated,1);

plot3(fullTrajectory_closed_rotated(idx,1), ...
      fullTrajectory_closed_rotated(idx,2), ...
      fullTrajectory_closed_rotated(idx,3), ...
      'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5);

% 起点 / 终点
plot3(fullTrajectory_closed_rotated(1,1), ...
      fullTrajectory_closed_rotated(1,2), ...
      fullTrajectory_closed_rotated(1,3), ...
      'mo', 'MarkerSize', 18, 'MarkerFaceColor', 'm');

plot3(fullTrajectory_closed_rotated(end,1), ...
      fullTrajectory_closed_rotated(end,2), ...
      fullTrajectory_closed_rotated(end,3), ...
      'mo', 'MarkerSize', 18, 'MarkerFaceColor', 'm');

title('最终焊接轨迹在点云上的显示', 'FontSize', 16);
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
legend({'点云', '焊接轨迹', '轨迹采样点', '起点 / 终点'}, 'Location', 'best');

grid on; axis equal; view(-30, 30);
saveas(gcf, fullfile(output_dir, '10_final_trajectory_on_point_cloud.png'));

%% 重采样函数
function traj_resampled = resampleTrajectory(traj, numPoints)
    if size(traj, 1) < 2
        traj_resampled = traj;
        return;
    end
    dist_cum = [0; cumsum(sqrt(sum(diff(traj).^2, 2)))];
    dist_norm = dist_cum / dist_cum(end);
    [dist_norm_unique, ia] = unique(dist_norm, 'stable');
    traj_unique = traj(ia, :);
    dist_uniform = linspace(0, 1, numPoints)';
    traj_resampled = interp1(dist_norm_unique, traj_unique, dist_uniform, 'linear', 'extrap');
end