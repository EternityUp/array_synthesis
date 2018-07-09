% 根据阵型仿真阵列形状
clear all
close all
clc

%% 均匀线阵
M = 6; % 阵元个数
d = 0.05; % 阵元间距
array_loc = zeros(M,3);
for i = 1 : M
    array_loc(i,1) = ( i - 1 ) * d;
end

plot3( array_loc(:, 1), array_loc(:, 2), array_loc(:, 3), 'b-o', ...
    'markerfacecolor', 'g')
axis tight

%% 均匀圆阵
M = 6; % 阵元个数
r = 1; % 圆阵半径
theta = 360 / M; % 相邻阵元的圆心角
array_loc = zeros(M,3);
for i = 1 : M
    array_loc(i,1) = r * cosd( ( i - 1 ) * theta );
    array_loc(i,2) = r * sind( ( i - 1 ) * theta );
end

plot3( array_loc(:, 1), array_loc(:, 2), array_loc(:, 3), 'bo', ...
    'markerfacecolor', 'g')
axis tight









