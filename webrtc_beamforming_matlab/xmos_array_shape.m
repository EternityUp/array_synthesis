clear all
close all
clc

M = 6; % 阵元个数
r = 43.37e-2; % 圆阵半径
theta = 360 / M; % 相邻阵元的圆心角
array_loc = zeros(M,3);
for i = 1 : M
    array_loc(i,1) = r * cosd( ( i - 1 ) * theta );
    array_loc(i,2) = r * sind( ( i - 1 ) * theta );
end

plot3( array_loc(:, 1), array_loc(:, 2), array_loc(:, 3), 'bo', ...
    'markerfacecolor', 'g')
axis tight


array_loc_xmos = zeros(M-2,3);
array_loc_xmos(1,:) = array_loc(1,:);  % mic1 
array_loc_xmos(2,:) = array_loc(5,:);  % mic3
array_loc_xmos(3,:) = array_loc(4,:);  % mic4 
array_loc_xmos(4,:) = array_loc(2,:);  % mic6
save array_loc_xmos.mat array_loc_xmos






