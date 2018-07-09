% �������ͷ���������״
clear all
close all
clc

%% ��������
M = 6; % ��Ԫ����
d = 0.05; % ��Ԫ���
array_loc = zeros(M,3);
for i = 1 : M
    array_loc(i,1) = ( i - 1 ) * d;
end

plot3( array_loc(:, 1), array_loc(:, 2), array_loc(:, 3), 'b-o', ...
    'markerfacecolor', 'g')
axis tight

%% ����Բ��
M = 6; % ��Ԫ����
r = 1; % Բ��뾶
theta = 360 / M; % ������Ԫ��Բ�Ľ�
array_loc = zeros(M,3);
for i = 1 : M
    array_loc(i,1) = r * cosd( ( i - 1 ) * theta );
    array_loc(i,2) = r * sind( ( i - 1 ) * theta );
end

plot3( array_loc(:, 1), array_loc(:, 2), array_loc(:, 3), 'bo', ...
    'markerfacecolor', 'g')
axis tight









