% SL_HW5 V2_19示例数据处理
% V0.0.1 20220831 tkalpha
clear;clc

a = readmatrix("dtt1_a.txt");
b = readmatrix("dtt1_b.txt");
c = readmatrix("dtt1_c.txt");
d = readmatrix("dtt1_d.txt");
s = readmatrix("dtt1_s.txt");
q = readmatrix("dtt1_q.txt");
qv = readmatrix("dtt1_qv.txt");
qa = readmatrix("dtt1_qa.txt");

% curvature velocity limit
curve = zeros(length(q), 1);
for i = 1 : length(q)
    j = max(min(i, length(q)-1),2);
    curve_dist = [-norm(q(j,:)-q(j-1,:));0;norm(q(j,:)-q(j+1,:))];
    curve_mat = [ones(3,1), curve_dist, curve_dist.^2];
    curve_x = [q(j-1,1); q(j,1); q(j+1,1)];
    curve_y = [q(j-1,2); q(j,2); q(j+1,2)];
    curve_px = curve_mat \ curve_x;
    curve_py = curve_mat \ curve_y;
    curve(i) = -2.0*(curve_px(3)*curve_py(2)-curve_py(3)*curve_px(2))...
        /(curve_px(2)^2+curve_py(2)^2)^(3/2);
end
acc_max = 1.0;
v_curve_max = sqrt(acc_max./abs(curve));

% origin method velocity limit
v_origin_max = ones(length(q), 1)*1e20;
v_origin_min = zeros(length(q), 1);
for i = 1:length(q)-1
    if(qa(i,1) > 0)
        v_origin_max(i) = min(v_origin_max(i), sqrt(max((acc_max-qv(i,1)*a(i))/qa(i,1),0)));
    else
        v_origin_min(i) = max(v_origin_min(i), sqrt(max((acc_max-qv(i,1)*a(i))/qa(i,1),0)));
    end
    if(qa(i,2) > 0)
        v_origin_max(i) = min(v_origin_max(i), sqrt(max((acc_max-qv(i,2)*a(i))/qa(i,2),0)));
    else
        v_origin_min(i) = max(v_origin_min(i), sqrt(max((acc_max-qv(i,2)*a(i))/qa(i,2),0)));
    end
end


figure(220831001); clf(220831001);
plot(q(:,1), q(:,2), '-b');
axis equal; grid on; xlabel('X (m)'); ylabel('Y (m)');
title('global trajectory');

figure(220831002); clf(220831002);
plot(s, v_origin_max, '-', 'LineWidth', 4, 'Color', [0.5,0.5,0.5]);
hold on;
plot(s, v_origin_min, ':', 'LineWidth', 3.5, 'Color', [0.5,0.5,0.5]);
plot(s, sqrt(b), '-b', 'LineWidth', 2.0);
grid on; xlabel('trajectory length (m)'); ylabel('line speed (m/s)');
ylim([0, 7]);
yyaxis right;
plot(s(1:end-1), a, '-r', 'LineWidth', 1.5);
ylabel('line acceleration (m/s^2)');
legend({'Upper speed bound','Lower speed bound','Planned speed', 'Planned acceleration'});