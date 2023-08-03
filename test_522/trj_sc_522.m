glvs
ts = 0.01;       % sampling interval
avp0 = [[0;0;0]; [0;0;0]; glv.pos0]; % init avp
% trajectory segment setting 轨迹如下：  
% 表示协调左转弯，先横滚方向左转4s，再整体转45s,最后横滚方向右转4s.考虑到了航向与横滚
%1.直线运动(27s)：  静止20s;     5m/s^2加速度加速5s;    25m/s速度匀速直线2s;
%2."8"字形圆周运动（52s)：表示协调左转弯，先横滚方向左转4s，20度/s 转角变化向左圆周运动，最后横滚方向右转4s.  
%（考虑到了航向与横滚）表示协调右转弯，先横滚方向右转4s，20度/s 转角变化向右圆周运动，最后横滚方向右转右s.
%3.加速与上下坡：25m/s速度匀速直线2s,7m/s^2加速度加速3s;46m/s匀速直线运动2s;5度上坡5s，持续2s，5度下坡5s;
%4.拐弯运动：
%5.
xxx = [];
seg = trjsegment(xxx, 'init',         0);   %表示轨迹结构数组的初始化
seg = trjsegment(seg, 'uniform',      20);  %表示保持上一状态20s
seg = trjsegment(seg, 'accelerate',   5, xxx, 5);
seg = trjsegment(seg, 'uniform',      2); 
seg = trjsegment(seg, '8turn', [], 20, [], 4);
% seg = trjsegment(seg, '8turn', [], w, [], rolllasting); 
seg = trjsegment(seg, 'uniform',      2);
seg = trjsegment(seg, 'accelerate',   7, xxx, 3);
seg = trjsegment(seg, 'uniform',      2);
seg = trjsegment(seg, 'climb',        5, 5, xxx, 2);
seg = trjsegment(seg, 'uniform',      1);
seg = trjsegment(seg, 'coturnright', 3, 30, [], 1); 
seg = trjsegment(seg, 'uniform',      1);
seg = trjsegment(seg, 'turnright', 3, 30);
% seg = trjsegment(seg, 'coturnright', 3, 30, [], 1); 
seg = trjsegment(seg, 'uniform',      1);
seg = trjsegment(seg, 'descent', 2, 30, [], 2);
seg = trjsegment(seg, 'uniform',      1);
seg = trjsegment(seg, 'deaccelerate', 3, [], 10);
seg = trjsegment(seg, 'uniform',      1);
seg = trjsegment(seg, 'coturnright', 3, 30, [], 1); 
seg = trjsegment(seg, 'coturnleft', 3, 30, [], 1);
seg = trjsegment(seg, 'accelerate',   3, xxx, 12);
seg = trjsegment(seg, 'accelerate',   3, xxx, 15);
seg = trjsegment(seg, 'uniform',      5);
seg = trjsegment(seg, 'turnright', 3, 30);
seg = trjsegment(seg, 'deaccelerate', 1, [], 18);
seg = trjsegment(seg, 'uniform',      1);

trj_sc522  = trjsimu(avp0, seg.wat, ts, 1);
trjfile('trj_sc522.mat', trj_sc522);
insplot(trj_sc522.avp);
imuplot(trj_sc522.imu);

%%
figure(101); scatter3(trj_sc522.avp(:, 7), trj_sc522.avp(:, 8), trj_sc522.avp(:, 9));
title("trajectory in ENU frame");

