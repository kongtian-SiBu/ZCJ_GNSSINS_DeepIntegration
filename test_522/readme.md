# 文件说明

( 注：需提前安装好PSINS，可直接百度搜索PSINS并下载 )

## 1. test_522.m
该文件用于生成接收机的运行轨迹，可获取载体的位置(BLH)、速度(ENU)和姿态信息(ENU, 右前上)(国际单位)，并根据轨迹反演得到IMU的数据真值(使用所谓的真值进行纯惯导推算得到的误差极小)。轨迹长度为120秒。

test_522.txt中存储了载体的位置、速度、加速度和加加速度信息，不过是在ECEF坐标系下的。将test_522.txt拷贝到导航信号模拟器(由实验室提供)中，可生成与轨迹相对应的GNSS信号。本实验仅使用了GPS L1CA信号。如此一来IMU和GPS时间天然对准了且IMU无钟差钟漂。

GNSS数据的采样率、中频频率、数据存储格式是三个最重要的参数，要根据实际情况而定。

## 2. test_522.RSIM(M2B1-GPS_L1)GnssObs(20230522-1550).dat.txt
该文件是由模拟器输出的，里面包含了给定时刻的伪距、伪距率和载波相位等真实值。可作为参考来测试程序的正确性。
