import loadmap
import myslam
import simuradar
import matplotlib.pyplot as plt

mp = loadmap.Map()
mp.load_map('map.dat')
rd = simuradar.SimuRadar()
ms = myslam.MySlam()
plt.ion()
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111, polar=True)
fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111)
for ms.stepit in range(0, mp.allstep):
    rd.runradar(mp, ms)
    ms.draw_orgdata(ms, ax1)
    #plt.pause(0.001)
    ms.get_feature()
    ms.draw_feature(ms, ax1)
    #plt.pause(0.001)
    if ms.stepit > 0:
        ms.line_match()
        ms.renew_robot()
        print ms.stepit, ms.robotpos
    ms.trans_feature()
    ms.draw_map(ms, ax2)
    plt.pause(0.001)
    continue
