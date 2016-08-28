import numpy as np
import math

class SimuRadar:
    def __init__(self):
        self.angvec = []
        self.rho = []
        self.theta = []
        return

    def runradar(self, mp, ms):
        ms.orgrho = []
        ms.orgtheta = []
        tx0 = mp.pos[ms.stepit][0]
        ty0 = mp.pos[ms.stepit][1]
        angstart = mp.pos[ms.stepit][2]-30.0
        angend = angstart+240.0
        angnow = angstart
        angnum=0
        while angnow <= angend:
            self.angvec.append(angnow)
            angnum += 1
            angnow += ms.angper
        for j in range(0, angnum):
            ttheta = self.angvec[j] + np.random.normal(0, 0.01, 1)[0]
            tdis = ms.rmax
            cpnum = 0
            cpx = []
            cpy = []
            cplen = []
            for k in range(0, mp.linenum):
                tx1 = mp.nodevec[mp.linevec[k][0]-1][0]
                ty1 = mp.nodevec[mp.linevec[k][0]-1][1]
                tx2 = mp.nodevec[mp.linevec[k][1]-1][0]
                ty2 = mp.nodevec[mp.linevec[k][1]-1][1]
                flag, tcpx, tcpy = self.get_ray_line_cp(tx0, ty0, ttheta, tx1, ty1, tx2, ty2)
                if flag == 1:
                    cpx.append(tcpx)
                    cpy.append(tcpy)
                    cplen.append(math.sqrt((cpx[cpnum]-tx0)**2+(cpy[cpnum]-ty0)**2))
                    cpnum += 1
            if cpnum > 0:
                tdis = min(cplen)
            tdis = tdis+np.random.normal(0, 10, 1)[0]
            ms.orgrho.append(tdis)
            ms.orgtheta.append(ttheta - mp.pos[ms.stepit][2])
        return

    @staticmethod
    def get_ray_line_cp(x0, y0, theta, x1, y1, x2, y2):
        pi = 3.1415926535
        rmax = 5000
        flag = 0
        cpx = rmax*math.cos(theta*pi/180.0)
        cpy = rmax*math.sin(theta*pi/180.0)
        if abs(theta-90.0) < 1e-3:
            if abs(x2-x1) < 1e-3:
                flag = 0
                return flag, cpx, cpy
            elif x0 < min(x1, x2) or x0 > max(x1, x2):
                flag = 0
                return flag, cpx, cpy
            else:
                cpx = x0
                cpy = y1+(y2-y1)/(x2-x1)*(x0-x1)
                tm = (cpx-x0)*math.cos(theta*pi/180)+(cpy-y0)*math.sin(theta*pi/180)
                if tm <= 0:
                    flag = 0
                    return flag, cpx, cpy
                else:
                    flag = 1
                    return flag, cpx, cpy
        else:
            tk1 = math.tan(theta*pi/180)
            if abs(x2-x1) < 1e-3:
                cpx = x1
                cpy = tk1*(x1-x0)+y0
                if cpy < min(y1, y2) or cpy > max(y1, y2):
                    flag = 0
                    return flag, cpx, cpy
                else:
                    tm = (cpx-x0)*math.cos(theta*pi/180)+(cpy-y0)*math.sin(theta*pi/180)
                    if tm <= 0:
                        flag = 0
                        return flag, cpx, cpy
                    else:
                        flag = 1
                        return flag, cpx, cpy
            else:
                tk2 = (y2-y1)/(x2-x1)
                cpx = ((tk1*x0-y0)-(tk2*x1-y1))/(tk1-tk2)
                cpy = tk1*(cpx-x0)+y0
                tm = (cpx-x0)*math.cos(theta*pi/180)+(cpy-y0)*math.sin(theta*pi/180)
                if tm <= 0:
                    flag = 0
                    return flag, cpx, cpy
                tm = (x1-cpx)*(x2-cpx)+(y1-cpy)*(y2-cpy)
                if tm >= 0:
                    flag = 0
                    return flag, cpx, cpy
                flag = 1
                return flag, cpx, cpy