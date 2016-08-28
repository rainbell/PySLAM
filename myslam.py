# encoding=utf8
import math
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

class MySlam:
    dmax = 200
    tmax = 5
    rmax = 5000
    rhothreshold = 5000
    pcntthreshold = 50
    pdisthreshold = 100
    angper = 360.0/1024
    errcontrl = [50, 5]
    robotpos = [0, 0, 0]

    def __init__(self):
        self.orgrho = []
        self.orgtheta = []
        self.rho = []
        self.theta = []
        self.brkrho = []
        self.brktheta = []
        self.brkcnt = 0
        self.seprho = []
        self.septheta = []
        self.linecnt = 0
        self.glbline = []
        self.glbrcdata = []
        self.fittedline = []
        self.fittedrcdata = []
        self.matchedline = []
        self.stepit = 0

    def get_feature(self):
        self.rho_filtration()
        self.break_rho()
        self.break_polyline()
        self.fit_line()
        self.get_fittedrcdata()
        return

    def line_match(self):
        self.matchedline = []
        glblinesize = len(self.glbline)
        fittedlinesize = len(self.fittedline)
        matchedit = 0
        for i in range(0, glblinesize):
            dismin = -1
            itmin = -1
            for j in range(0, fittedlinesize):
                distmp = math.sqrt((self.fittedrcdata[j][1]-self.glbrcdata[i][1])**2 + (self.fittedrcdata[j][2]-self.glbrcdata[i][2])**2)
                if dismin < 0 or distmp < dismin:
                    dtrho = abs(self.fittedrcdata[j][0]-self.glbrcdata[i][0])
                    dttheta = abs(self.fittedrcdata[j][3]-self.glbrcdata[i][3])
                    dttheta2 = abs(abs(self.fittedrcdata[j][3]-self.glbrcdata[i][3])-360.0)
                    dttheta = min(dttheta, dttheta2)
                    #dttheta = min(abs(self.fittedrcdata[j][3]-self.glbrcdata[i][3]), abs(self.fittedrcdata[j][3]-self.glbrcdata[i][3]))
                    if dtrho < self.errcontrl[0] and dttheta < self.errcontrl[1]:
                        dismin = distmp
                        itmin = j
            if itmin > -1:
                self.matchedline.append([i, itmin])
                matchedit += 1
        if matchedit == 1:
            checkwhy = 1
        return

    def renew_robot(self):
        matchedlinesize = len(self.matchedline)
        dxythetait = 0
        dxythetavec = []
        for i in range(0, matchedlinesize):
            tA1 = self.glbline[self.matchedline[i][0]][0]
            tB1 = self.glbline[self.matchedline[i][0]][1]
            for j in range(i+1, matchedlinesize):
                tA2 = self.glbline[self.matchedline[j][0]][0]
                tB2 = self.glbline[self.matchedline[j][0]][1]
                touterproduct = abs(tA1*tB2-tA2*tB1)/(math.sqrt(tA1**2+tB1**2)*math.sqrt(tA2**2+tB2**2))
                if touterproduct < 0.5:
                    continue
                i1 = self.matchedline[i][0]
                i2 = self.matchedline[j][0]
                j1 = self.matchedline[i][1]
                j2 = self.matchedline[j][1]
                dx,dy,dtheta = self.cal_coortranspara(i1,i2,j1,j2)
                dxythetavec.append([dx, dy, dtheta])
                dxythetait += 1
                if False:
                    self.robotpos[0] = dx
                    self.robotpos[1] = dy
                    self.robotpos[2] = dtheta
                    return
        l1 = [tmp[0] for tmp in dxythetavec]
        l2 = [tmp[1] for tmp in dxythetavec]
        self.robotpos[0] = np.mean(np.array([tmp[0] for tmp in dxythetavec]))
        self.robotpos[1] = np.mean(np.array([tmp[1] for tmp in dxythetavec]))
        self.robotpos[2] = np.mean(np.array([tmp[2] for tmp in dxythetavec]))
        return

    def trans_feature(self):
        fittedlinesize = len(self.fittedline)
        for i in range(0, fittedlinesize):
            xai = self.fittedline[i][3]
            yai = self.fittedline[i][4]
            xaii = self.fittedline[i][5]
            yaii = self.fittedline[i][6]

            xbi = xai*math.cos(self.robotpos[2]*math.pi/180.0) - yai*math.sin(self.robotpos[2]*math.pi/180.0)+self.robotpos[0]
            ybi = xai*math.sin(self.robotpos[2]*math.pi/180.0) + yai*math.cos(self.robotpos[2]*math.pi/180.0)+self.robotpos[1]
            xbii = xaii*math.cos(self.robotpos[2]*math.pi/180.0) - yaii*math.sin(self.robotpos[2]*math.pi/180.0)+self.robotpos[0]
            ybii = xaii*math.sin(self.robotpos[2]*math.pi/180.0) + yaii*math.cos(self.robotpos[2]*math.pi/180.0)+self.robotpos[1]
            tA = (ybii-ybi)
            tB = -(xbii-xbi)
            tC = -xbi*(ybii-ybi)+ybi*(xbii-xbi)
            self.fittedline[i][0] = tA
            self.fittedline[i][1] = tB
            self.fittedline[i][2] = tC
            self.fittedline[i][3] = xbi
            self.fittedline[i][4] = ybi
            self.fittedline[i][5] = xbii
            self.fittedline[i][6] = ybii
        matchedlinesize = len(self.matchedline)
        if True:
            for i in range(0, matchedlinesize):
                Mi = self.matchedline[i][0]
                Mj = self.matchedline[i][1]
                xsi = self.glbline[Mi][3]
                ysi = self.glbline[Mi][4]
                xei = self.glbline[Mi][5]
                yei = self.glbline[Mi][6]
                xsii = self.fittedline[Mj][3]
                ysii = self.fittedline[Mj][4]
                xeii = self.fittedline[Mj][5]
                yeii = self.fittedline[Mj][6]
                lsi = math.sqrt((xei-xsii)**2+(yei-ysii)**2)
                lsii = math.sqrt((xeii-xsii)**2+(yeii-ysii)**2)
                if lsi > lsii:
                    self.fittedline[Mj][5] = xei
                    self.fittedline[Mj][6] = yei
                lei = math.sqrt((xeii-xsi)**2+(yeii-ysi)**2)
                leii = math.sqrt((xeii-xsii)**2+(yeii-ysii)**2)
                if lei > leii:
                    self.fittedline[Mj][3] = xsi
                    self.fittedline[Mj][4] = ysi
        self.glbline = []
        self.glbrcdata = []
        for i in range(0, fittedlinesize):
            self.glbline.append(self.fittedline[i])
            self.glbrcdata.append(self.fittedrcdata[i])
        self.fittedline = []
        self.fittedrcdata = []
        return

    def draw_orgdata(self, ms, ax):
        ax.cla()
        theta = np.array(ms.orgtheta)*math.pi/180
        rho = np.array(ms.orgrho)
        ax.plot(theta, rho, 'b+', linewidth=1)
        return

    def draw_feature(self, ms, ax):
        fittedlinesize = len(ms.fittedline)
        for i in range(0, fittedlinesize):
            tmplinepara = ms.fittedline[i]
            rhotmp = []
            thetatmp = []
            tx1 = tmplinepara[3]
            ty1 = tmplinepara[4]
            rhotmp.append(math.sqrt(tx1**2+ty1**2))
            if tx1 >= 0 and ty1 >= 0:
                thetatmp.append(math.asin(ty1/rhotmp[0]))
            elif tx1 < 0 and ty1 >= 0:
                thetatmp.append(math.pi-math.asin(ty1/rhotmp[0]))
            elif tx1 < 0 and ty1 < 0:
                thetatmp.append(math.pi-math.asin(ty1/rhotmp[0]))
            else:
                thetatmp.append(2*math.pi+math.asin(ty1/rhotmp[0]))
            tx2 = tmplinepara[5]
            ty2 = tmplinepara[6]
            rhotmp.append(math.sqrt(tx2**2+ty2**2))
            if tx2 >= 0 and ty2 >= 0:
                thetatmp.append(math.asin(ty2/rhotmp[1]))
            elif tx2 < 0 and ty2 >= 0:
                thetatmp.append(math.pi-math.asin(ty2/rhotmp[1]))
            elif tx2 < 0 and ty2 < 0:
                thetatmp.append(math.pi-math.asin(ty2/rhotmp[1]))
            else:
                thetatmp.append(2*math.pi+math.asin(ty2/rhotmp[1]))
            ax.plot(thetatmp,rhotmp,'r-',linewidth=1)
        return

    def draw_map(self, ms, ax):
        glblinesize = len(ms.glbline)
        for i in range(0, glblinesize):
            ax.plot([ms.glbline[i][3],ms.glbline[i][5]],[ms.glbline[i][4],ms.glbline[i][6]],'b-',linewidth=3)
            ax.plot(ms.robotpos[0],ms.robotpos[1],'r*')
        return

    def rho_filtration(self):
        self.rho = []
        self.theta = []
        orgrhosize = len(self.orgrho)
        rhoit=0
        for i in range(0, orgrhosize):
            if self.orgrho[i] < self.rmax:
                self.rho.append(self.orgrho[i])
                self.theta.append(self.orgtheta[i])
                rhoit += 1
        return

    def break_rho(self):
        self.brkrho = []
        self.brktheta = []
        lastrho = self.rho[0]
        lasttheta = self.theta[0]
        rhosize=len(self.rho)
        self.brkrho.append(lastrho)
        self.brktheta.append(lasttheta)
        brkit=1
        brkcnt=1
        for i in range(1, rhosize):
            trho = self.rho[i]
            ttheta = self.theta[i]
            dis = abs(trho - lastrho)
            dtheta = abs(ttheta - lasttheta)
            if dis>=self.dmax or dtheta>=self.tmax:
                self.brkrho.append(-1)
                self.brktheta.append(1000.0)
                brkit += 1
                brkcnt += 1
            self.brkrho.append(trho)
            self.brktheta.append(ttheta)
            brkit += 1
            lastrho = trho
            lasttheta = ttheta
        self.brkrho.append(-1)
        self.brktheta.append(1000.0)
        brkit += 1
        return

    def break_polyline(self):
        self.seprho = []
        self.septheta = []
        pointcnt = 0
        linecnt = 0
        X = []
        Y = []
        rhocopy = []
        thetacopy = []
        brkrhosize = len(self.brkrho)
        sepit=0
        for i in range(0, brkrhosize):
            trho = self.brkrho[i]
            ttheta = self.brktheta[i]
            if trho < 0:
                if pointcnt > self.pcntthreshold:
                    cornercnt = 0
                    cornerindex = []
                    self.find_corners(cornerindex, X, Y, 0, pointcnt, self.pdisthreshold)
                    #sorted(cornerindex)
                    cornercnt = len(cornerindex)
                    if cornercnt > 1:
                        cornerindex.sort()
                    if cornercnt == 0:
                        linecnt += 1
                        for j in range(0, pointcnt):
                            self.seprho.append(rhocopy[j])
                            self.septheta.append(thetacopy[j])
                            sepit += 1
                        self.seprho.append(-1)
                        self.septheta.append(1000.0)
                        sepit += 1
                    else:
                        tmpindex = 0
                        for j in range(0, pointcnt):
                            self.seprho.append(rhocopy[j])
                            self.septheta.append(thetacopy[j])
                            sepit += 1
                            if j == cornerindex[tmpindex]:
                                self.seprho.append(-1)
                                self.septheta.append(1000.0)
                                sepit += 1
                                linecnt += 1
                                if tmpindex < cornercnt-1:
                                    tmpindex += 1
                        self.seprho.append(-1)
                        self.septheta.append(1000.0)
                        sepit += 1
                        linecnt += 1
                X = []
                Y = []
                pointcnt = 0
                rhocopy = []
                thetacopy = []
            else:
                X.append(trho*math.cos(ttheta*math.pi/180.0))
                Y.append(trho*math.sin(ttheta*math.pi/180.0))
                pointcnt += 1
                rhocopy.append(trho)
                thetacopy.append(ttheta)

    def fit_line(self):
        self.fittedline = []
        X=[]
        Y=[]
        pointcnt = 0
        seprhosize = len(self.seprho)
        fittedit = 0
        for i in range(0, seprhosize):
            trho = self.seprho[i]
            ttheta = self.septheta[i]
            if trho < 0:
                if pointcnt < 20:
                    pointcnt=0
                    X = []
                    Y = []
                    continue
                tmplinepara = [0]*7
                if max(X)-min(X) < 100:
                    tmplinepara[0] = (Y[pointcnt-1]-Y[0])
                    tmplinepara[1] = -(X[pointcnt-1]-X[0])
                    tmplinepara[2] = -X[0]*(Y[pointcnt-1]-Y[0]) + Y[0]*(X[pointcnt-1]-X[0])	#X+C=0
                    tmplinepara[3] = X[0]
                    tmplinepara[4] = Y[0]
                    tmplinepara[5] = X[pointcnt-1]
                    tmplinepara[6] = Y[pointcnt-1]
                else:
                    npX = np.array(X)
                    npY = np.array(Y)
                    p = np.polyfit(npX, npY, 1)
                    tmplinepara[0] = p[0]
                    tmplinepara[1] = -1
                    tmplinepara[2] = p[1]	                                #kX-Y+b=0
                    tmplinepara[3] = X[0]
                    tmplinepara[4] = tmplinepara[0]*X[0]+tmplinepara[2]
                    tmplinepara[5] = X[pointcnt-1]
                    tmplinepara[6] = tmplinepara[0]*X[pointcnt-1]+tmplinepara[2]
                self.fittedline.append(tmplinepara)
                fittedit += 1
                pointcnt = 0
                X = []
                Y = []
            else:
                X.append(trho*math.cos(ttheta*math.pi/180.0))
                Y.append(trho*math.sin(ttheta*math.pi/180.0))
                pointcnt += 1

    def find_corners(self, cornerindex, X, Y, pointsrt, pointcnt, eps):
        maxdisind = self.poly_contourfit(X, Y, pointcnt, eps)
        if maxdisind == 0:
            return
        else:
            cornerindex.append(pointsrt + maxdisind)
            self.find_corners(cornerindex, X[0:maxdisind], Y[0:maxdisind], pointsrt, maxdisind, eps)
            self.find_corners(cornerindex, X[maxdisind:pointcnt], Y[maxdisind:pointcnt], pointsrt+maxdisind, pointcnt-maxdisind, eps)

    def poly_contourfit(self, x, y, n, eps):
        if n == 1:
            return 0
        dis = math.sqrt((x[0]-x[n-1])**2+(y[0]-y[n-1])**2)
        costheta = (x[n-1]-x[0])/dis
        sintheta = -(y[n-1]-y[0])/dis
        maxdis = 0
        maxdisind = -1
        for i in range(0, n):
            dbdis = abs((y[i]-y[0])*costheta+(x[i]-x[0])*sintheta)
            if dbdis > maxdis:
                maxdis = dbdis
                maxdisind = i
        if maxdis > eps:
            return maxdisind
        return 0

    def get_fittedrcdata(self):
        self.fittedrcdata = []
        fittedlinesize = len(self.fittedline)
        for i in range(0, fittedlinesize):
            tA = self.fittedline[i][0]
            tB = self.fittedline[i][1]
            tC = self.fittedline[i][2]
            trcdata = [0]*4
            trcdata[0] = abs(tC/math.sqrt(tA**2+tB**2))
            tx0 = -(tA*tC)/(tA**2+tB**2)
            ty0 = -(tB*tC)/(tA**2+tB**2)
            trcdata[1] = tx0
            trcdata[2] = ty0
            if tx0 >= 0 and ty0 >= 0:
                trcdata[3] = math.asin(ty0/math.sqrt(tx0**2+ty0**2))/math.pi*180.0
            elif tx0 <0 and ty0 >= 0:
                trcdata[3] = 180.0-math.asin(ty0/math.sqrt(tx0**2+ty0**2))/math.pi*180.0
            elif tx0 < 0 and ty0 <= 0:
                trcdata[3] = 180.0-math.asin(ty0/math.sqrt(tx0**2+ty0**2))/math.pi*180.0
            else:
                trcdata[3] = 360.0+math.asin(ty0/math.sqrt(tx0**2+ty0**2))/math.pi*180.0
            self.fittedrcdata.append(trcdata)

    def cal_coortranspara(self, i1, i2, j1, j2):
        dtheta1 = self.glbrcdata[i1][3]-self.fittedrcdata[j1][3]
        if dtheta1 > self.errcontrl[1]:
            dtheta1 = dtheta1 - 360.0
        elif dtheta1 < -self.errcontrl[1]:
            dtheta1 = dtheta1 + 360.0
        dtheta2 = self.glbrcdata[i2][3] - self.fittedrcdata[j2][3]
        if dtheta2 > self.errcontrl[1]:
            dtheta2 = dtheta2 - 360.0
        elif dtheta2 < -self.errcontrl[1]:
            dtheta2 = dtheta2 + 360.0
        dtheta = self.robotpos[2] + (dtheta1 + dtheta2)/2
        tA = [0]*2
        tB = [0]*2
        tC = [0]*2
        tA[0] = self.glbline[i1][0]
        tB[0] = self.glbline[i1][1]
        tC[0] = self.glbline[i1][2]
        tA[1] = self.glbline[i2][0]
        tB[1] = self.glbline[i2][1]
        tC[1] = self.glbline[i2][2]
        Xw =  (tC[1]*tB[0]-tC[0]*tB[1])/(tA[0]*tB[1]-tA[1]*tB[0])
        Yw = -(tC[1]*tA[0]-tC[0]*tA[1])/(tA[0]*tB[1]-tA[1]*tB[0])
        tA[0] = self.fittedline[j1][0]
        tB[0] = self.fittedline[j1][1]
        tC[0] = self.fittedline[j1][2]
        tA[1] = self.fittedline[j2][0]
        tB[1] = self.fittedline[j2][1]
        tC[1] = self.fittedline[j2][2]
        Xr =  (tC[1]*tB[0]-tC[0]*tB[1])/(tA[0]*tB[1]-tA[1]*tB[0])
        Yr = -(tC[1]*tA[0]-tC[0]*tA[1])/(tA[0]*tB[1]-tA[1]*tB[0])
        dx = Xw - math.cos(dtheta*math.pi/180.0)*Xr + math.sin(dtheta*math.pi/180)*Yr
        dy = Yw - math.sin(dtheta*math.pi/180.0)*Xr - math.cos(dtheta*math.pi/180)*Yr
        return dx, dy, dtheta