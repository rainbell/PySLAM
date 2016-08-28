class Map:
    def __init__(self):
        self.nodenum = 0
        self.nodevec = []
        self.linenum = 0
        self.linevec = []
        self.wpnum = 0
        self.wpstep = []
        self.wpvec = []
        self.pos = []
        self.allstep = 0

    def load_map(self, filename):
        fin = open(filename, 'rt')
        tnodenum = map(int, fin.readline().split())
        self.nodenum = tnodenum[0]
        for i in range(0, self.nodenum):
            self.nodevec.append(map(float, fin.readline().split()))
        tlinenum = map(int, fin.readline().split())
        self.linenum = tlinenum[0]
        for i in range(0, self.linenum):
            self.linevec.append(map(int, fin.readline().split()))
        twpnum = map(int, fin.readline().split())
        self.wpnum = twpnum[0]
        for i in range(0, self.wpnum):
            twpstep = map(int, fin.readline().split())
            self.wpstep.append(twpstep[0])
            self.wpvec.append(map(float, fin.readline().split()))
        fin.close()
        posnum = 0
        for i in range(0, self.wpnum-1):
            tstep = self.wpstep[i+1] - self.wpstep[i]
            for j in range(0, tstep):
                tpos = [0]*3
                tpos[0] = self.wpvec[i][0] + (self.wpvec[i+1][0] - self.wpvec[i][0])/tstep*j
                tpos[1] = self.wpvec[i][1] + (self.wpvec[i+1][1] - self.wpvec[i][1])/tstep*j
                tpos[2] = self.wpvec[i][2] + (self.wpvec[i+1][2] - self.wpvec[i][2])/tstep*j
                self.pos.append(tpos)
            self.allstep += tstep
        return