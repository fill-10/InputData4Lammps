
# bead.py
# bead
class bead:
    def __init__(self, **kwargs):
        self.name = 'H'
        self.coor = [0.,0.,0.]
        self.tpnumber = 1
        self.mass = 1.0
        for kw in kwargs.items():
            if kw[0] == 'mass':
                self.mass = kw[1]
            elif kw[0] == 'coor':
                self.coor = kw[1]
            elif kw[0] == 'typenumber':
                self.tpnumber = kw[1]
            elif kw[0] == 'name':
                self.name = kw[1]



if __name__ == '__main__':
    bd1 = bead(mass = 10, coor = [1.,1.,1.])
    print( bd1.mass)
    print( bd1.coor)
