'''
by sdust yuzhigang 20210.09.28 
'''
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
class MineParmeters:
    def __init__(self):
        '''
        读取参数
        :return:
        '''
        self.parsDict = {}  
        self.matrix=None 
        self.matriy=None 
        self.matriX=None 
        self.matriY=None 
        filename = 'pars.txt' 
        with open(filename, 'r') as file_to_read:
            while True:
                lines = file_to_read.readline() 
                if not lines:
                    break
                    pass
                p_tmp, E_tmp, C_tmp= [str(i) for i in lines.split('\t')]#
                self.parsDict[p_tmp]=float(E_tmp)
                pass
            pass
        print(self.parsDict['q'])
    
        self.parsDict['alpha']=self.parsDict['alpha']/180.0*math.pi 
        self.parsDict['beta'] = self.parsDict['beta'] / 180.0 * math.pi   
        self.parsDict['beta1'] = self.parsDict['beta1'] / 180.0 * math.pi  
        self.parsDict['thet'] = self.parsDict['thet'] / 180.0 * math.pi   
        self.parsDict['thet0'] = self.parsDict['thet0'] / 180.0 * math.pi  
        self.l=self.parsDict['D3']-self.parsDict['s3']-self.parsDict['s4'] 
        self.L=self.parsDict['D1']-self.parsDict['s1']-self.parsDict['s2']
        self.w0 = self.parsDict['m'] * self.parsDict['q'] *math.cos(self.parsDict['alpha']) 
        self.cxm=None
        self.cym=None
        print(self.w0,'w0')
    def BuldMatrix(self): #, r1, r2, r3, r4, l, D1):
        '''
     
        '''
        l=self.parsDict['D3']-self.parsDict['s3']-self.parsDict['s3']
        L=abs((self.parsDict['D1']-self.parsDict['s1']-self.parsDict['s2'])*math.sin(self.parsDict['thet0']+self.parsDict['alpha'])/math.sin(self.parsDict['thet0']))
        print(L,"L")
        matrixXL=round(self.parsDict['r1']*1.5+self.parsDict['r2']*1.5+l)
        matrixYL=round(self.parsDict['r3']*1.5+self.parsDict['r4']*1.5+L)
        #matrix=np.arange(0,matrixXL,1)
        matrixList=[[i for i in range(0,matrixXL)] for j in range(0,matrixYL)]
        matriyList = [[j for i in range(0,matrixXL)] for j in range(0,matrixYL)]
        matrix=np.array(matrixList)
        matriy = np.array(matriyList)
        self.matrix=matrix-round(self.parsDict['r1']*1.5) #?????
        self.matriy = matriy-round(self.parsDict['r3']*1.5)
        return self.matrix,self.matriy
    def caculateSemiInfinite(self, b, r, x):
        '''

        '''
        # w0 = m * q * cos(alpha)  
        Row, Col = x.shape
        wx=np.ones((Row,Col))
        ix=np.ones((Row,Col))
        kx=np.ones((Row,Col))
        ux=np.ones((Row,Col))
        ex=np.ones((Row,Col))
        m_exp=np.ones((Row,Col))
        for i in range(Row):
            for j in range(Col):
                m_exp[i,j] = math.exp(-math.pi * x[i,j]*x[i,j] / r / r)
                wx[i,j] = self.w0 / 2 * (math.erf(math.sqrt(math.pi) / r * x[i,j]) + 1)  #
                ix[i,j] = self.w0 / r * m_exp[i,j]  # 
                kx[i,j] = -(2 * math.pi * self.w0/ r / r / r * x[i,j] * m_exp[i,j])  # 
                ux[i,j] = b * r * ix[i,j]  # 
                ex[i,j] = b * r * kx[i,j]  # 
        return wx, ix, kx, ux, ex


    def caculateLimited1(self,b, r, x):
        '''
   
        '''
        l = self.l
        m_wx, m_ix, m_kx, m_ux, m_ex = MineParmeters.caculateSemiInfinite( self,b, r, x)
        m_wxl, m_ixl, m_kxl, m_uxl, m_exl = MineParmeters.caculateSemiInfinite(self, b, r, (x - l))
        w0x = m_wx - m_wxl
        i0x = m_ix - m_ixl
        k0x = m_kx - m_kxl
        u0x = m_ux - m_uxl
        e0x = m_ex - m_exl
        return w0x, i0x, k0x, u0x, e0x
    def caculateLimited2(self, b1, r1, b2, r2, y):
        '''

        '''
        cot = 1 / math.tan(self.parsDict['thet0'])
        m_wy, m_iy, m_ky, m_uy, m_ey = MineParmeters.caculateSemiInfinite(self, b1, r1, y)
        m_wyL, m_iyL, m_kyL, m_uyL, m_eyL = MineParmeters.caculateSemiInfinite(self, b2, r2, (y - self.L))
        w0y = m_wy - m_wyL
        i0y = m_iy - m_iyL
        k0y = m_ky - m_kyL
        u0y = m_uy + m_wy * cot - m_uyL - m_wyL * cot
        e0y = m_ey + m_iy * cot - m_eyL - m_iyL * cot

        return w0y, i0y, k0y, u0y, e0y
    def caculateLimitedAll(self,x, y):
        '''

        '''
        h_l = self.l / 2
        m_exp = math.exp(-math.pi * h_l * h_l / self.parsDict['r'] / self.parsDict['r'])
        w0x1 = self.w0 / 2 * (math.erf(math.sqrt(math.pi) / self.parsDict['r'] * h_l) + 1)  # 
        h_L = self.L / 2
        m_exp = math.exp(-math.pi * h_L * h_L / self.parsDict['r']  / self.parsDict['r'] )
        w0y1 = self.w0 / 2 * (math.erf(math.sqrt(math.pi) / self.parsDict['r']  * h_L) + 1)  # 
        self.cxm = w0x1 / self.w0  # 
        self.cym = w0y1 / self.w0  # 
        print(self.cxm,self.cym,'self.cxm')

        wxt3, ixt3, kxt3, uxt3, ext3 = MineParmeters.caculateSemiInfinite( self,self.parsDict['b3'], self.parsDict['r3'], x)
        wxlt4, ixlt4, kxlt4, uxlt4, exlt4 = MineParmeters.caculateSemiInfinite( self,self.parsDict['b4'], self.parsDict['r4'], x - self.l)
        wyt1, iyt1, kyt1, uyt1, eyt1 = MineParmeters.caculateSemiInfinite(self, self.parsDict['b1'], self.parsDict['r1'], y)
        wyLt2, iyLt2, kyLt2, uyLt2, eyLt2 = MineParmeters.caculateSemiInfinite(self, self.parsDict['b2'], self.parsDict['r2'], y - self.L)
        w0xAll = self.cym * (wxt3 - wxlt4)
        i0xAll = self.cym * (ixt3 - ixlt4)
        k0xAll = self.cym * (kxt3 - kxlt4)
        u0xAll = self.cym * (uxt3 - uxlt4)
        e0xAll = self.cym * (ext3 - exlt4)

        w0yAll = self.cxm * (wyt1 - wyLt2)
        i0yAll = self.cxm * (iyt1 - iyLt2)
        k0yAll = self.cxm * (kyt1 - kyLt2)
        u0yAll = self.cxm * (uyt1 - uyLt2)
        e0yAll = self.cxm * (eyt1 - eyLt2)

        return  w0xAll, i0xAll, k0xAll, u0xAll, e0xAll,w0yAll, i0yAll, k0yAll, u0yAll, e0yAll
    def caculate_any_xy(self,x, y, phi):
        '''

        '''
        phi=phi/180.0*math.pi  #
        w0xAll, i0xAll, k0xAll, u0xAll, e0xAll, w0yAll, i0yAll, k0yAll, u0yAll, e0yAll = MineParmeters.caculateLimitedAll(self,x, y)

        #wxy = w0xAll* self.cym
        w01=1/self.w0  #
        wxy = w01*w0xAll * w0yAll
        print(wxy.shape,'wxy.shape')
        ixy_phi=w01*(i0xAll*w0yAll*math.cos(phi)+w0xAll*i0yAll*math.sin(phi))
        kxy_phi =w01*(k0xAll*w0yAll*math.cos(phi)*math.cos(phi)+k0yAll*w0xAll*math.sin(phi) * math.sin(
            phi) +i0xAll * i0yAll * math.sin(2 * phi))
        uxy_phi =ixy_phi*self.parsDict['b']*self.parsDict['r']
        exy_phi=kxy_phi*self.parsDict['b']*self.parsDict['r']


        #ixy_phi = i0xAll * self.cym * math.cos(phi) + i0yAll * self.cxm * math.sin(phi)  # iax he  iay??????
        #kxy_phi = k0xAll * self.cym * math.cos(phi) * math.cos(phi) + k0yAll * self.cxm * math.sin(phi) * math.sin(
        #    phi) + i0xAll * i0yAll / self.w0 * math.sin(2 * phi)
        #uxy_phi = u0xAll * self.cym * math.cos(phi) + u0yAll * self.cxm * math.sin(phi)
        #exy_phi = e0xAll * self.cym * math.cos(phi) * math.cos(phi) + e0yAll * self.cxm * math.sin(phi) * math.sin(phi) + (
                    #u0xAll * i0yAll + u0yAll * i0xAll) / self.w0 * math.sin(phi) * math.cos(phi)
        return wxy,ixy_phi,kxy_phi,uxy_phi,exy_phi
    def planimage(self,wx, ix, kx, ux, ex ):
        '''
        :return:
        '''
        plt.figure(figsize=(8, 8))
        plt.subplot(231)
        plot1 = plt.imshow(wx, cmap=plt.cm.jet)
        plt.colorbar()
        plt.title('w-chenjiang')
        plt.subplot(232)
        plot1 = plt.imshow(ix, cmap=plt.cm.jet)
        plt.colorbar()
        plt.title('ix-qingxie')
        plt.subplot(233)
        plot1 = plt.imshow(kx, cmap=plt.cm.jet)
        plt.colorbar()
        plt.title('kx-qulv')
        plt.subplot(234)
        plot1 = plt.imshow(ux, cmap=plt.cm.jet)
        plt.colorbar()
        plt.title('ux-shuipingyidong')
        plt.subplot(235)
        plot1 = plt.imshow(ex, cmap=plt.cm.jet)
        plt.colorbar()
        plt.title('ex-shuipingbianxing')
        plt.subplot(236)
        plot1 = plt.imshow(matrix, cmap=plt.cm.jet)
        plt.colorbar()
        plt.title('x-l')
        plt.show()
    def image3D(self,wx, ix, kx, ux, ex ):
        '''
        :return:
        '''
        x,y=matrix.shape
        xline=self.matrix.reshape(x*y)
        yline=self.matriy.reshape(x*y)
        wxline=wx.reshape(x*y)

        fig = plt.figure()
        ax3 = plt.axes(projection='3d')
        YY = np.arange(-round(self.parsDict['r1'] * 1.5),x-round(self.parsDict['r1'] * 1.5),1)
        XX = np.arange(-round(self.parsDict['r3'] * 1.5), y - round(self.parsDict['r3'] * 1.5), 1)
        # YY = np.arange(-0, x, 1)
        # XX = np.arange(0, y, 1)
        X, Y = np.meshgrid(XX, YY)
        Z = -wx
        minz=np.min(Z)

        ax3.plot_surface(X, Y, Z, cmap='rainbow')
        cset=ax3.contour(X, Y, Z,8,zdim='z', offset=minz, cmap='rainbow')  # 
        #plt.clabel(cset, inline=1, fontsize=12,inline_spacing=0.1,colors='red')
        plt.clabel(cset, inline=1, fontsize=12, colors='red')
        plt.show()

    def coord_trans(self,int):
        '''

        '''

        deltaY=self.parsDict['Y1']-self.parsDict['Y0']
        deltaX=self.parsDict['X1']-self.parsDict['X0']
        if deltaX==0 :
            if deltaY>0 :
                thetxX=0
            elif deltaY<0:
                thetxX = math.pi
        elif deltaY==0 :
            if deltaX>0 :
                thetxX=math.pi/2
            elif deltaX<0:
                thetxX = math.pi / 2+math.pi
        else:
            thetxX=math.atan(deltaY/deltaX)

        if int==0 :      #
            self.matriX = self.matrix+self.parsDict['X0']  # 
            self.matriY = self.matriy+self.parsDict['X0']  #
            self.matriX = self.matriX*math.cos(thetxX) -self.matriY*math.sin(thetxX)  #
            self.matriY = self.matriX*math.sin(thetxX) + self.matriY*math.cos(thetxX)  # 
        if int==0 :    #
            self.matriXl = self.matriX-self.parsDict['X0']  # 
            self.matriYl = self.matriY-self.parsDict['X0']  # 
            self.matrixl = self.matriXl * math.cos(thetxX) + self.matriYl * math.sin(thetxX)  # 
            self.matriyl = self.matriYl * math.cos(thetxX) - self.matriXl * math.sin(thetxX)  # 
        return

if __name__ == '__main__':
    myMine = MineParmeters()
    matrix, matriy=myMine.BuldMatrix()
    # plt.figure(figsize=(8, 8))
    # plt.subplot(121)
    # plot1 = plt.imshow(matrix, cmap=plt.cm.jet)
    # plt.colorbar()
    # plt.title('matrix')
    # plt.subplot(122)
    # plot1 = plt.imshow(matriy, cmap=plt.cm.jet)
    # plt.colorbar()
    # plt.title('matriy')
    # plt.show()
    #wx, ix, kx, ux, ex=myMine.caculateSemiInfinite( myMine.parsDict['b'], myMine.parsDict['r'], myMine.matrix)
    #wx, ix, kx, ux, ex = myMine.caculateLimited1( myMine.parsDict['b'],myMine.parsDict['r'], myMine.matrix)
    #wx, ix, kx, ux, ex =myMine.caculateLimited2(myMine.parsDict['b1'],myMine.parsDict['r1'],myMine.parsDict['b2'],
                                                #myMine.parsDict['r2'], matriy)
    #wx, ix, kx, ux, ex =myMine.caculateLimitedAll(matrix, matriy)
    wx, ix, kx, ux, ex =myMine.caculate_any_xy(matrix, matriy,0)
    myMine.planimage(wx, ix, kx, ux, ex)
    myMine.image3D(wx, ix, kx, ux, ex)
    myMine.coord_trans(0)


    plt.figure(figsize=(8, 8))
    plt.subplot(121)
    plot1 = plt.imshow(myMine.matriX, cmap=plt.cm.jet)
    plt.colorbar()
    plt.subplot(122)
    plot1 = plt.imshow(myMine.matriY, cmap=plt.cm.jet)
    plt.colorbar()
    plt.show()


