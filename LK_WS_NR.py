# -*- coding:utf-8 -*-
u"""
        EQUAÇÃO DE ESTADO GENERALIZADA DE LEE-KESLER:
Para a utlização do módulo, instãncias das classes Saturacao (mistura
líquido-vapor) e Prop_Term (liquido comprimido e vapor super-aquecido)
devem ser criadas, respeitando-se as respectivas assinaturas dos seus
construtores, como pode se ver nos exemplos abaixo.
As instâncias da classe Test_Sat podem ser utilizadas para se obter seja a
pressão reduzida de saturação, conhecidos tr e w, ou a temperatura reduzida
de saturação, conhecidos respectivamente pr e w. Novamente, analise com cuidado
os exemplos a seguir. Observe que os resultados estão sendo comparados com
os obtidos através do software fechado CATT2 (Van Wilen & Sonntag (1996)).
20-04-2009
"""

from math import e,log,exp
import math
from robustNR_args import robustNewton
from numpy import array,around


#**********************************************************************************

class Lee_Kesler(object):
    u"""
    Classe com a equacao de estado de Lee and Kesler (Lee-Kesler, 1975)
    em funcao de Tr e vr'. As constantres da equacao de estado sao armazenadas
    no construtor em forma de listas onde o endereço e representado por id,
    assim [id=0] sao as constantes para o fluido simples e [id=1] sao as constantes
    do fluido de referencia (octano). Os coeficientes viriais B, C e D estao
    dispostos em metodos. A equacao de estado completa e calculada no metodo Z
    em funcao de Tr, vr' e id (endereco das listas com as constants do fluido
    """
    def __init__(self):
        
        u""" As constantes abaixo dispostas na forma de lista são para fluido simples
            e de referencia respectivamente usadas no calculo de propriedades
            generalizadas atraves da equacao de estado de Lee and Kesler."""
         
        self.b1   = [0.1181193,    0.2026579]
        self.b2   = [0.265728,    0.331511]
        self.b3   = [0.154790,    0.027655]
        self.b4   = [0.030323,    0.203488]
        self.c1   = [0.0236744,   0.0313385]
        self.c2   = [0.0186984,   0.0503618]
        self.c3   = [0.0,         0.016901]
        self.c4   = [0.042724,    0.041577]
        self.d1   = [0.155488e-4, 0.48736e-4]
        self.d2   = [0.623689e-4, 0.0740336e-4]
        self.beta = [0.65392,     1.226]
        self.gama = [0.060167,    0.03754]
      

    def B(self,Tr,id):
        u"""Coeficiente virial B
        B(self,Tr,id) > B
        Tr = Temperatura reduzida
        id = Endereco da constante.
             Se id = 0 > Fluido Simples
             Se id = 1 > Fluido de Referencia (Octano)"""
        return self.b1[id]-(self.b2[id]/Tr)-(self.b3[id]/pow(Tr,2))-(self.b4[id]/pow(Tr,3))
    
    def C(self,Tr,id):
        u"""Coeficiente virial C
        C(self,Tr,id) > C
        Tr = Temperatura reduzida
        id = Endereco da constante.
             Se id = 0 > Fluido Simples
             Se id = 1 > Fluido de Referencia (Octano)"""
        C=(self.c1[id])-(self.c2[id]/Tr)+(self.c3[id]/pow(Tr,3))
        return C
    
    def D(self,Tr,id):
        u"""Coeficiente virial D
        D(self,Tr,id) > D
        Tr = Temperatura reduzida
        id = Endereco da constante.
             Se id = 0 > Fluido Simples
             Se id = 1 > Fluido de Referencia (Octano)"""
        D=(self.d1[id])+(self.d2[id]/Tr)
        return D
    
    def Z(self, Tr, vr,id):
        u""" funcao que calcula o fator de compressibilidade Z em funcao de Tr,vr,id
        Z(self, Tr, vr, id) > Z
        Tr = Temperatura reduzida
        vr = volume reduzida
        id = Endereco da constante.
             Se id = 0 > Fluido Simples
             Se id = 1 > Fluido de Referencia (Octano)"""
        
        B,C,D=self.B(Tr, id),self.C(Tr, id),self.D(Tr, id)

        Z=1+(B/vr)+(C/pow(vr,2))+(D/pow(vr,5))+(self.c4[id]/(pow(Tr,3)*pow(vr,2)
        ))*(self.beta[id]+(self.gama[id]/pow(vr,2)))*math.exp(-self.gama[id]/pow(vr,2))
        return Z

#*******************************************************************************************

class WuStiel(object):
    """Esta classe será utilizada para a complementação do
procedimento de Lee-Kesler para os fluidos polares, por meio do
método sugerido por Wu-Stiel. Este utiliza a água como referência,
um fator Y de correção e a equação de estado de Keenan para a água.
    """

    def __init__(self): # Constantes para a equação de estado de Keenan
        
        self.A = array((
    [29.492937,     -5.198586,  6.833535,   -0.1564104, -6.397241,  -3.966140,  -0.6904855],
    [-132.13917,    7.777918,   -26.149751, -0.7254611, 26.409282,  15.453061,  2.7407416],
    [274.64632,     -33.301902, 65.326396,  -9.2734289, -47.740374, -29.14247,  -5.1028070],
    [-360.93828,    -16.254622, -26.181978, 4.3125840,  56.323130,  29.568796,  3.9636085],
    [342.18431,     -177.31074, 0.,         0.,         0.,         0.,         0.],
    [-244.50042,    127.48742,  0.,         0.,         0.,         0.,         0.],
    [155.18535,     137.46153,  0.,         0.,         0.,         0.,         0.],
    [5.972849,      155.97836,  0.,         0.,         0.,         0.,         0.],
    [-410.30848,    337.31180,  -137.46618, 6.7874983,  136.87317,  79.847970,  13.041253],
    [-416.05860,    -209.88866, -733.96848, 10.401717,  645.81880,  399.17570,  71.531353]
                        ))

        self.Ta = array((
        [1.544912,    2.5,    2.5,    2.5,    2.5,    2.5,    2.5]
                    ))

        self.Ra = ((
        [0.634, 1., 1., 1., 1., 1., 1.]
                ))

    def __call__(self,Z,Tr,Pr,tol=1e-8):       
        T = 647.29*Tr   
        P = 22.088*Pr # CUIDADO! Em MPa no original...
        t = 1000./T   # 1/K
	
        def difZw(Z):               
            rw = P/(.41615*Z*T)
            Q,DQ,DQT = 0.,0.,0.

            for j in range(7):
                EX = e**(-4.8*rw)*(self.A[8,j] + self.A[9,j]*rw)
                DEX = e**(-4.8*rw)*(-4.8*(self.A[8,j] + 
                            self.A[9,j]*rw) + self.A[9,j])

                QS,DQS = 0.,0.

                for i in range(8):
                    QS += self.A[i,j]*(rw - self.Ra[j])**(i)
                    DQS += (i)*self.A[i,j]*(rw - self.Ra[j])**(i-1)

                Q +=  (t - self.Ta[j])**(j-1)*(QS +   EX)*(t - self.Ta[0])
                DQ += (t - self.Ta[j])**(j-1)*(DQS + DEX)*(t - self.Ta[0])
                DQT += ((t - self.Ta[j])**(j-1)*(QS + EX) + 
                    (j-1)*(t - self.Ta[j])**(j-2)*(QS + EX)*(t - self.Ta[0]) )

            return Z -(1 + rw*Q + rw*rw*DQ),Q,DQT
      
        if Z <= 0.1: #0.2347 'liquido'
            Zin = 0.001
        else:
            Zin = 1.1
            
        Z3 = robustNewton(lambda z,args: difZw(z)[0],Zin)[0]
        self.Zw = Z3
        
        Q,DQT = difZw(Z3)[1:]
        rw = P/(.41615*self.Zw*T)
        self.rw = rw
        self.Hw = -self.Zw + 1. - rw*t*DQT   # (h*-h)/RTc
        self.Sw = -log(self.Zw) + rw*Q - rw*t*DQT  # (s*-s)/R

    
#***********************************************************************************************

class newton_2(object):
    u"""classe que calcula valores vr', dados pr e tr, atraves de um metodo numerico de
    convergencia de equacoes.A funcao Z_ e a funcao de onde se obtem o valor de vr',
    onde o termo -vr'*pr/tr é o lado esquerdo da equacao de estado. O metodo newton_raphson
    é aplicado o metodo de convergencia de vr', que posteriormente deve ser usado na obtencao do
    fator de compressibilidade e outras propriedades termodinamicas a partir da equacao de estado
    de Lee and Kesler. Observa-se que o valor inicial usado no metodo dita quase que completamente
    a fase em que o fluido se encontra.
    """

    def __init__(self,tr,pr,id):
        self.lee_kesler = Lee_Kesler()
        
        """tr = Temperatura reduzida
        vr = volume reduzido
        id = Endereco da constante"""

        self.id=id
        self.pr = pr
        self.tr = tr
  