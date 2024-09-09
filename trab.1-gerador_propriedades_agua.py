# -*- coding: cp1252 -*-

# Autores: R�binson Erazo e Agmar Pereira
# RAs: 200711491 e 200712331
# Contato: robinson.erazo@hotmail.com ou agmar_filho@hotmail.com

"""Programa para gerar tabelas de propriedades termodin�micas da �gua
   nas regi�es saturadas e fora da satura��o. Foi usado para esse pro-
   grama a equa��o de estado generalizada de Lee-Kesler e rela��es en-
   tre as propriedades termodin�micas."""

from math import *
import LK_WS_NR as LK


"Constantes importantes da �gua e companhia."

Tc = 647.3    # Kelvin
Pc = 22.09 * pow(10,6) #Pa
Ttri = 273.16 # Kelvin
Ptri = 0.6113 * pow(10,3) #Pa
MM = 18.015 # g/mol massa molecular
w = 0.344   # Fator ac�ntrico da �gua
R = 8.3145 #Constante universal dos gases em J/K.mol  



"--------- FUN��ES QUE CALCULAM AS PROPRIEDADES TERMODIN�MICAS --------------"


"Calcula a antidiferencial de Cp(T) para calcularmos integral SCp(T)dT."

def Scp(T): # Cp(T) � dado em J/mol.K
            
    cp0 =  3.224 * pow(10, 1)  # Coeficientes do calor espec�fico para �gua
    cp1 =  1.924 * pow(10,-3)
    cp2 =  1.055 * pow(10,-5) 
    cp3 = -3.596 * pow(10,-9) 

    return cp0*T + cp1*pow(T,2)/2  +  cp2*pow(T,3)/3  +  cp3*pow(T,4)/4    




'---------------------------------------------------------------------------'

"Calcula a antidiferencial de Cp(T)/T para calcularmos SCp(T)/TdT."

def Scp_per_T(T):

    cp0 =  3.224 * pow(10, 1)  # Coeficientes do calor espec�fico
    cp1 =  1.924 * pow(10,-3)
    cp2 =  1.055 * pow(10,-5) 
    cp3 = -3.596 * pow(10,-9) 

    return cp0*log(T) + cp1*pow(T,1)  +  cp2*pow(T,2)/2  +  cp3*pow(T,3)/3   




'----------------------------------------------------------------------------'

"Fun��o que calcula volume espec�fico da �gua"

def v( T , P = 1.0 , Fase_sat = 'nada'): 
 
    Tr = T/Tc # Temperaturas e press�es reduzidas
    Pr = P/Pc

    """O volume espec�fico � determinado pelo fator de compressibilidade
       Z da �gua em um determinado estado. Com Z, usamos ele na equa��o
       de estado do g�s ideal pv = RT ."""

    if Fase_sat == 'nada': 

        Z = LK.H_S(Tr = Tr , Pr = Pr , w = 0.344).prop['Z']
        v =  (Z*R*T) * 1000 / (MM * P) #m3/kg
        return v

    elif Fase_sat == 'liq':

        Z = LK.H_S(Tr = Tr , x = 0.5 , w = 0.344).prop['Zl']
        Psat = LK.H_S(Tr = Tr , x = 0.5 , w = 0.344).prop['Pr'] * Pc #Pa
        
        v = (Z*R*T) * 1000 / (MM * Psat) #m3/kg
        return v

    elif Fase_sat == 'vap':

        Z = LK.H_S(Tr = Tr , x = 0.5 , w = 0.344).prop['Zv']
        Psat = LK.H_S(Tr = Tr , x = 0.5 , w = 0.344).prop['Pr'] * Pc #Pa
        
        v = (Z*R*T) * 1000 / (MM * Psat) #m3/kg
        return v




'----------------------------------------------------------------------------'

"""Fun��o que calcula ENTALPIA da �gua na satura��o, quando � l�quido
   comprimido ou vapor superaquecido. Deve ser especificada a fase no
   caso de regi�o saturada."""

def h( T , P = 1.0 , Fase_sat = 'nada'): 
 
    Tr = T/Tc # Temperaturas e press�es reduzidas
    Pr = P/Pc

    """    Pelo princ�pio de igualdade de varia��o de entalpia de bases
           diferentes, pode se calcular a entalpia da �gua para qualquer
           temperatura e press�o, lembrando que a refer�ncia para essa en-
           talpia � o l�quido saturado no ponto triplo da �gua, arbitra-
           riamente estabelecido como zero. Tamb�m se percebe que � ne-
           cess�rio dimensionalizar os desvios pelo termo R*Tc."""
    
    if Fase_sat == 'nada':  
        
        Dh = LK.H_S(Tr = Tr , Pr = Pr , w = w).prop['h']
        # Desvio de entalpia para temperatura e press�o solicitadas

        Dh_tri = LK.H_S(Tr = Ttri/Tc , Pr = Ptri/Pc , w = w).prop['h']
        # Desvio de entalpia do l�quido saturado no ponto triplo da �gua

        SCpdt = Scp(T)-Scp(Ttri) #Integral

         
        h = (Dh_tri - Dh) * R * Tc + SCpdt #J/mol
        return h / MM # J/g ou kJ/kg


    elif Fase_sat == 'liq':

        Dh = LK.H_S(Tr = Tr , x = 0.5 ,w = w).prop['hl']
        # Desvio de entalpia do l�quido para temperatura e press�o solicitadas

        Dh_tri = LK.H_S(Tr = Ttri/Tc , x = 0.5 , w = w).prop['hl']
        # Desvio de entalpia para l�quido saturado no ponto triplo da �gua

        SCpdt = Scp(T)-Scp(Ttri) #Integral

        h = (Dh_tri - Dh) * R * Tc + SCpdt #J/mol
        return h / MM # J/g ou kJ/kg

    elif Fase_sat == 'vap':

        Dh = LK.H_S(Tr = Tr , x = 0.5 ,w = w).prop['hv']
        Dh_tri = LK.H_S(Tr = Ttri/Tc , x = 0.5 , w = w).prop['hl']
        
        SCpdt = Scp(T)-Scp(Ttri)

        h = (Dh_tri - Dh) * R * Tc + SCpdt #J/mol
        return h / MM #J/g ou kJ/kg



"-----------------------------------------------------------------------------"

"""Fun��o que calcula ENTROPIA da �gua na satura��o, quando � l�quido
   comprimido ou vapor superaquecido. Deve ser especificada a fase no
   caso de regi�o saturada."""

def s( T , P = 1.0 , Fase_sat = 'nada'): 
 
    Tr = T/Tc # Temperaturas e press�es reduzidas
    Pr = P/Pc

    """Pelas princ�pio de igualdade de varia��o de entropia de bases
           diferentes, pode se calcular a entropia da �gua para qualquer
           temperatura e press�o, lembrando que a refer�ncia para a  en-
           tropia � valor no l�quido saturado do ponto triplo da �gua,
           arbitrariamente estabelecido como zero. Tamb�m se percebe que
           � necess�rio dimensionalizar os desvios pelo termo R*Tc.""" 

    if Fase_sat == 'nada':

        Ds = LK.H_S(Tr = Tr , Pr = Pr , w = w).prop['s']
        # Desvio de entropia para temperatura e press�o solicitadas

        Ds_tri = LK.H_S(Tr = Ttri/Tc , Pr = Ptri/Pc , w = w).prop['s']
        # Desvio de entropia para l�quido saturado no ponto triplo da �gua                   

        SCpTdt = Scp_per_T(T) - Scp_per_T(Ttri)
        
        s  = (Ds_tri - Ds) * R + ( SCpTdt ) -  R * log(P/Ptri) #J/K.mol
        return s / MM # J/g.K ou kJ/kg.K


    elif Fase_sat == 'liq':

        Ds = LK.H_S(Tr = Tr , x = 0.5 ,w = w).prop['sl']
        # Desvio de entalpia para temperatura e press�o solicitadas

        Ds_tri = LK.H_S(Tr = Ttri/Tc , x = 0.5 , w = w).prop['sl']
        # Desvio de entalpia para l�quido saturado no ponto triplo da �gua

        Psat = LK.H_S(Tr = Tr , x = 0.5 , w = 0.344).prop['Pr'] * Pc
        #Press�o de satura��o para tal temperatura em Pa

        SCpTdt = Scp_per_T(T) - Scp_per_T(Ttri)

        s  = (Ds_tri - Ds) * R + ( SCpTdt ) -  R * log(Psat/Ptri) #J/K.mol
        return s / MM # J/g.K ou kJ/kg.K

    elif Fase_sat == 'vap':

        Ds = LK.H_S(Tr = Tr , x = 0.5 ,w = w).prop['sv']
        Ds_tri = LK.H_S(Tr = Ttri/Tc , x = 0.5 , w = w).prop['sl']
        SCpTdt = Scp_per_T(T) - Scp_per_T(Ttri)
        Psat = LK.H_S(Tr = Tr , x = 0.5 , w = 0.344).prop['Pr'] * Pc 

        s  = (Ds_tri - Ds) * R + ( SCpTdt ) -  R * log(Psat/Ptri) #J/K.mol
        return s / MM #J/g.K ou kJ/kg.K



"------------------------------------------------------------------------------"

"""Fun��o que calcula ENERGIA INTERNA da �gua na satura��o, quando � l�quido
   comprimido ou vapor superaquecido. Deve ser especificada a fase no
   caso de regi�o saturada."""

def u( T , P = 1.0 , Fase_sat = 'nada'): 
 
    Tr = T/Tc # Temperaturas e press�es reduzidas
    Pr = P/Pc

    """Para calcular energia interna usamos a defini��o de entalpia
       H = U + PV, no caso, usando como refer�ncia o l�quido saturado
       no ponto triplo da �gua, onde a s, h e u s�o zero. Assim, a e-
       nergia interna espec�fica � calculada por
       u(T,P)=h(T,P)-(P*v(T,P)-Ptri*vtri)."""

    if Fase_sat == 'nada': 

        h_ = h(T = T , P = P)  # J/g
        v_ = v(T = T , P = P)  # m3/kg
        vtri =  v(T = Ttri , Fase_sat = 'liq') #Vol. espec�fico do l�q. m3/kg
                                       
        u = 1000 * h_ - (P * v_ - Ptri*vtri) #J/kg
        return u / 1000 #kJ/kg


    elif Fase_sat == 'liq':

        h_ = h(T = T , Fase_sat = 'liq')  # J/kg
        v_ = v(T = T , Fase_sat = 'liq')  # m3/kg
        vtri =  v(T = Ttri , Fase_sat = 'liq') #Vol. espec�fico do l�q. m3/kg
        Psat = LK.H_S(Tr = Tr , x = 0.5 , w = 0.344).prop['Pr'] * Pc #Pa

        u = 1000 * h_ - (Psat * v_ - Ptri*vtri) #J/kg
        return u / 1000 #kJ/kg

    elif Fase_sat == 'vap':

        h_ = h(T = T , Fase_sat = 'vap')  # J/kg
        v_ = v(T = T , Fase_sat = 'vap')  # m3/kg
        vtri =  v(T = Ttri , Fase_sat = 'liq') #Vol. espec�fico do l�q. m3/kg
        Psat = LK.H_S(Tr = Tr , x = 0.5 , w = 0.344).prop['Pr'] * Pc #Pa

        u = 1000 * h_ - (Psat * v_ - Ptri*vtri) #J/kg
        return u / 1000 #kJ/kg



"----------------------------------TABELAS-------------------------------------------------------------------------------------------------------------------------"


"""Essa fun��o gera a tabela das propriedades termodin�micas
   da �gua saturada em fun��o da temperatura"""

def tab_ag_sat_temp():
    
    print '\t\t\t\t\t\t Tabela 1 - Agua Saturada em funcao da temperatura. \n'
    print 'Temp.     Pressao|   Volume Especifico (m3/kg)    |      Energia Interna (kJ/kg)      |         Entalpia (kJ/kg)          |        Entropia (kJ/kg.K)' 
    print '(Celsius)  (kPa) | Liquido sat.        Vapor sat. | Liquido sat.   Evap.   Vapor sat. | Liquido sat.   Evap.   Vapor sat. | Liquido sat.   Evap.   Vapor sat'
    print '-' * 155
    
    #Constru��o das temperaturas de entrada

    lista_T = range(5,371,5) #Lista homog�nea range(5,370,5)

    #Inserindo valores n�o dados pelo range()

    lista_T.insert(0,0.01)   #Temperatura do ponto triplo
    lista_T.append(374.13)   #Temperatura quase cr�tica para n�o
                             #dar indefini��o no m�dulo LK


    for T in lista_T : #Criando a tabela . . . ufa!

        

        # Propriedades do vapor saturado
        Vv = v(T = T+273.15, Fase_sat='vap')
        Uv = u(T = T+273.15, Fase_sat='vap')
        Hv = h(T = T+273.15, Fase_sat='vap')
        Sv = s(T = T+273.15, Fase_sat='vap')

        # Propriedades do l�quido saturado
        Vl = v(T = T+273.15, Fase_sat='liq')
        Ul = u(T = T+273.15, Fase_sat='liq')
        Hl = h(T = T+273.15, Fase_sat='liq')
        Sl = s(T = T+273.15, Fase_sat='liq')
        
        # Propriedades para evapora��o
        Ulv = Uv - Ul
        Hlv = Hv - Hl
        Slv = Sv - Sl

        if T == 0.01:

            #Press�o de satura��o em kPa
            Psat = LK.H_S(Tr = (T+273.15)/Tc, w=w, x=.5 ).prop['Pr']*Pc / 1000

            print '%.2f\t%.5f   \t%.6f \t %.3f    \t %.2f   %.1f   %.1f   \t%.2f       %.1f    %.1f    \t%.4f    %.4f    %.4f' %(T,Psat,Vl,Vv,Ul,Ulv,Uv,Hl,Hlv,Hv,Sl,Slv,Sv)

        if T == 374.13:

            #Press�o de satura��o em MPa
            Psat = LK.H_S(Tr = (T+273.15)/Tc, w=w, x=.5 ).prop['Pr']*Pc / 1000000

            print '%.2f\t%.5f   \t%.6f \t %.3f    \t %.2f   %.1f   %.1f   \t%.2f       %.1f    %.1f    \t%.4f    %.4f    %.4f' %(T,Psat,Vl,Vv,Ul,Ulv,Uv,Hl,Hlv,Hv,Sl,Slv,Sv)

        if T < 100 and T != 0.01:
 
            #Press�o de satura��o em kPa
            Psat = LK.H_S(Tr = (T+273.15)/Tc, w=w, x=.5 ).prop['Pr']*Pc / 1000

            print '%.0f\t%.5f   \t%.6f \t %.3f    \t %.2f   %.1f   %.1f   \t%.2f       %.1f    %.1f    \t%.4f    %.4f    %.4f' %(T,Psat,Vl,Vv,Ul,Ulv,Uv,Hl,Hlv,Hv,Sl,Slv,Sv)

        if T >= 100 and T != 374.13:

            #Press�o de satura��o em MPa
            Psat = LK.H_S(Tr = (T+273.15)/Tc, w=w, x=.5 ).prop['Pr']*Pc / 1000000

            if T == 100 or T == 200:  

                print '\n' 
                print 'Temp.     Pressao|   Volume Especifico (m3/kg)    |      Energia Interna (kJ/kg)      |         Entalpia (kJ/kg)          |        Entropia (kJ/kg.K)' 
                print '(Celsius)  (MPa) | Liquido sat.        Vapor sat. | Liquido sat.   Evap.   Vapor sat. | Liquido sat.   Evap.   Vapor sat. | Liquido sat.   Evap.   Vapor sat'
                print '-' * 155

            print '%.0f\t%.5f   \t%.6f \t %.3f    \t %.2f   %.1f   %.1f   \t%.2f       %.1f    %.1f    \t%.4f    %.4f    %.4f' %(T,Psat,Vl,Vv,Ul,Ulv,Uv,Hl,Hlv,Hv,Sl,Slv,Sv)    

                                       




"--------------------------------------------------------------------------------------------------------------------------------------------------------------"

"""Essa fun��o gera a tabela das propriedades termodin�micas
   da �gua saturada em fun��o da press�o"""

def tab_ag_sat_press():

    print '\t\t\t\t\t\t Tabela 2 - Agua Saturada em fun��o da pressao. \n'
    print 'Pressao   Temp.  |   Volume Especifico (m3/kg)    |      Energia Interna (kJ/kg)      |         Entalpia (kJ/kg)          |        Entropia (kJ/kg.K)' 
    print '(kPa)   (celsius)| Liquido sat.        Vapor sat. | Liquido sat.   Evap.   Vapor sat. | Liquido sat.   Evap.   Vapor sat. | Liquido sat.   Evap.   Vapor sat'
    print '-' * 155


    #Constru��o das press�es de entrada

    lista_P = [(0.5 + x * 0.5)/1000 for x in range(1,15)] #Tudo em MPa
    lista_P = lista_P + [(5. + x * 5)/1000 for x in range(1,15)] # At� 75 kPa
    lista_P = lista_P + [0.075 + x * 0.025 for x in range(1,13)]#At� 0.375 MPa
    lista_P = lista_P + [0.35 + x * 0.05 for x in range(1,14)]#De 0.40 a 1.0 MPa
    lista_P = lista_P + [1.10,1.20,1.30,1.4,1.50,1.75,2.0,2.25,2.5,3.0,3.5]
    lista_P = lista_P + [ x for x in range(4,23)] #De 4 a 22 MPa
    
    #Inserindo valores n�o dados pelos range()
    lista_P.insert(0,0.0006113) #Press�o em MPa no ponto triplo
    lista_P.append(22.08)    #Press�o quase cr�tica para n�o
                             #dar indefini��o no m�dulo LK


    for P in lista_P : #Criando a segunda tabela . . . yes!

        #Temperatura de satura��o Tsat = T em K
        T = LK.H_S(Pr = (1000000.0*P)/Pc, w=w, x=.5 ).prop['Tr']*Tc #T = Tsat em K

        # Propriedades do vapor saturado
        Vv = v(T = T , Fase_sat='vap')
        Uv = u(T = T , Fase_sat='vap')
        Hv = h(T = T , Fase_sat='vap')
        Sv = s(T = T , Fase_sat='vap')

        # Propriedades do l�quido saturado
        Vl = v(T = T , Fase_sat='liq')
        Ul = u(T = T , Fase_sat='liq')
        Hl = h(T = T , Fase_sat='liq')
        Sl = s(T = T , Fase_sat='liq')
        
        # Propriedades para evapora��o
        Ulv = Uv - Ul
        Hlv = Hv - Hl
        Slv = Sv - Sl

        if P == lista_P[0]:

            #Temperatura de satura��o em K
            T = LK.H_S(Pr = (1000000.0 * P)/Pc, w=w, x=.5 ).prop['Tr']*Tc 
            P_ = P*1000 #Press�o em kPa

            print '%.4f\t%.2f    \t%.6f \t %.3f    \t %.2f   %.1f   %.1f   \t%.2f       %.1f    %.1f    \t%.4f    %.4f    %.4f' %(P_,T-273.15,Vl,Vv,Ul,Ulv,Uv,Hl,Hlv,Hv,Sl,Slv,Sv)

        if P == 22.08:

            #Temperatura de satura��o em K
            T = LK.H_S(Pr = (1000000.0*P)/Pc, w=w, x=.5 ).prop['Tr']*Tc

            print '%.2f\t%.2f   \t%.6f \t %.6f    \t %.2f    %.1f   %.1f   \t%.2f      %.1f      %.1f    \t%.4f    %.4f    %.4f' %(P,T-273.15,Vl,Vv,Ul,Ulv,Uv,Hl,Hlv,Hv,Sl,Slv,Sv)

        if P < 0.100 and P != lista_P[0]:
 
            #Temperatura de satura��o em K
            T = LK.H_S(Pr = (1000000.0*P)/Pc, w=w, x=.5 ).prop['Tr']*Tc
            P_= P*1000 #Press�o em kPa

            print '%.1f\t%.2f   \t%.6f \t %.3f    \t %.2f   %.1f   %.1f   \t%.2f       %.1f    %.1f    \t%.4f    %.4f    %.4f' %(P_,T-273.15,Vl,Vv,Ul,Ulv,Uv,Hl,Hlv,Hv,Sl,Slv,Sv)

        if P >= 0.100 and P != 22.08:

            #Temperatura de satura��o em K
            T = LK.H_S(Pr = (1000000.0*P)/Pc, w=w, x=.5 ).prop['Tr']*Tc

            if P == 0.100 or P == 1.4:  

                print '\n' 
                print 'Pressao   Temp.  |   Volume Especifico (m3/kg)    |      Energia Interna (kJ/kg)      |         Entalpia (kJ/kg)          |        Entropia (kJ/kg.K)' 
                print '(MPa)   (celsius)| Liquido sat.        Vapor sat. | Liquido sat.   Evap.   Vapor sat. | Liquido sat.   Evap.   Vapor sat. | Liquido sat.   Evap.   Vapor sat'
                print '-' * 155

            print '%.3f\t%.2f   \t%.6f \t %.6f    \t %.2f   %.1f   %.1f   \t%.2f       %.1f    %.1f    \t%.4f    %.4f    %.4f' %(P,T-273.15,Vl,Vv,Ul,Ulv,Uv,Hl,Hlv,Hv,Sl,Slv,Sv)                                                       




"--------------------------------------------------------------------------------------------------------------------------------------------------------"

"""Essa fun��o gera a tabela das propriedades termodin�micas
   do vapor d'�gua superaquecido"""

def tab_vapor():

    #Constru��o das press�es de entrada

    lista_P = [0.01,0.05] #Tudo em MPa
    lista_P = lista_P + [x * 0.10 for x in range(1,21)] # At� 2.00 MPa
    lista_P = lista_P + [2.50,3.00,3.50,4.0,4.5,5.0,6.0,7.0,8.0,9.0,10.0]#MPa
    lista_P = lista_P + [12.5,15.0,17.5,20.0,25.0,30.0,35.0,40.0,50.0,60.0]#MPa

    print '\t\t\t\t\t\t Tabela 3 - Vapor de agua superaquecido. \n'

    for p in lista_P:
        
        #Temperatura de satura��o em K
        Tsat = LK.H_S(Pr = (1000000.0*p)/Pc, w=w, x=.5 ).prop['Tr']*Tc

        
        print '<<PRESSAO>> = %.2f MPa || Temperatura(Celsius) | Volume Especifico (m3/kg) |  Energia Interna (kJ/kg) |  Entalpia (kJ/kg) | Entropia (kJ/kg.K)' %p
        print '-' * 141

        for T in [Tsat-273.15]+range(375,1301,25):

            #Propriedades do vapor superaquecido
            S = s(T = T+273.15, P = p * 1000000.) #MPa
            H = h(T = T+273.15, P = p * 1000000.)
            V = v(T = T+273.15, P = p * 1000000.)
            U = u(T = T+273.15, P = p * 1000000.)

            print '\t\t\t\t%.0f \t\t\t%.6f \t\t\t%.1f \t\t\t%.1f \t\t%.4f' %(T, V, U, H, S)

    












"--------------------------------------------------------------------------------------------------------------------------------------------------------"
    
"""Essa fun��o gera a tabela das propriedades termodin�micas
   da �gua l�quida comprimida"""

def tab_liq_compr():


    #Constru��o das press�es de entrada

    lista_P = [5.0,10,15,20,30] #Tudo em MPa

    print '\t\t\t\t\t\t Tabela 4 - Liquido comprimido. \n'

    for p in lista_P:
        
        #Temperatura de satura��o em K
        Tsat = LK.H_S(Pr = (1000000.0*p)/Pc, w=w, x=.5 ).prop['Tr']*Tc

        
        print '<<PRESSAO>> = %.2f MPa || Temperatura(Celsius) | Volume Especifico (m3/kg) |  Energia Interna (kJ/kg) |  Entalpia (kJ/kg) | Entropia (kJ/kg.K)' %p
        print '-' * 141

        for T in [Tsat-273.15]+range(0,381,20):

            #Propriedades do vapor superaquecido
            S = s(T = T+273.15, P = p * 1000000.) #MPa
            H = h(T = T+273.15, P = p * 1000000.)
            V = v(T = T+273.15, P = p * 1000000.)
            U = u(T = T+273.15, P = p * 1000000.)

            print '\t\t\t\t%.0f \t\t\t%.6f \t\t\t%.1f \t\t\t%.1f \t\t%.4f' %(T, V, U, H, S)






'---------------------------------- PRINCIPAL ---------------------------------'

if __name__ == '__main__':

    print """\t***Programa para gerar tabelas de propriedades termodin�micas da �gua
             nas regi�es saturadas e fora da satura��o. Foi usado para esse pro-
             grama a equa��o de estado generalizada de Lee-Kesler e rela��es en-
             tre as propriedades termodin�micas."""

    print """\n\nDigite o numero da tabela que voce quer gerar?

                    1) Agua saturada - temperatura
                    2) Agua saturada - pressao
                    3) Vapor Super-aquecido
                    4) Liquido Comprimido

                Obs: a tabela 3 demora uns 10min. pra ser feita!
            """
    
    opc = raw_input()

    if opc == '1':
        tab_ag_sat_temp()

    if opc == '2':
        tab_ag_sat_press()

    if opc == '3':
        tab_vapor()

    if opc == '4':
        tab_liq_compr()






        
        
