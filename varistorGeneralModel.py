'''
Varistor - Transmission Line Time Domain

This code uses the Transmission Line Modelling Method to analyse two-wire
transmission lines, any frequency. It includes also source models,
load models, linear and non-linear devices models.

The user input parameters and choose the node to get results in the time domain.

Model developed by Mauro Fazion Filho. Code Copyright: Mauro Fazion Filho.

****************************************************************

Definitions:

K - number os instants - iterations - from 1 to KT
KT - total number of instants - iterations
j - Node number
N - total number of Nodes (transmission line nodes)
RD - Transmission line distributed Resistance
LD - Transmission line distributed Inductance
CD - Transmission line distributed Capacitance
GD - Transmission line distributed Conductance
C - Capacitance by line segment
L - Inductance by line segment
G - Conductance by line segment
R - Resistance by line segment
RL - Load Resistance
RS - Source Resistance
ZL - Load Impedance
Z0 - Characteristic Impedance of line segment
IL[K] - Load Current at instant K
VRL[K] - Voltage reflected for the Load Inductance
VS[K] - Source Voltage at instant K
VL[K] - Voltage at Load Inductance, instant K
VIL[K] - Incident Voltage at Load Inductance, instant K
I[j][K] - Current at node j, instant K
V[j][K] - Voltage at node j, instant K
VD[j][K] - Voltage at node j, Right side, instant K
VE[j][K] - Voltage at node j, Left side, instant K
VID[j][K] - Right side Incident Voltage at node j, instant K
VIE[j][K] - Left side Incident Voltage at node j, instant K
VRD[j][K] - Right side Reflected Voltage from node j, instant K
VRE[j][K] - Left side Reflected Voltage from node j, instant K
W - Source sine Angular frequency
FREQ - Source sine frequency
Deltat - Time step
Deltax - transmission line segment length
COMP - total transmission line length
Tempo - Time

CAP - Capacitance Filter (on the load)
ZC - Filter Impedance

****************************************************************       

'''
import time
from math import *

startTime = time.clock()

def fileOut (V, I, NO, K, Tempo, file):

    inFile = open(file +'.txt', 'r')
    contents = inFile.read()
    cont = contents

    outFile = open(file +'.txt', 'w')
    outFile.write(str(cont)+'\n'+str(Tempo)+','+str(V[NO-1][K])+','+str(V[NO][K])) #writing 'cont' as a string
    outFile.close()


def procuraVAR(ZC, CONST, ALFA,V1,VRES):

    # procuraVAR calculates the value of reflected voltage VRES from the Varistor

    V2 = [0.0]*1500
    ERRO = [0.0]*1500
    i=0
    BAUX = log10(ZC*CONST)

    V2[1]= V1-1.
    OP=1
    ERRO[0]=1.
    SUP=V1
    INF=V2[1]

    if V1==0.0:
        V2[0]=0.0
    elif V1 < 4.0:
        V2[i]=V1
    else:
        while ERRO[i] > 1e-5:
            i=i+1

            A = log10(V1 - V2[i])
            B = BAUX + ALFA*(log10(V1 + V2[i]))
            ERRO[i] = fabs(A-B)
            # print (ERRO[i], V1, V2[i])

            if ((ERRO[i] > ERRO[i-1]) and (i == 2)):
                OP=3

            if OP == 1:
                if fabs(A) < fabs(B):
                    V2[i+1]=V2[i] - 1.
                    OP=1
                else:
                    SUP=V2[i-1]
                    INF=V2[i]
                    V2[i+1]=(SUP+INF)/2
                    OP=2

            elif OP == 2:
                if fabs(A) < fabs(B):
                    SUP=V2[i]
                    INF = INF
                else:
                    SUP = SUP
                    INF = V2[i]
                V2[i+1]=(SUP+INF)/2

            elif OP == 3:
                if fabs(A) < fabs(B):
                    SUP=SUP
                    INF = V2[i]
                else:
                    SUP = V2[i]
                    INF = V2[i-1]
                V2[i+1]=(SUP+INF)/2

            if V2[i] == V2[i-1]:
                ERRO[i] = 1e-6

    VRES = V2[i]
    return VRES



def sourceVoltage(KT, OP, Deltat, VS):

    #  defines the Source Voltage
    
    if OP == 1:
        VSO = float(input('What is the Boost Voltage (V)? '))
        VS[0]=VSO
    elif OP == 2:
        VSO = float(input('What is the Maximum Pulse Voltage (V)? '))
        LARG = float(input('What is the Pulse width (seconds) ?'))
        for K in range(KT):
            ONDA = K*Deltat
            if ONDA < LARG:
                VS[K]=VSO
            else:
                VS[K] = 0
    elif OP ==3:
        VSO = float(input('What is the Maximum Surge Voltage (V)? '))
        for K in range(KT):
            TAUX = K*Deltat
            if TAUX < 1.2e-6:
                VS[K] = (VSO/1.2e-6)*TAUX
            elif TAUX < 100e-6:
                VS[K] = (VSO*((100e-6)-TAUX))/98.8e-6
            else:
                VS[K] = 0
    elif OP ==4:
        VSO = float(input('What is the Maximum Surge Voltage (V)? '))
        for K in range(KT):
            TAUX = K*Deltat
            if TAUX < 250e-6:
                VS[K] = (VSO/250e-6)*TAUX
            elif TAUX < 2500e-6:
                VS[K] = (VSO*((2500e-6)-TAUX))/2500e-6
            else:
                VS[K] = 0
    elif OP ==5:
        VSO = float(input('What is the Maximum Sine Voltage (V)? '))
        FREQ = float(input('What is the Sine Frequency (Hz)? '))
        W = 2*FREQ*pi
        for K in range(KT):
            TAUX = K*Deltat
            VS[K]=VSO*sin(W*TAUX)
    return VS



def main():
    print('')
    print('Inputting data: Transmission line parameters')
    
    KT = int(input('Number of iterations: ',))
    COMP = 400. #float(input('Line length (meters): ', ))
    N = 51 #int(input('Number of Nodes (transmission line nodes): ',))
    RD = 1.0e-5 #float(input('Transmission line distributed Resistance (ohm/m): ',))
    GD = 0.0 #float(input('Transmission line distributed Conductance (siemen/m): ',))
    CD = 1.0e-10 #float(input('Transmission line distributed Capacitance (farad/m): ',))
    LD = 2.5e-7 #float(input('Transmission line distributed Inductance (henry/m): ',))

    print('')
    print('Inputting data: Load and source parameters')

    RL = 50. #float(input('Load Resistance (ohm): ',))
    RS = 0.0 #float(input('Source Resistance (ohm): ',))


    # varistor parameters

    print('')

    print('Inputting data: Varistor parameters (non linear device) ')

    CAP = 5.2e-9 #float(input('Varistor Capacitance (farad): ',))
    RON = 10. # float(input('Varistor Resistance - online - conducting (ohm): ',))
    ROFF = 1.e9 #float(input('Varistor Resistance - offline (ohm): ',))
    
    Vlim1 = 22.0 #float(input('Voltage at I = 1 mA (V): ',))
    Ilim1 = 1.e-6
    Va = 28. #36.0 #float(input('Va to calculate Varistor parameters (V)): ',))
    Ia = 0.01 # 1.0 # float(input('Ia to calculate Varistor parameters (A): ',))
    Vb = 36. #36.0 #float(input('Vb to calculate Varistor parameters (V): ',))
    Ib = 0.5 # 1.0 # float(input('Ib to calculate Varistor parameters (A): ',))
 
    ALFA = (log10(Ib) - log10(Ia))/(log10(Vb) - log10(Va))
    CONST =(Ib/(Vb**ALFA))

    print(ALFA,'   ',CONST)

    #ALFA = 11.55 #float(input('Varistor Alpha: ',))   << S05K11 Siemens >>
    #CONST = 6.58e-18 #float(input('Varistor K constant: ',))
    

    # transmission line parameters
 
    Deltax = COMP/(N-1)
    L = LD*Deltax
    C = CD*Deltax
    R = RD*Deltax
    G = GD*Deltax
    Z0 = sqrt(L/C)
    Deltat = sqrt(L*C)

    ZC = Deltat/(2*CAP)

    # selecting the input voltage signal
    # selecting the node to get results
    
    print('Select the input voltage signal from the following:')
    print('')
    print(' 1 - Boost (voltage of just one time step)')
    print(' 2 - Square pulse')
    print(' 3 - Atmospheric Surge (1.2 x 50 microseconds)')
    print(' 4 - Switching Surge (250 x 2500 microseconds)')
    print(' 5 - Sine wave')
    OP = int(input('What is the option ? 1/2/3/4 or 5 ? '))

    while OP != 1 and OP != 2 and OP != 3 and OP != 4 and OP != 5:
        print ('Incorreto valor de OP')
        OP = int(input('What is the option ? 1/2/3/4 or 5 ?')) 

    VS = [0.]*(KT)

    sourceVoltage(KT, OP, Deltat, VS)
        
    NO = 50 #int(input('What is the node to get results of Voltage and Current? '))

    fileName = input("Enter the name of the results file: ")
    file = str(fileName) #transform the name in string
    outFile = open(file +'.txt', 'w')  # creates the file for results
  
    # Iterative process.
    # Constituted by incidence and reflection (scattering) computation for each time
    # step, and connection to the following time step.
    # References: see Mauro Faccioni Filho Thesis.


    # At first instant, K=0, all currents and incident voltages are equal Zero
    IL = [0.]*(KT)
    VIC = [0.]*(KT)
    VRC = [0.]*(KT)
    VC = [0.]*(KT)
    IVAR = [0.]*(KT)
    
    VID = []
    VIE = []
    I = []
    V = []
    VD = []
    VE = []
    VRD = []
    VRE = []

    for j in range(N):
        line1 = []   # creates a list 
        for K in range (KT):
            line1.append(0.0)  # building the matrix: zero in all nodes
        VID.append(line1)
    for j in range(N):
        line2 = [] 
        for K in range (KT):
            line2.append(0.0)
        VIE.append(line2)
    for j in range(N):
        line3 = []
        for K in range (KT):
            line3.append(0.0)
        I.append(line3)
    for j in range(N):
        line4 = []
        for K in range (KT):
            line4.append(0.0)
        V.append(line4)
    for j in range(N):
        line5 = []
        for K in range (KT):
            line5.append(0.0)
        VD.append(line5)
    for j in range(N):
        line6 = []
        for K in range (KT):
            line6.append(0.0)
        VE.append(line6)
    for j in range(N):
        line7 = []
        for K in range (KT):
            line7.append(0.0)
        VRD.append(line7)
    for j in range(N):
        line8 = []
        for K in range (KT):
            line8.append(0.0)
        VRE.append(line8)
    VRES=0.0


    Tempo = 0.0

    
    #### START TIME ITERATION

    for K in range(KT):

        for j in range(N):

            if j == 0:

                ## excites first node
                if RS == 0.0:
                    V[j][K]= VS[K]
                else:
                    V[j][K] = ((VS[K]/RS) + (2*VID[j][K]/(R+Z0)))/(1/RS + 1/(R+Z0))
                I[j][K] = (V[j][K] - 2*VID[j][K])/(R+Z0)

                VE[j][K] = V[j][K]   #### MAYBE NOT NECESSARY
                VD[j][K] = 2*VID[j][K] + I[j][K]*Z0  

                ## first node replies
                VRD[j][K] = VD[j][K] - VID[j][K]
                VRE[j][K] = VE[j][K] - VIE[j][K]  #### MAYBE NOT NECESSARY

                ## connection to the next timestep
                if (K+1) != KT:
                    VIE[j][K+1] = 0.0
                    VID[j][K+1] = VRE[j+1][K]

            elif j == (N-1):

                ##  region =1, V[j-1][K] <= Vlim1
                ##  region =2, Vlim1 < V[j-1][K] < Vlim2
                ##  region =3, Vlim2 <= V[j-1][K]

                ## excites last node with VARISTOR
                if V[j-1][K] <= Vlim1:
                    V[j][K] = ((2*VIE[j][K])/Z0 + (2*VIC[K])/ZC)/(1/Z0 + 1/RL + 1/ZC)
                    IL[K] = V[j][K]/RL
                    VC[K]=V[j][K]
                elif V[j-1][K] > Vlim1: # and V[j-1][K] < Vlim2: # testing with no RON
                    V[j][K] = ((2*VIE[j][K])/Z0 + (2*VIC[K])/ZC)/(1/Z0 + 1/RL + 1/ZC)
                    IL[K] = V[j][K]/RL
                    # I[j][K] = V[j][K]/RL    # must be the IL
                    VC[K]=V[j][K]
                    IVAR[K]=(V[j][K]-2*VIC[K])/ZC
                else:
                    V[j][K] = (2*VIE[j][K]/Z0)/(1/Z0 + 1/RL + 1/RON)
                    IL[K] = V[j][K]/RL
                    # I[j][K] = V[j][K]/RL    # must be the IL
                    IVAR[K]= V[j][K]/RON
                    VC[K]= 0.0                 

                ## last node replies

                if V[j-1][K] <= Vlim1:
                    VRC[K] = VC[K] - VIC[K]
                    VRE[j][K] = V[j][K] - VIE[j][K]             
                elif V[j-1][K] > Vlim1: # and V[j-1][K] < Vlim2: #testing with no RON
                    VRE[j][K] = V[j][K] - VIE[j][K]
                    VRC[K] = VC[K] - VIC[K]            
                else:
                    VRE[j][K] = V[j][K] - VIE[j][K]
                    VRC[K] = 0.0

                ## last node connects to the next timestep
                if (K+1) != KT:

                    if V[j-1][K] <= Vlim1:
                        VIC[K+1] = VRC[K]
                        VIE[j][K+1] = VRD[j-1][K]
                        
                    elif V[j-1][K] > Vlim1: # and V[j-1][K] < Vlim2: #testing with no RON
                        VIE[j][K+1] = VRD[j-1][K]
                        V1=VRC[K]
                        VIC[K+1]= procuraVAR(ZC, CONST, ALFA,V1,VRES)                                 
                    else:
                        VIE[j][K+1] = VRD[j-1][K]
                        VIC[K+1]= -VRC[K]                        

            else:

                ## excites each node
                V[j][K] = (2*VIE[j][K]/Z0 + 2*VID[j][K]/(R+Z0))/(1/Z0 + 1/(R+Z0) + G)
                I[j][K] = (V[j][K] - (2*VID[j][K]))/(R+Z0)

                VE[j][K] = V[j][K]
                VD[j][K] = 2*VID[j][K] + I[j][K]*Z0

                ## each node replies
                VRD[j][K] = VD[j][K] - VID[j][K]
                VRE[j][K] = VE[j][K] - VIE[j][K]

                ## connection to the next timestep
                if (K+1) != KT:
                    VIE[j][K+1] = VRD[j-1][K]
                    VID[j][K+1] = VRE[j+1][K]


        Tempo= K*Deltat
                      
        ###### # write the results in a txt file
        fileOut (V, I, NO, K, Tempo, file)
        ## testing
        ## print(Tempo,'   ',V[0][K],'  ',V[NO][K])


main()
print (time.clock() - startTime, "seconds")
