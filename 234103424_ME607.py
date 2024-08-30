# Assumption-each lamina is considered to be of the same material (Glass-Epoxy)
# Thickness of each lamina is assumed to be same=0.125 mm (line 30 and 60)
# input Theta for each lamina starting from top(line 61)
# input delta_T and delta_C ( line 58 and 59) and mechanical load and moment (line 54)
# input stiffness values (line 75 to 78) and ultimate strengths (line 55)
# input the values of CTEs in material axes (line 128)
# input the values of CHEs in material axes (line 134)


import numpy as np
import math as mt


def mat_prod(matrix1,matrix2):        #mat_prod function returns the product of 2 matrices
    r1=matrix1.shape[0]
    c1=matrix1.shape[1]
    r2=matrix2.shape[0]
    c2=matrix2.shape[1]
    matrix3=np.zeros([r1,c2])
    for i in range(r1):
            for k in range(c2):
                sum_product=0
                for j in range(r2):
                    sum_product+=matrix1[i][j]*matrix2[j][k]
                matrix3[i][k]=(sum_product)  
    return matrix3


def ABBD_formation(Q_bar,n):          # retruns ABBD matrix based on updated input of Q_bar
    tk=0.125*10**(-3)
    A=np.zeros([3,3])
    B=np.zeros([3,3])
    D=np.zeros([3,3])
    zk_minus_1=-tk*(n//2)   
    for i in range(n):
        A+=Q_bar[i]*tk
        zk=zk_minus_1+tk   
        B+=0.5*Q_bar[i]*(zk**2-zk_minus_1**2)
        D+=(1/3)*Q_bar[i]*(zk**3-zk_minus_1**3)
        zk_minus_1+=tk
    ABBD=np.zeros([6,6])
    for i in range(3):
        for j in range(3):
            ABBD[i][j]=A[i][j]
        for j in range(3,6):
            ABBD[i][j]=B[i][j-3] 
    for i in range(3,6):
        for j in range(3):
            ABBD[i][j]=B[i-3][j]
        for j in range(3,6):
            ABBD[i][j]=D[i-3][j-3] 
    return ABBD

N_M_mech=np.array([100000,0,0,0,0,0]).reshape(6,1)   #input the mechanical force and moment per unit lenth in SI units
su=np.array([1062*10**6,610*10**6,31*10**6,118*10**6,72*10**6]).reshape(5,1)    #input the ultimate strength in the 
                                                                                #order Longitudinal,transverse, shear 
                                                                                #followed sorting by tensile and compressive
delta_T=0                                              # Delta T
delta_c=0
tk=0.125*10**(-3)
Theta=np.array(np.radians([0,45,-45,90,90,-45,45,0]))     
n=len(Theta)
symm=1
for i in range(n//2):                                  # to check if the given laminate is symmetric
    if Theta[i]!=Theta[n-1-i]:
        symm=0
        break
if symm==1:
    print("This is a symmetric laminate with",n,"plies")
else:
    print("This is not a symmetric laminate with",n,"plies")
    
Q=np.zeros([3,3])
#stiffness in longitudinal, trnasverse and shear modulus and poisson's ratio
E1=38.6*10**9
E2=8.27*10**9
v12=0.28
G12=4.14*10**9
v21=E2*v12/E1
#print("v21=",v21)
Q[0][0]=E1/(1-v12*v21)
Q[1][1]=E2/(1-v12*v21)
Q[0][1]=(v12*E2)/(1-v12*v21)
Q[1][0]=(v12*E2)/(1-v12*v21)
Q[2][2]=G12
Q[0][2]=0
Q[2][0]=0
Q[1][2]=0
Q[2][1]=0
#print(Q)
#Theta=np.loadtxt("angles.txt",dtype=float)

Q_bar=np.zeros([n,3,3])
for i in range(n):
    Q_bar[i][0][0]=Q[0][0]*(np.cos(Theta[i]))**4+2*(Q[0][1]+2*Q[2][2])*((np.sin(Theta[i]))*(np.cos(Theta[i])))**2+Q[1][1]*np.sin(Theta[i])**4
    Q_bar[i][1][1]=Q[0][0]*(np.sin(Theta[i]))**4+2*(Q[0][1]+2*Q[2][2])*((np.sin(Theta[i]))*(np.cos(Theta[i])))**2+Q[1][1]*np.cos(Theta[i])**4
    Q_bar[i][0][1]=(Q[0][0]+Q[1][1]-4*Q[2][2])*((np.sin(Theta[i]))*(np.cos(Theta[i])))**2+Q[0][1]*(np.sin(Theta[i])**4+np.cos(Theta[i])**4)
    Q_bar[i][1][0]=(Q[0][0]+Q[1][1]-4*Q[2][2])*((np.sin(Theta[i]))*(np.cos(Theta[i])))**2+Q[0][1]*(np.sin(Theta[i])**4+np.cos(Theta[i])**4)
    Q_bar[i][2][2]=(Q[0][0]+Q[1][1]-2*Q[0][1]-2*Q[2][2])*((np.cos(Theta[i]))*(np.sin(Theta[i])))**2+Q[2][2]*(np.sin(Theta[i])**4+np.cos(Theta[i])**4)
    Q_bar[i][0][2]=(Q[0][0]-Q[0][1]-2*Q[2][2])*(np.sin(Theta[i]))*(np.cos(Theta[i]))**3-(Q[1][1]-Q[0][1]-2*Q[2][2])*(np.cos(Theta[i]))*(np.sin(Theta[i]))**3
    Q_bar[i][2][0]=(Q[0][0]-Q[0][1]-2*Q[2][2])*(np.sin(Theta[i]))*(np.cos(Theta[i]))**3-(Q[1][1]-Q[0][1]-2*Q[2][2])*(np.cos(Theta[i]))*(np.sin(Theta[i]))**3
    Q_bar[i][1][2]=(Q[0][0]-Q[0][1]-2*Q[2][2])*(np.cos(Theta[i]))*(np.sin(Theta[i]))**3-(Q[1][1]-Q[0][1]-2*Q[2][2])*(np.sin(Theta[i]))*(np.cos(Theta[i]))**3
    Q_bar[i][2][1]=(Q[0][0]-Q[0][1]-2*Q[2][2])*(np.cos(Theta[i]))*(np.sin(Theta[i]))**3-(Q[1][1]-Q[0][1]-2*Q[2][2])*(np.sin(Theta[i]))*(np.cos(Theta[i]))**3
#print(Q_bar) 
T=np.zeros([n,3,3])
for i in range(n):
    T[i][0][0]=(np.cos(Theta[i]))**2
    T[i][0][1]=(np.sin(Theta[i]))**2
    T[i][0][2]=2*(np.sin(Theta[i]))*np.cos(Theta[i])
    T[i][1][0]=(np.sin(Theta[i]))**2
    T[i][1][1]=(np.cos(Theta[i]))**2
    T[i][1][2]=-2*(np.sin(Theta[i]))*np.cos(Theta[i])
    T[i][2][0]=-1*(np.sin(Theta[i]))*np.cos(Theta[i])
    T[i][2][1]=(np.sin(Theta[i]))*np.cos(Theta[i])
    T[i][2][2]=(np.cos(Theta[i]))**2-(np.sin(Theta[i]))**2
Active_ply=n
Failed_ply=[]
counter=0
FPF_val=0                                   # to check if first failure has happened
while Active_ply>0:
    counter+=1
    ABBD = ABBD_formation(Q_bar,n)
    #print(ABBD)
    ABBD_inv=np.linalg.inv(ABBD)
    #print(ABBD_inv)
    e_k=np.zeros([6,1])
    #Thermal load
    alpha_12=np.array([8.6*10**(-6),22.1*10**(-6),0/2]).reshape(3,1)   #insert the values of CTEs in material axes
    alpha_xy=np.zeros([n,3,1])
    for i in range(n):
        alpha_xy[i]=mat_prod(np.linalg.inv(T[i]),alpha_12)
        alpha_xy[i][2][0]=alpha_xy[i][2][0]/2
    #Hygroscopic load
    beta_12=np.array([0,0.60,0/2]).reshape(3,1)    # values of beta taken for glass epoxy
    beta_xy=np.zeros([n,3,1])
    for i in range(n):
        beta_xy[i]=mat_prod(np.linalg.inv(T[i]),beta_12)
        beta_xy[i][2][0]=beta_xy[i][2][0]/2
    N_T=np.zeros([3,1])
    N_H=np.zeros([3,1])
    M_T=np.zeros([3,1])
    M_H=np.zeros([3,1])
    zk_minus_1=-tk*(n//2)
    for i in range(n):
        N_T+=mat_prod(Q_bar[i],alpha_xy[i])*tk*delta_T
        N_H+=mat_prod(Q_bar[i],beta_xy[i])*tk*delta_c
        zk=zk_minus_1+tk   
        M_T+=0.5*(zk**2-zk_minus_1**2)*delta_T*mat_prod(Q_bar[i],alpha_xy[i])
        M_H+=0.5*(zk**2-zk_minus_1**2)*delta_c*mat_prod(Q_bar[i],beta_xy[i])
        zk_minus_1+=tk

    N_M=np.zeros([6,1])
    for i in range(3):
        N_M[i][0]=N_T[i][0]+N_H[i][0]+N_M_mech[i][0]
        N_M[i+3][0]=M_T[i][0]+M_H[i][0]+N_M_mech[i+3][0]
    e_k=mat_prod(ABBD_inv,N_M)
    #print(e_k)
    epsilon_xy_0=np.zeros([3,1])
    k=np.zeros([3,1])
    for i in range(3):
        epsilon_xy_0[i][0]=e_k[i][0]
        k[i][0]=e_k[i+3][0]
    #print(epsilon_xy_0)
    zk_minus_1=-tk*(n//2)
    epsilon_xy=np.zeros([n,3,1])
    for i in range(n):
        if zk_minus_1!=0:
            epsilon_xy[i]=epsilon_xy_0+zk_minus_1*k
        else:
            zk_minus_1+=tk
            epsilon_xy[i]=epsilon_xy_0+zk_minus_1*k
            #epsilon_xy[i][2][0]/=2
        zk_minus_1+=tk
    #print(epsilon_xy)
    sigma_xy=np.zeros([n,3,1])
    for i in range(n):
        sigma_xy[i]=mat_prod(Q_bar[i],epsilon_xy[i])

    sigma_12=np.zeros([n,3,1])
     #Free thermal strain in lamina
    free_epsilon_T=np.zeros([n,3,1])
    free_epsilon_C=np.zeros([n,3,1])
    for i in range(n):
        free_epsilon_T[i]=delta_T*alpha_xy[i]
        free_epsilon_C[i]=delta_c*beta_xy[i]

    #Residual Strains calculation
    epsilon_xy_residual=np.zeros([n,3,1])
    for i in range(n):
        epsilon_xy_residual[i] = epsilon_xy[i]-free_epsilon_T[i]-free_epsilon_C[i]
    #print("\nResidual Strains :",epsilon_xy_residual)

    #Residual stress calculation  
    sigma_xy_residual=np.zeros([n,3,1])
    for i in range(n):
        sigma_xy_residual[i]=mat_prod(Q_bar[i],epsilon_xy_residual[i])
    #print("\nResidual Stress xy :",sigma_xy_residual)

    sigma_12_residual=np.zeros([n,3,1])
    for i in range(n):
        sigma_12_residual[i]=mat_prod(T[i],sigma_xy_residual[i])
    #print("\nResidual Stress 12 :",sigma_12_residual)
    sigma_12_net=np.zeros([n,3,1])
    SR=np.zeros([n,4,1])
    Tsai_Hill=np.zeros(n)
    F_index=[]
    Strength_ratio = 0
    Strength_ratio_2=0
    for i in range(n):
        sigma_12[i]=mat_prod(T[i],sigma_xy[i])
        if sigma_12[i][0][0]>=0:
            SR[i][0][0]=sigma_12[i][0][0]/su[0]
        if sigma_12[i][0][0]<0:
            SR[i][0][0]=sigma_12[i][0][0]/su[1]
        if sigma_12[i][1][0]>=0:
            SR[i][1][0]=sigma_12[i][1][0]/su[0]
        if sigma_12[i][1][0]<0:
            SR[i][1][0]=sigma_12[i][1][0]/su[1]
        if sigma_12[i][1][0]>=0:
            SR[i][2][0]=sigma_12[i][1][0]/su[2]
        if sigma_12[i][1][0]<0:
            SR[i][2][0]=sigma_12[i][1][0]/su[3]
        SR[i][3][0]=sigma_12[i][2][0]/su[4]
        Tsai_Hill[i]=SR[i][0][0]**2-SR[i][0][0]*SR[i][1][0]+SR[i][2][0]**2+SR[i][3][0]**2
        print("For the Layer",i+1,"from top for the iteration",counter,"Tsai hill equation evaluates to",Tsai_Hill[i])
        if Tsai_Hill[i]>=1:
            F_index.append(i)
            Failed_ply.append(i+1)
            Active_ply-=1
            Strength_ratio = max(Strength_ratio,max(SR[i]))
        Strength_ratio_2 = max(Strength_ratio_2,max(SR[i]))
        
    #print(sigma_12)
    
    #print(Strength_ratio)
    if Active_ply == n-2 and symm==1 and FPF_val==0:
        FPF = N_M_mech/Strength_ratio
        print("FPF=",FPF)
        FPF_val=1
        
    elif Active_ply==n-1 and symm==0 and FPF_val==0:
        FPF = N_M_mech/Strength_ratio
        print("FPF=",FPF)    
        FPF_val=1
    if Active_ply == 0:
        LPF =N_M_mech/Strength_ratio
        print("LPF",LPF)
    if Strength_ratio==0:
        N_M_mech=N_M_mech/Strength_ratio_2
    else:
        N_M_mech=N_M_mech/Strength_ratio
    #LPF Calculation
    for k in F_index: 
        for i in range(3):
            for j in range(3):
                Q_bar[k][i][j]=0
print("Order of failure of plies in the laminate from top",Failed_ply)



        



    
    


 