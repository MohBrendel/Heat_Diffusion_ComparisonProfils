from numpy import*
from matplotlib.pyplot import*
from mpl_toolkits import mplot3d
from scipy import integrate
from scipy.integrate import quad
from os import getcwd
from os import mkdir

NFichier = input('Give your repetory name: ')
mkdir(str(NFichier)) #Déclaration d'un dossier
mkdir(str(NFichier)+'/Courbes 2D') #Déclaration d'un dossier
mkdir(str(NFichier)+'/Courbes 3D') #Déclaration d'un dossier
mkdir(str(NFichier)+'/Solutions') #Déclaration d'un dossier

                        #╔════════════════════════════════════════════════════════╗
                        #║               Structures et variables                  ║
                        #╙────────────────────────────────────────────────────────╜

T1=50; TM=400; L=0.1; t_final=360              #Définition de la température T(0)=T(L)=T1 et la température moyenne à T(x,0)
                                               #La longueur de la plaque L, et le temps auquel s'arrête l'algorithme

dt=float(input('Give time step, dt= '))
NbTemps = int(t_final/dt)                       #Nombre d'élément du vecteur temps
t = linspace(0,t_final,NbTemps)                 #Vecteur Temps

D = 4.25*10**(-6)                               #Coefficient de diffusion thermique
gamma = 1/6                                     #gamma = D*dt/dx²
dx = sqrt(D*dt/gamma)                           #Calcul de dx via le calcul de gamma
xmin = 0; xmax = L                              #Bornes du vecteur X
NbPoints = int((xmax-xmin)/dx)                  #Nombre d'élément du vecteur X
x = linspace(xmin,xmax,NbPoints)                #Vecteur X
Const = ones(NbPoints)*50                       #Vecteur de 50

NbCourbe = 15                                   #Déclare le nombre de courbe coloré (à différent instant t)
f = NbTemps/NbCourbe                            #Variable de "Coincidence des courbes" explicité ci-après
                           

ListeDisc=[0]                                   #Liste des discontinuités 
Nom_Fonc=[]                                     #Stockage des noms des fonctions
Tanaly = []                                     #Stockage des temps analytiques
Tnum = []                                       #Stockage des temps analytiques


                        #╔════════════════════════════════════════════════════════╗
                        #║           Structure code série de Fourier              ║
                        #╙────────────────────────────────────────────────────────╜

def Fourier(x,tf,F,T1,a,b):                                                         #Fonction qui prend comme élément, le vecteur x, le temps (à un instant tf) et F la fonction à développer
        T=2*L                                                                   #Période de la fonction
        w=float(2*pi/T)                                                         #Pulsation de la fonction
        N = 50                                                                 #Nombre d'harmoniques
        def f(x):                                                               #Définition de la fonction en évaluant la fonction F
                return eval(F)-T1
        def bn(x,n,w):                                                          #Définition du coefficient de Fourier impaire
                return f(x)*sin(n*x*w)
        BnfImpaire = 0
        for n in range(0,N,1):                                                  #Sommation de N harmoniques
                IntbImpaire = quad(lambda x: bn(x,n,w), a, b)                 #Transformé en série de Fourier de f(x) de 0 à T/2
                BnfImpaire += IntbImpaire[0]*sin(n*x*w)*exp(-D*(n*pi/L)**2*tf)  #Calcul de la diffusion de chaleur par série de Fourier
        def FourierImpaire(x):                                                  #Normalisation de la série de Fourier (diviser par 4/T)
                return (4/T)*BnfImpaire
        R = FourierImpaire(x)                                                   #Ajout de la condition initiale T(0)=T(L)=50°C
        return R

                        #╔════════════════════════════════════════════════════════╗
                        #║        Fonction précalculer à copier coller            ║
                        #╙────────────────────────────────────────────────────────╜

#►[Fonction types exponentielles]◄
a1=36; b1=1778; c1=96.9246486585  #Tentes: a1 + b1*exp(-c1*abs(x-L/2))
a2=50; b2=1396.298781906; c2=5000 #Gaussienne: a2 + b2*exp(-c2*(x-L/2)**2)
#►[Fonction types périodiques]◄
a3=T1; b3=175*pi #sin: a3 + b3*sin(x*pi/L)
a4=T1; b4=700    #sin²:a4 + b4*sin(x*pi/L)**2
#►[Fonction types polynomiales]◄
a5=750; b5=14000                         #valeur absolue: a5 - b5*abs(x-L/2)
a6=575; b6=210000                        #polynôme 2nd degré: a6 - b6*(x-L/2)**2
a7=49.3001; b7=70744.2181642; c7=4*10**7 #Lorentzienne: a7 + b7*1/(c7*(x-L/2)**2+1)
#►[Notre fonction]◄
a8=T1; b8=420000        # a8 + b8*x**2
                        # a8 + b8*(x-L)**2

                        #╔════════════════════════════════════════════════════════╗
                        #║                                                        ║
                        #║                   Coeur du programme                   ║
                        #║                                                        ║
                        #╙────────────────────────────────────────────────────────╜

def Laplacien(N):                                                                       #Définition du laplacien de dimension N
        K = gamma*diag(ones(N-1),1) + gamma*diag(ones(N-1),-1) + (1-2*gamma)*eye(N,N)
        return  K

Nb = int(input('How many Fonctions ? '))
K=zeros(Nb)                                             #Déclaration du stockage des valeurs des temps numériques moyens
FoncM = zeros((NbTemps-1,Nb))                           #Déclaration du stockage des fonctions (la même) à des instants t différents



for ite in range(0,Nb):                                 #Répéter le programme autant de fois qu'il y a de fonction
        F=zeros(NbPoints)                               #Déclaration du stockage de la fonction que l'utilisateur va saisir
        
        Disc = input('Fonction n°'+str(ite+1)+' discontinue ? (y/n): ')         #Demande à l'utilisateur si sa fonction est discontinue

                        #╔════════════════════════════════════════════════════════╗
                        #║                  Fonction discontinue                  ║
                        #╙────────────────────────────────────────────────────────╜

        if (Disc == 'y'):       #Si la fonction est discontinue
                FoncNom=[]                                                                      #Déclaration de la liste des fonctions
                Nom_Fonc.append(input('Give the function name: '))
                NbDisc = int(input('Give the number of discontinuties: '))+2                    #Nombre de discontinuité plus x=0 et x=L
                Division = NbDisc-1                                                              #Division par le nombre de fonction
                Fdisc = zeros((NbPoints,NbDisc-1))                                              #Stockage des fonctions discontinues
                Gdisc = zeros((NbPoints,NbDisc-1))
                ListeDisc.append(eval(input('Give the abscissa n° 1: ')))
                print('for x € [',0,',',ListeDisc[1],']:')
                FoncNom.append(input('Give your function, T0(x)= '))

                Fdisc[:,0]=eval(FoncNom[0])                                                     #Stockage de la première fonction sur la ligne 0
                Fdisc[int(ListeDisc[1]/dx):NbPoints,0]=0                                        #Suppression des valeurs si l'on dépasse la première discontinuité
                F+=Fdisc[:,0]                                                                   #Ajout des valeurs de la fonction discontinue dans l'intervalle de la première discontinuité

                for i in range(1,NbDisc-2):                                                     #S'applique si le nombre de discontinuité est >1

                        ListeDisc.append(eval(input('Give the abscissa n° '+str(i+1)+': ')))
                        print('for x € [',ListeDisc[i],',',ListeDisc[i+1],']:')
                        FoncNom.append(input('Give your function, T0(x)= '))
                        Fdisc[:,i]=eval(FoncNom[i])
                        Fdisc[0:int(ListeDisc[i]/dx),i]=0; Fdisc[int(ListeDisc[i+1]/dx):NbPoints,i]=0
                        F+=Fdisc[:,i]

                ListeDisc.append(L)
                print('for x € [',ListeDisc[NbDisc-2],',',L,']:')
                FoncNom.append(input('Give your function, T0(x)= '))
                Fdisc[:,NbDisc-2]=eval(FoncNom[NbDisc-2])
                Fdisc[0:int(ListeDisc[NbDisc-2]/dx),NbDisc-2]=0
                F+=Fdisc[:,NbDisc-2]

                def Gog1(x,it_t,FoncNom):                                                       #Série de fourier type 1 (intégrale de 0 à L pour toutes les fonctions)
                        G = zeros(NbPoints)                                                     #Déclaration de G
                        Gdisc[:,0]=Fourier(x,it_t,FoncNom[0],T1,0,L)                            #Série de Fourier de la première fonction
                        Gdisc[int(ListeDisc[1]/dx+dx):NbPoints,0]=0                             #Supression de tout ce qui est après la première discontinuité
                        G+=Gdisc[:,0]                                                           #Stockage de la série de Fourier de la première fonction
                        for i in range(1,NbDisc-2):
                                Gdisc[:,i]=Fourier(x,it_t,FoncNom[i],T1,0,L)
                                Gdisc[0:int(ListeDisc[i]/dx+dx),i]=0; Gdisc[int(ListeDisc[i+1]/dx):NbPoints,i]=0
                                G+=Gdisc[:,i]
                        Gdisc[:,NbDisc-2]=Fourier(x,it_t,FoncNom[NbDisc-2],T1,0,L)
                        Gdisc[0:int(ListeDisc[NbDisc-2]/dx),NbDisc-2]=0
                        G+=Gdisc[:,NbDisc-2]
                        for j in range(2,NbPoints-2):
                                G[j]=(G[j-2]+G[j-1]+G[j]+G[j+1]+G[j+2])/5                       #Moyennation de la Séride de Fourier
                        G+=Const                                                                #Ajout de 50 °C
                        return G

                def Gog2(x,it_t,FoncNom):                                                       #Série de fourier type 2 (intégrale de a à b pour toutes les fonctions)
                        G = zeros(NbPoints)
                        Gdisc[:,0]=Fourier(x,it_t,FoncNom[0],T1,0,ListeDisc[1])
                        Gdisc[int(ListeDisc[1]/dx+dx):NbPoints,0]=0
                        G+=Gdisc[:,0]
                        for i in range(1,NbDisc-2):
                                Gdisc[:,i]=Fourier(x,it_t,FoncNom[i],T1,ListeDisc[i],ListeDisc[i+1])
                                Gdisc[0:int(ListeDisc[i]/dx+dx),i]=0; Gdisc[int(ListeDisc[i+1]/dx):NbPoints,i]=0
                                G+=Gdisc[:,i]
                        Gdisc[:,NbDisc-2]=Fourier(x,it_t,FoncNom[NbDisc-2],T1,ListeDisc[NbDisc-2],L)
                        Gdisc[0:int(ListeDisc[NbDisc-2]/dx),NbDisc-2]=0
                        G+=Gdisc[:,NbDisc-2]
                        for j in range(2,NbPoints-2):
                                G[j]=(G[j-2]+G[j-1]+G[j]+G[j+1]+G[j+2])/5
                        G+=Const
                        return G

                aa = 1000; bb = 0; it_t=0                                                       #Déclaration des bornes pour un algorithme de recherche par dichotomie
                for i in range(0,50):                                                           #On répète l'algorithme 50 fois
                        it_t=(aa+bb)/2
                        if (mean((Gog1(x,it_t,FoncNom)+Gog2(x,it_t,FoncNom))/2)-200<=0):          #Si la moyenne est inférieur à 0 le temps prend la valeur de la borne supérieur
                                aa=it_t
                        if (mean((Gog1(x,it_t,FoncNom)+Gog2(x,it_t,FoncNom))/2)-200>0):           #Si la moyenne est supérieur à 0 le temps prend la valeur de la borne inférieur
                                bb=it_t
                        print((mean((Gog1(x,it_t,FoncNom)+Gog2(x,it_t,FoncNom))/2)),'/',200,' for ',Nom_Fonc[ite])
                Tanaly.append(it_t)                                     #On sauvergarde le temps analytique

                        #╔════════════════════════════════════════════════════════╗
                        #║                   Fonction Continue                    ║
                        #╙────────────────────────────────────────────────────────╜

        if (Disc == 'n'):       #Si la fonction est continue
                FoncNom=[]                                                      #Déclaration de la liste des fonctions
                Division=1                                                      #Déclaration variable
                Nom_Fonc.append(input('Give the function name: '))
                FoncNom.append(input('Give your function, T0(x)= '))
                F=eval(FoncNom[0])                                              #On évalue la fonction saisie par l'utilisateur
                G=Fourier(x,0,FoncNom[0],T1,0,L)+T1                             #Transformé en Série de Fourier de la fonction
                aa = 1000; bb = 0; it_t=0                                        
                for i in range(0,50):                                   #On répète l'algorithme 50 fois
                        it_t=(aa+bb)/2
                        if (mean(Fourier(x,it_t,FoncNom[0],T1,0,L)+T1)-200<=0):
                                aa=it_t
                        if (mean(Fourier(x,it_t,FoncNom[0],T1,0,L)+T1)-200>0):
                                bb=it_t
                        print(mean(Fourier(x,it_t,FoncNom[0],T1,0,L)+T1),'/',200,' for',Nom_Fonc[ite])
                Tanaly.append(it_t)

                        #╔════════════════════════════════════════════════════════╗
                        #║                 Résolution numérique                   ║
                        #╙────────────────────────────────────────────────────────╜

        Forme = F                                       #Stockage des valeurs de la fonction de l'utilisateur
        MatTheo = Forme                                 #Initialisation des valeurs
        Lap = Laplacien(NbPoints)                       #Définitions du laplacien

        Fonc = zeros((NbPoints,NbTemps))                #Déclaration du stockage des courbes colorées
        DD = zeros((NbTemps,NbPoints))                  #Déclaration du stockage de MatTheo après l'avoir multiplié par le laplacien
        T12M = TM/2                                     #Déclaration de la moitié de la valuer moyenne
        DimM = [mean(Forme)]                            #Stockage des valeurs moyennes de T(x,0)
        k=0; j=0

        for nt in range(1,NbTemps-1):                   #calcul sur le tout le temps
                MatTheo = matmul(Lap,MatTheo)           #laplacien * fonction (multiplication matricielle)
                MatTheo[0]=T1; MatTheo[NbPoints-1]=T1   #Conditions de Dirichlet
                DD[nt,:] = MatTheo                      #Sauvergarde des valeurs calculées
                DimM.append(mean(DD[nt,:]))             #Calcul de la moyenne de la fonction à un instant 'nt'
                if DimM[nt]>=T12M:                      #Si la valeur moyenne numérique est supérieur ou égale à la valeur TM/2, on sauvegarde les valeur (de sorte à ce que dès que la valeur moyenne numérique est infèrieur à TM/2 on arrête de sauvegarder les valeurs)
                        k = nt                          #Sauvegarde du temps (jusqu'à atteindre t1/2)
                        K[ite]=k                        #Sauvegarde du temps t1/2
                if mod(nt,f) == 0:                                                      #si nt est divisible (et sans reste) par f (ce qui se passe autant de fois qu'il y a de courbe à tracer)
                        Fonc[:,j] = MatTheo                                             #Sauvegarde des courbes
                        print('substep ',(j+1),'/', NbCourbe-1,'on ',(ite+1),'/',Nb)
                        j += 1
        Tnum.append(K[ite])
        FoncM[:,ite] = DimM             #Sauvegarde de la courbe de la diminution de la moyenne de cette fonction

                        #╔════════════════════════════════════════════════════════╗
                        #║                 Tracer des graphiques                  ║
                        #╙────────────────────────────────────────────────────────╜

        colors = cm.rainbow(linspace(0, 1, j+1))                #Création du vecteur couleur
        fig = subplots(figsize=(10, 10))                        #Création d'une figure
        plot(x, Forme,'o',markevery=15,markersize=4, c='black')         #Tracer de la forme initiale (avec les points)
        plot(x, Forme, c='black', label='Initial')
        for m in range(0,j):                                            #Répéter l'opération autant qu'il y a de fonction
                plot(x,Fonc[:,m],'o',markevery=15,markersize=4,color=colors[j-m])                                               #Tracer de la forme (à un instant t)
                plot(x,Fonc[:,m], linewidth=2, color=colors[j-m], label=(f'Time {round((m+1)*t_final/NbCourbe):.0f} sc.'))
                legend(loc = 'upper right')
                xlim(xmin,xmax)
        xlabel('x',fontsize=18)
        ylabel('Temperature (°C)',fontsize=18)
        title('Heat Courbe '+str(Nom_Fonc[ite]),fontsize=18)
        savefig(str(NFichier)+'/Solutions/Heat Courbe'+str(Nom_Fonc[ite])+'.png')                 #Sauvegarde du graphique dans le répertoire

        X, T = meshgrid(x,t)
        fig = subplots(figsize=(10, 10))
        contf = contourf(X[1:len(X)-1], T[1:len(T)-1], DD[1:len(DD)-1], 100, cmap='hot')
        colorbar(contf, label='Temperature (K)')
        xlabel("Distance (x)")
        ylabel("time (t)")
        title('azimutal heat Graphic '+str(Nom_Fonc[ite]),fontsize=18)
        savefig(str(NFichier)+'/Courbes 2D/azimutal heat Graphic'+str(Nom_Fonc[ite])+'.png')

        fig = subplots(figsize=(15, 15))
        ax = axes(projection='3d')
        surf3D = ax.plot_surface(X[1:len(X)-1], T[1:len(T)-1], DD[1:len(DD)-1],cmap='hot', edgecolor='none', alpha=1)
        ax.set_title("Heat diffusion")
        colorbar(surf3D, label='Température (°C)')
        xlabel("Distance (x)")
        ylabel("Time (s)")
        title('3D heat Graphic '+str(Nom_Fonc[ite]),fontsize=18)
        savefig(str(NFichier)+'/Courbes 3D/3D heat Graphic'+str(Nom_Fonc[ite])+'.png')

V123=[]
for j in range(0,len(Tanaly)):
        V123.append(j)

fig = subplots(figsize=(10, 10))
bar(V123,K*dt)
xticks(V123,Nom_Fonc)
ylabel('$Time T1/2M (s)$',fontsize=18)
title('numerical time',fontsize=18)
savefig(str(NFichier)+'/Solutions/numerical time.png')

colors = cm.rainbow(linspace(0, 1, Nb+1))
fig = subplots(figsize=(10, 10))
for m in range(Nb):
        plot(t[1:NbTemps] , FoncM[:,m]/TM, color=colors[Nb-m] , label=('Mean regression', Nom_Fonc[m]))      #Tracer de la diminution de la moyenne des fonction (normalisée)
xlabel('$Time (s)$',fontsize=18)
ylabel('$Temperature (°C)$',fontsize=18)
title('Mean diminution',fontsize=18)
legend()
savefig(str(NFichier)+'/Solutions/Mean diminution.png')


fig = subplots(figsize=(10, 10))
bar(V123,Tanaly)
xticks(V123,Nom_Fonc)
ylabel('$Time T1/2M (s)$',fontsize=18)
title('analytique solution',fontsize=18)
savefig(str(NFichier)+'/Solutions/analytique solution.png')

fich = open(str(NFichier)+'/Solutions/résults.txt', "w")

for i in range(0,Nb):
        fich.write('The analytique time is '+str(Tanaly[i])+' sc and '+str(K[i]*dt)+' sc for numerical time for '+str(Nom_Fonc[i])+'\n') #Sauvegarde des résultats en .txt
fich.close()