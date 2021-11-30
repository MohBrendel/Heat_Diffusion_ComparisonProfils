                        #╔════════════════════════════════════════════════════════╗
                        #║        Fonction précalculer à copier coller T0(x)      ║
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

All the function are starting from 400 to 200 °C.
It's an exemple of Heat diffusion with varying profile and comparing (with different profile) it with the analitycal (with numerical computation yes...).
