def rk4_um_passo(f,x,t,h):
    # Implementa o algoritmo Runge-Kutta de 4ta ordem
    # dotx = f(t,x)
    # x0 = numpy.array([x1,...,xn]),
    # t0 : tempo inicial
    # tf : tempo final
    # h : passo de integração
    # as saídas são:
    # t : o vetor tempo 
    # x : o vetor de estados
#N = absolute(floor((tf-t0)/h)).astype(int)
    #x = zeros((N+1,x0.size))
    #t = zeros(N+1)
    #x[0,:] = x0
    #t[0] = t0
#   for i in range(0,N):
    k1 = f(t,x)
    k2 = f(t+h/2,x+(h*k1)/2)
    k3 = f(t+h/2,x+(h*k2)/2)
    k4 = f(t+h,x+h*k3)
    x_out = x+(h/6)*(k1+2*k2+2*k3+k4)
    t_out = t+h
    return t_out, x_out

def rk4u(f,x0,t0,tf,h, gfunc, m):
    # Implementa o algoritmo Runge-Kutta de 4ta ordem
    # dotx = f(t,x,u)
    # x0 = numpy.array([x1,...,xn]),
    # t0 : tempo inicial
    # tf : tempo final
    # h : passo de integração
    # as saídas são:
    # t : o vetor tempo 
    # x : o vetor de estados
    # 
    # gfunc é a função de controle
    # m é o número de entradas de controle
    from numpy import zeros, absolute, floor
    N = absolute(floor((tf-t0)/h)).astype(int)
    x = zeros((N+1, x0.size))
    t = zeros(N+1)
    u = zeros((N+1, m))
    x[0, :] = x0
    t[0] = t0
    for i in range(0, N):
        k1 = f(t[i], x[i,:])
        k2 = f(t[i]+h/2, x[i,:]+(h*k1)/2)
        k3 = f(t[i]+h/2, x[i,:]+(h*k2)/2)
        k4 = f(t[i]+h, x[i,:]+h*k3)
        x[i+1, :] = x[i, :]+(h/6)*(k1+2*k2+2*k3+k4)
        t[i+1] = t[i]+h
        u[i+1,:] = gfunc(t[i], x[i,:])
    return t, x, u

def rk4(f,x0,t0,tf,h):
    # Implementa o algoritmo Runge-Kutta de 4ta ordem
    # dotx = f(t,x)
    # x0 = numpy.array([x1,...,xn]),
    # t0 : tempo inicial
    # tf : tempo final
    # h : passo de integração
    # as saídas são:
    # t : o vetor tempo 
    # x : o vetor de estados
    from numpy import zeros, absolute, floor
    N = absolute(floor((tf-t0)/h)).astype(int)
    x = zeros((N+1, x0.size))
    t = zeros(N+1)
    x[0, :] = x0
    t[0] = t0
    for i in range(0, N):
        k1 = f(t[i], x[i])
        k2 = f(t[i]+h/2, x[i]+(h*k1)/2)
        k3 = f(t[i]+h/2, x[i]+(h*k2)/2)
        k4 = f(t[i]+h, x[i]+h*k3)
        x[i+1, :] = x[i, :]+(h/6)*(k1+2*k2+2*k3+k4)
        t[i+1] = t[i]+h
    return t, x
#
def rkf45a(f,a,b,ya,M,tol):
    # Implementa o algoritmo Runge-Kutta-Fehlberg de 4(5)ta ordem
    # Input
    #   f, a função 
    #   a e b são os pontos a esquerda e direita
    #   ya, é a condição inicial
    #   M é o número de passos
    #   tol é a torlância
    # Saída
    # T, Y, onde T é o vetor das abscissas, Y é o vetor das ordenadas
    #
    from numpy import zeros, absolute, floor, power
    from numpy.linalg import norm
    # Os coeficientes
    a2=1/4
    b2=1/4
    a3=3/8
    b3=3/32
    c3=9/32
    a4=12/13
    b4=1932/2197
    c4=7200/2197
    d4=7296/2197
    a5=1
    b5=439/216
    c5=-8
    d5=3680/513
    e5=-845/4104
    a6=1/2
    b6=-8/27
    c6=2
    d6=-3544/2565
    e6=1859/4104
    f6=-11/40
    r1=1/360
    r3=-128/4275
    r4=-2197/75240
    r5=1/50
    r6=2/55
    n1=25/216
    n3=1408/2565
    n4=2197/4104
    n5=-1/5
    big=1e15
    h=(b-a)/M
    hmin=h/64
    hmax=64*h
    max1=200
    N = absolute(floor((b-a)/hmin)).astype(int)
    Y=zeros((N+1,ya.size))
    Y[0,:]=ya
    T=zeros(N+1)
    T[0]=a
    j=0
    br=b-0.00001*absolute(b)
    while(T[j]<b):
        if((T[j]+h)>br):
            h=b-T[j]
        k1=h*f( T[j], Y[j,:] )
        y2=Y[j,:]+b2*k1
        if(big<norm(y2)):
            break
        k2=h*f(T[j]+a2*h,y2)
        y3=Y[j,:]+b3*k1+c3*k2
        if(big<norm(y3)):
            break
        k3=h*f(T[j]+a3*h,y3)
        y4=Y[j,:]+b4*k1+c4*k2+d4*k3
        if(big<norm(y4)):
            break
        k4=h*f(T[j]+a4*h,y4)
        y5=Y[j,:]+b5*k1+c5*k2+d5*k3+e5*k4
        if(big<norm(y5)):
            break
        k5=h*f(T[j]+a5*h,y5)
        y6=Y[j,:]+b6*k1+c6*k2+d6*k3+e6*k4+f6*k5
        if(big<norm(y6)):
            break
        k6=h*f(T[j]+a6*h,y6)
        err=norm(r1*k1+r3*k3+r4*k4+r5*k5+r6*k6)
        ynew=Y[j,:]+n1*k1+n3*k3+n4*k4+n5*k5
        # Erro a controle do tamanho do passo
        if((err<tol)or(h<2*hmin)):
            Y[j+1,:]=ynew
            if((T[j]+h)>br):
                T[j+1]=b
            else:
                T[j+1]=T[j]+h
            j=j+1
        if(err==0):
            s=0
        else:
            s=power(0.84*(tol*h/err),1/4)
        if((s<0.75)and(h>2*hmin)):
            h=h/2
        if((s>1.50)and(2*h<hmax)):
            h=2*h
        #if((big<norm(Y[j,:]))or(max1==j)):
        #    break
        #M=j
        #if(b>T[j]):
        #    M=j+1
        #else:
        #    M=j
    # truncando os vetores     
    # Posso retornar uma aproximação linear o cúbica nos instantes de tempo solicitados
    Y = Y[:j,:]
    T = T[:j]
    return T,Y
#
def rkf45(f,x0,t0,tf,h):
    # Implementa o algoritmo Runge-Kutta-Fehlberg de 4(5)ta ordem
    # dotx = f(t,x)
    # x0 = numpy.array([x1,...,xn]),
    # t0 : tempo inicial
    # tf : tempo final
    # h : passo de integração
    # as saídas são:
    # t : o vetor tempo 
    # x : o vetor de estados
    from numpy import zeros, absolute, floor, power
    from numpy.linalg import norm
    tol = 2e-5
    h_min = 1e-4
    N = absolute(floor((tf-t0)/h_min)).astype(int)
    x = zeros((N+1,x0.size))
    y = zeros((N+1,x0.size))
    t = zeros(N+1)
    hs = zeros(N+1)
    erro = zeros(N+1)
    x[0,:] = x0
    y[0,:] = x0
    t[0] = t0
    hs[0] = h
    tempo = t0
    i = 0
    print('tempo =', tempo)
    print('tf =', tf)
    while (tempo<=tf):
        k1 = h*f( t[i], x[i] )
        k2 = h*f( t[i] + (1/4)*h, x[i] + (1/4)*k1 )
        k3 = h*f( t[i] + (3/8)*h, x[i] + (3/32)*k1 + (9/32)*k2 )
        k4 = h*f( t[i] + (12/13)*h, x[i] + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3 )
        k5 = h*f( t[i] + h, x[i] + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4)
        k6 = h*f( t[i] + (1/2)*h, x[i] - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5)
        #
        x[i+1,:] = x[i,:] + (25/216)*k1 + (1408/2565)*k3 + (2197/4101)*k4 - (1/5)*k5
        y[i+1,:] = x[i,:] + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6
        #
        #print('erro = ',norm(y[i+1,:]-x[i+1,:]), 'h =',h)
        erro[i] = norm(((1/360)*k1)-((128/4275)*k3)-((2197/75240)*k4)+((1/50)*k5)+((2/55)*k6))
        #print('erro[i]=',erro[i])
        hs[i] = h
        if erro[i] >= tol:
            s = power(tol/(2*norm(y[i+1,:]-x[i+1,:])),1/4)
            print('s=',s)
            h = h*s
        #elif erro[i] < tol:
        #    h = 2*h
        if h<h_min:
            h=h_min            
        t[i+1]=t[i]+h
        tempo = tempo + h
        i=i+1
    # truncando os vetores     
    x = x[:i,:]
    t = t[:i]
    erro = erro[:i]
    hs = hs[:i] 
    return t, x, erro, hs
#
def rkf45b(f,x0,t0,tf,h):
    # Implementa o algoritmo Runge-Kutta-Fehlberg de 4(5)ta ordem
    # dotx = f(t,x)
    # x0 = numpy.array([x1,...,xn]),
    # t0 : tempo inicial
    # tf : tempo final
    # h : passo de integração
    # as saídas são:
    # t : o vetor tempo 
    # x : o vetor de estados
    from numpy import zeros, absolute, floor, power
    from numpy.linalg import norm
    tol = 2e-5
    h_min = 1e-4
    N = absolute(floor((tf-t0)/h_min)).astype(int)
    x = zeros((N+1,x0.size))
    y = zeros((N+1,x0.size))
    t = zeros(N+1)
    hs = zeros(N+1)
    erro = zeros(N+1)
    x[0,:] = x0
    y[0,:] = x0
    if (tf>t0):
        t[0] = t0
    else:
        t[0] = tf
    hs[0] = h
    #tempo = t0
    flag_error = False
    i = 0
    for i in range(N):
        k1 = h*f( t[i], x[i] )
        k2 = h*f( t[i] + (1/4)*h, x[i] + (1/4)*k1 )
        k3 = h*f( t[i] + (3/8)*h, x[i] + (3/32)*k1 + (9/32)*k2 )
        k4 = h*f( t[i] + (12/13)*h, x[i] + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3 )
        k5 = h*f( t[i] + h, x[i] + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4)
        k6 = h*f( t[i] + (1/2)*h, x[i] - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5)
        #
        x[i+1,:] = x[i,:] + (25/216)*k1 + (1408/2565)*k3 + (2197/4101)*k4 - (1/5)*k5
        y[i+1,:] = x[i,:] + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6
        #
        #print('erro = ',norm(y[i+1,:]-x[i+1,:]), 'h =',h)
        erro[i] = norm(((1/360)*k1)-((128/4275)*k3)-((2197/75240)*k4)+((1/50)*k5)+((2/55)*k6))
        #print('erro[i]=',erro[i])
        hs[i] = h
        if erro[i] >= tol:
            s = power(tol/(2*norm(y[i+1,:]-x[i+1,:])),1/4)
            flag_error = True
            h = h*s
        #elif erro[i] < tol:
        #    h = 2*h
        if h<h_min:
            h=h_min            
        t[i+1]=t[i]+h
        #tempo = tempo + h
    # truncando os vetores     
    x = x[:i,:]
    t = t[:i]
    erro = erro[:i]
    hs = hs[:i] 
    if flag_error:
        print('Warning: Erro maior que Tolerância...')
    return t, x, erro, hs
#
def rkf45c(f,x0,t0,tf,h):
    # data: 08/05/2021
    # status: OK
    # permite Forward e Backward integration
    # produz relatório do erro e do tempo utilizado
    # Implementa o algoritmo Runge-Kutta-Fehlberg de 4(5)ta ordem
    # dotx = f(t,x)
    # x0 = numpy.array([x1,...,xn]),
    # t0 : tempo inicial
    # tf : tempo final
    # h : passo de integração
    # as saídas são:
    # t : o vetor tempo 
    # x : o vetor de estados
    from numpy import zeros, absolute, floor, power, mean, std, min, max
    from numpy.linalg import norm
    import time
    tic = time.time()
    N = absolute(floor((tf-t0)/h)).astype(int)
    x = zeros((N+1,x0.size))
    y = zeros((N+1,x0.size))
    t = zeros(N+1)
    erro = zeros(N+1)
    x[0,:] = x0
    y[0,:] = x0
    t[0] = t0
    bck_flag = False
    if (tf>t0):
        print('Forward intergration')
        h = absolute(h)
    else:
        h = -absolute(h)
        print('Backward intergration')
        bck_flag = True
    i = 0
    for i in range(N):
        k1 = h*f( t[i], x[i] )
        k2 = h*f( t[i] + (1/4)*h, x[i] + (1/4)*k1 )
        k3 = h*f( t[i] + (3/8)*h, x[i] + (3/32)*k1 + (9/32)*k2 )
        k4 = h*f( t[i] + (12/13)*h, x[i] + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3 )
        k5 = h*f( t[i] + h, x[i] + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4)
        k6 = h*f( t[i] + (1/2)*h, x[i] - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5)
        #
        x[i+1,:] = x[i,:] + (25/216)*k1 + (1408/2565)*k3 + (2197/4101)*k4 - (1/5)*k5
        y[i+1,:] = x[i,:] + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6
        #
        erro[i] = norm(((1/360)*k1)-((128/4275)*k3)-((2197/75240)*k4)+((1/50)*k5)+((2/55)*k6))
        t[i+1]=t[i]+h
        #print (i/N*100, end="\r") <<<---- REVISAR 08/05/2021
    # truncando os vetores     
    x = x[:i,:]
    t = t[:i]
    erro = erro[:i]
    if (bck_flag):
        t = t[::-1]
        x = x[::-1]
        erro = erro[::-1]
    print('Média do erro   :', "{:.2e}".format(mean(erro)))
    print('Des. padrão     :', "{:.2e}".format(std(erro)))
    print('Erro min.       :', "{:.2e}".format(min(erro)))
    print('Erro max.       :', "{:.2e}".format(max(erro)))
    print('Tempo utilizado :', "{:.2f}".format(time.time() - tic), 'segundos')
    return t, x, erro
#
#
def controlavel(A,B):
    # Verifica a controlabilidade do par (A B)
    # Verifica se dim(A) = n x n, dim(B) = n x m
    # computa a matriz de controlabilidade C = [B A*B ... A^(n-1)*B]
    # Output True se posto(C)=n, False se não
    from numpy.linalg import matrix_power, matrix_rank
    from numpy import nan, concatenate, matmul
    rowsA = A.shape[0]
    colsA = A.shape[1]
    rowsB = B.shape[0]
    colsB = B.shape[1]
    if (rowsA!=colsA) or (colsA!=rowsB):
        return nan
    else:
        n = rowsA
        m = colsB
        C = B
        #print('n =',n)
        #print('C =',C)
        for i in range(1,n):
            C = concatenate((C,matmul(matrix_power(A,i),B)),axis=1)
            # print('i =',i)
            # print('C =',C)
    posto = matrix_rank(C)
    if posto==n:
        return True
    else:
        return False
#    
def P_analitico(tf,t,A,B,F,R,Q):
    # Calcula a solução analítica da Equação de Riccati
    from numpy import matmul, concatenate, diag, sort, argsort, ndarray
    from numpy.linalg import inv, eig
    from scipy.linalg import expm
    # Verifica se R é matriz o escalar
    if (not isinstance(R, ndarray)):
        # E = B*R^(-1)*B^T
        E = matmul(B*(1/R),B.transpose())
    else:
        # E = B*R^(-1)*B^T
        E = matmul(matmul(B, inv(R)),B.transpose())
    # Delta = [[A,-E],[-Q,-A^T]]
    Delta = concatenate((concatenate((A,-E),axis=1),concatenate((-Q,-A.transpose()),axis=1)))
    # calcula os autovalores D, e os autovetores W, não estão ordenados
    D_unsorted, W_unsorted = eig(Delta)
    # ordena os autovalores de menor a maior
    D = sort(D_unsorted)
    # matriz dos autovalores positivos
    M = diag(D[2:4])
    # ordena os autovetores W[:,i] na mesma sequência, correspondente a eig D[i]
    W = W_unsorted[:,argsort(D_unsorted)]
    # W = [[W11,W12],[W21,W22]]
    W11 = W[0:2,0:2]
    W12 = W[0:2,2:]
    W21 = W[2:,0:2]
    W22 = W[2:,2:]
    # T(tf) = -[W22-F*W12]^(-1)*[W21-F*W11]
    Ttf = -matmul(inv(W22-matmul(F,W12)),(W21-matmul(F,W11)))
    # T(t) = e^(-M(tf-t))*T(tf)*e^(-M(tf-t))
    T = matmul(matmul(expm(-M*(tf-t)),Ttf),expm(-M*(tf-t)))
    # P(t) = [W21 + W22*T(t)]*[W11 + W12*T(t)]^(-1)
    return matmul((W21+matmul(W22,T)),inv(W11+matmul(W12,T)))
#
def P_tf_inf(A,B,R,Q):
    # Calcula a solução analítica da equação de Riccati para Tf-> Infty
    from numpy import matmul, concatenate, diag, sort, argsort, ndarray
    from numpy.linalg import inv, eig
    #from scipy.linalg import expm
    # Verifica se R é matriz o escalar
    if (not isinstance(R, ndarray)):
        # E = B*R^(-1)*B^T se R é escalar
        E = matmul(B,B.transpose())*(1/R)
    else:
        # E = B*R^(-1)*B^T
        E = matmul(matmul(B, inv(R)),B.transpose())
    # Delta = [[A ,-E],[-Q,-A^T]]
    Delta = concatenate((concatenate((A,-E),axis=1),concatenate((-Q,-A.transpose()),axis=1)))
    # calcula os autovalores D, e os autovetores W, não estão ordenados
    D_unsorted, W_unsorted = eig(Delta)
    # ordena os autovalores de menor a maior
    D = sort(D_unsorted)
    # ordena os autovetores W[:,i] na mesma sequência, correspondente a eig D[i]
    W = W_unsorted[:,argsort(D_unsorted)]
    # W = [[W11,W12],[W21,W22]]
    W11 = W[0:2,0:2]
    W21 = W[2:,0:2]
    # P(t->Infty) = [W21]*[W11]^(-1)
    return matmul(W21,inv(W11))
#
def SDRE(A,B,R,Q,t0,tf,h,x0):
    # Trajetória sub-ótima SDRE
    # Dada a condição inicial x(t0)=x0
    # * Passo 1: tn=t0, 
    # * Passo 2: Verificar controlabilidade do par (A(x),B)
    # * Passo 3: Calcular a solução P(tn) da EAR com P(tf)=0 e x(tn)
    # * Passo 4: Calcular a entrada de controle sub-ótima u(tn)
    # * Passo 5: Integrar dx(t)/dt em malha fechada um passo de integração e calcular x(tn+h)
    # * Passo 6: Incrementar o tempo num passo de integração tn <- tn+h,
    # * Passo 7: Voltar ao passo 2
    # -----------------------------------------------------------------------------
    from numpy import matmul, floor, zeros
    from numpy.linalg import inv
    #
    def f_malha_fechada(A,B,R,P,t,x):
        return matmul(A-matmul(matmul(matmul(B,inv(R)),B.transpose()),P),x)
    #
    def f1(t,x):
        return f_malha_fechada(A(x),B,R,P,t,x)
    # Passo 1
    N = floor((tf-t0)/h).astype(int)
    x = zeros((N+1,x0.size))
    u = zeros((N+1,B.shape[1]))
    t = zeros(N+1)
    t[0] = t0
    x[0,:] = x0
    i=0
    #
    while (t[i]<=tf):
        #print('i=',i)
        # Passo 2:
        if controlavel(A(x[i,:]),B):
            # Passo 3:
            P = P_tf_inf(A(x[i,:]),B,R,Q)
            # Passo 4:
            u[i,:] = -matmul(matmul(matmul(inv(R),B.transpose()),P),x[i,:])
            # Passo 5: (Runge-Kutta 4)e passo 6
            #print('t[i]=',t[i])
            if (i+1<=N):
                t[i+1],x[i+1,:] = rk4_um_passo(f1,x[i,:],t[i],h)
                # posso utilizar aqui um RKF45 de um passo e verificar a aproximação de 5ta ordem do erro
                i=i+1
            else:
                break
        else:
            print('Par A(x),B não controlável em x=',x[i,:])
            break
    return t,x,u

def lqr_regulador(A, B, Q, R, h, t0, tf, dt, x0, message=False):
    # 13/05/2021
    #
    # LIMITAÇÕES: 
    # O código sõ funciona R escalar, entrada u escalar.
    # para sistema de 2a ordem.
    #
    # Calcula a resposta LQR para minimizar o funcional
    #
    #     J(x(t), p(t), u(t)) = int_0^tf g(x(t), p(t), u(t))dt
    #                         = int_0^tf (x^T*Q*x + u^T*R*u)dt
    #
    # sob a condição:
    #
    #              dx/dt = Ax(t) + Bu(t)
    #
    # com x(0) = x0.
    # Devolve:
    #
    # t, tempo de t0:dt:tf
    # x, as trajetórias ótimas
    # p, os coestados,
    # u, o controle ótimo
    # K, a matriz solução da equação de Riccati, calculada com horizonte infinito
    # g, avaliação do argumento do funcional em cada instante de tempo de t0 a tf
    # J, o valor do funcional
    # H, o valor Hamiltoniano H(x(t), p(t), u(t))
    # -------------------------------------------
    
    import numpy as np
    from scipy.interpolate import interp1d
    from myfunctions import rkf45c, controlavel
    import scipy.integrate as integrate
    
    def myprint(str, message_flag):
        if message_flag: 
            print(str)
            return
        
    # Passo 1: verificar a controlabilidade do par A, B
    myprint('Calculando controlabilidade do par A, B:', message)
    if not(controlavel(A,B)):
        myprint('Par A, B não controlável. Abortando.', message)
        return 
    else:
        myprint('Par A, B controlável.', message)
        
        # Passo 2: calculo da solução da equação diferencial de Ricatti
        myprint('\nCalculando a solução da equação diferenecial de Riccati.', message)
        
        def rhsK(t, K):
            from numpy import array, matmul
            k11 = K[0]
            k12 = K[1]
            k21 = K[2]
            k22 = K[3]
            Ka = array([[k11, k12],[k21, k22]])
            dKdt = -matmul(Ka, A) - matmul(A.transpose(), Ka) - Q + matmul(Ka, matmul(B, matmul(B.transpose(), Ka)))*(1/R)
            return dKdt.reshape(K.shape)
        
        tK, Ksol, erro_K = rkf45c(rhsK, h, tf, t0, dt)
        
        # Essa interpolão pode ser melhorada
        Kt = interp1d(tK, Ksol.transpose(), bounds_error=False, fill_value=(Ksol[0], Ksol[-1]))
        
        K = Ksol[0,:].reshape(A.shape)
        myprint('K =', message)
        myprint(K, message)
        
        # Passo 3: Testar os autovalores do sistema realimentado 
        # Kf = A - B*R^-1*B^T*K
        myprint('\nCalculando os autovalores do sistema realimentado Kf = A - B*R^-1*B^T*K', message)
        Kf = np.real(A - np.matmul(np.matmul(B,B.transpose()),K)*(1/R))
        w, v = np.linalg.eig(Kf)
        myprint('Autovalores de Kf =', message)
        myprint(w, message)
        
        # Passo 4: as trajetórias ótimas
        myprint('\nCalculando as trajetórias ótimas:',message)
        t, x, erro_x = rkf45c(lambda t, x: np.matmul(Kf, x), x0, t0, tf, dt)
        # definindo x1(t) e x2(t)
        xt = interp1d(t, x.transpose(), bounds_error=False, \
                      fill_value=(x[0], x[-1]))
        
        # Passo 5: o controle ótimo u^*(t)
        myprint('\nCalculando o controle ótimo u* = -R^{-1}*B^T*K*x', message)
        u = -(1/R)*np.matmul(np.matmul(B.transpose(),K),x.transpose()).transpose()
        ut = interp1d(t, u.transpose(), bounds_error=False, \
                      fill_value=(u[0], u[-1]))
                        
        # Passo 6: calcular o coestados
        myprint('\nCalculando os coestados p(t)=K(t)x(t)', message)
        pt = lambda t: np.matmul(Kt(t).reshape(A.shape), xt(t))
        
        # Passo 7: calculo de g(x, u)
        def gt(ts):
            xv = xt(ts).reshape(2,1)
            return (1/2)*np.matmul(np.matmul(xv.T, Q), xv).reshape(1) + (1/2)*ut(ts)**2*R

                     
        # Passo 8: calcular o Hamiltoniano    
        myprint('\nCalculando o Hamiltoniano', message)
        
        def Ha(th):
            pv = pt(th).reshape(2,1)
            xv = xt(th).reshape(2,1)
            return gt(th) + np.matmul(pv.T, np.matmul(A, xv)).reshape(1)\
        + np.matmul(pv.T, np.matmul(B, ut(th))).reshape(1)
        
        p = np.zeros((len(t), len(x0)))
        H = np.zeros(len(t))
        g = np.zeros(len(t))
        Kout = np.zeros((len(t), 2*len(x0)))
        
        for i in range(len(t)):
            p[i, :] = pt(t[i])
            g[i] = gt(t[i])
            H[i] = Ha(t[i])
            Kout[i, :] = Kt(t[i])
            
        # Passo 9: Calculando o valor do funcional   
        J, errorJ = integrate.quad(gt, t0, tf) 
        myprint('\nCalculando o valor do funcional', message)
        myprint('J =', message)
        myprint(J, message)
        return t, x, p, u, Kout, g, H, J
    
def lqr_regulador_1(A, B, Q, R, h, t0, tf, dt, x0, message=False):
    # 14/05/2021
    # Nesta versão o u(t) é limitado u_min <= u(t) <= u_max
    # Ainda não foi implementado
    # -----------------------------------------------------------
    # LIMITAÇÕES: 
    # O código sõ funciona R escalar, entrada u escalar.
    # para sistema de 2a ordem.
    #
    # Calcula a resposta LQR para minimizar o funcional
    #
    #     J(x(t), p(t), u(t)) = int_0^tf g(x(t), p(t), u(t))dt
    #                         = int_0^tf (x^T*Q*x + u^T*R*u)dt
    #
    # sob a condição:
    #
    #              dx/dt = Ax(t) + Bu(t)
    #
    # com x(0) = x0.
    # Devolve:
    #
    # t, tempo de t0:dt:tf
    # x, as trajetórias ótimas
    # p, os coestados,
    # u, o controle ótimo
    # K, a matriz solução da equação de Riccati, calculada com horizonte infinito
    # g, avaliação do argumento do funcional em cada instante de tempo de t0 a tf
    # J, o valor do funcional
    # H, o valor Hamiltoniano H(x(t), p(t), u(t))
    # -------------------------------------------
    
    import numpy as np
    from scipy.interpolate import interp1d
    from myfunctions import rkf45c, controlavel
    import scipy.integrate as integrate
    
    def myprint(str, message_flag):
        if message_flag: 
            print(str)
            return
     
    u_min = -1.0
    u_max = 0
    
    # Passo 1: verificar a controlabilidade do par A, B
    myprint('Calculando controlabilidade do par A, B:', message)
    if not(controlavel(A,B)):
        myprint('Par A, B não controlável. Abortando.', message)
        return 
    else:
        myprint('Par A, B controlável.', message)
        
        # Passo 2: calculo da solução da equação diferencial de Ricatti
        myprint('\nCalculando a solução da equação diferenecial de Riccati.', message)
        
        def rhsK(t, K):
            from numpy import array, matmul
            k11 = K[0]
            k12 = K[1]
            k21 = K[2]
            k22 = K[3]
            Ka = array([[k11, k12],[k21, k22]])
            dKdt = -matmul(Ka, A) - matmul(A.transpose(), Ka) - Q + matmul(Ka, matmul(B, matmul(B.transpose(), Ka)))*(1/R)
            return dKdt.reshape(K.shape)
        
        tK, Ksol, erro_K = rkf45c(rhsK, h, tf, t0, dt)
        
        # Essa interpolação pode ser melhorada
        Kt = interp1d(tK, Ksol.transpose(), bounds_error=False, fill_value=(Ksol[0], Ksol[-1]))
        
        K = Ksol[0,:].reshape(A.shape)
        myprint('K =', message)
        myprint(K, message)
        
        # Passo 3: Testar os autovalores do sistema realimentado 
        # Kf = A - B*R^-1*B^T*K
        myprint('\nCalculando os autovalores do sistema realimentado Kf = A - B*R^-1*B^T*K', message)
        Kf = np.real(A - np.matmul(np.matmul(B,B.transpose()),K)*(1/R))
        w, v = np.linalg.eig(Kf)
        myprint('Autovalores de Kf =', message)
        myprint(w, message)
        
        # Passo 4: as trajetórias ótimas
        myprint('\nCalculando as trajetórias ótimas:',message)
        t, x, erro_x = rkf45c(lambda t, x: np.matmul(Kf, x), x0, t0, tf, dt)
        # definindo x1(t) e x2(t)
        xt = interp1d(t, x.transpose(), bounds_error=False, \
                      fill_value=(x[0], x[-1]))
        
        # Passo 5: o controle ótimo u^*(t)
        myprint('\nCalculando o controle ótimo u* = -R^{-1}*B^T*K*x', message)
        u = -(1/R)*np.matmul(np.matmul(B.transpose(),K),x.transpose()).transpose()
        ut = interp1d(t, u.transpose(), bounds_error=False, \
                      fill_value=(u[0], u[-1]))
                        
        # Passo 6: calcular o coestados
        myprint('\nCalculando os coestados p(t)=K(t)x(t)', message)
        pt = lambda t: np.matmul(Kt(t).reshape(A.shape), xt(t))
        
        # Passo 7: calculo de g(x, u)
        def gt(ts):
            xv = xt(ts).reshape(2,1)
            return (1/2)*np.matmul(np.matmul(xv.T, Q), xv).reshape(1) + (1/2)*ut(ts)**2*R

                     
        # Passo 8: calcular o Hamiltoniano    
        myprint('\nCalculando o Hamiltoniano', message)
        
        def Ha(th):
            pv = pt(th).reshape(2,1)
            xv = xt(th).reshape(2,1)
            return gt(th) + np.matmul(pv.T, np.matmul(A, xv)).reshape(1)\
        + np.matmul(pv.T, np.matmul(B, ut(th))).reshape(1)
        
        p = np.zeros((len(t), len(x0)))
        H = np.zeros(len(t))
        g = np.zeros(len(t))
        Kout = np.zeros((len(t), 2*len(x0)))
        
        for i in range(len(t)):
            p[i, :] = pt(t[i])
            g[i] = gt(t[i])
            H[i] = Ha(t[i])
            Kout[i, :] = Kt(t[i])
            
        # Passo 9: Calculando o valor do funcional   
        J, errorJ = integrate.quad(gt, t0, tf) 
        myprint('\nCalculando o valor do funcional', message)
        myprint('J =', message)
        myprint(J, message)
        return t, x, p, u, Kout, g, H, J
