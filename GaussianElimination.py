__all__ = ["MatrizEscalonada", "Det", "Resolver"]

import numpy as np

def Pivoteo(M):
    """Reordena la matriz M de modo que abs(M[0,0]) sea el mayor de todas las
componentes abs(M[:,0]). Cuidado: Esta función cambia M. """

    # encuentra la fila que cuyo primer elemento es el mayor
    i = np.argmax( np.abs(M[:,0]) )

    requiere_permutacion =  ( i>0 )  # True=se requiere permutacion.
    
    # intercambia las filas i y 0 (solo si son distintos)
    if requiere_permutacion:
        M[[0,i]] = M[[i,0]] 

    # retorna (True/False, fila permutada)
    # si requiere o no una permutacion.
    # Esto es necesario para calcular determinantes.
    return requiere_permutacion, i   

        
def MatrizEscalonada(M):
    """Reescribe la matriz M como una matriz echelon (escalonada) como
    resultado de eliminacion Gaussiana: 
     
        M[i] = M[i] - M[0] * M[i,0] / M[0,0]

    En el proceso, se realiza pivoteo para que el primer elemento de
    la primera fila sea lo más grande posible. Para ello se requieren
    realizar permutaciones. Esta funcion retorna el número de
    permutaciones requeridas para el pivoteo.
    """

    numero_permutaciones = 0
    MatrizPermutacion = np.zeros( [len(M)]*2, int)
    
    for k in range(len(M)-1):
        m = M[k:, k:]

        # permuta dos filas si el pivote es pequeño;
        # retorna True/False si requiere permutacion
        # y retorna la fila permutada con la primera fila
        requiere_permutacion, fila = Pivoteo(m)
        
        numero_permutaciones += requiere_permutacion
        MatrizPermutacion[k, fila+k] = 1
        MatrizPermutacion[fila+k, k] = 1
        
        # reduce cada fila de la sub-matriz m
        for i in range(1, len(m)):
            m[i] = m[i] - m[0] * m[i,0] / m[0,0]

        # m[1:] = m[1:] - m[0] * m[1:,0,None] / m[0,0]

    return numero_permutaciones, MatrizPermutacion 
    

def Det(M):
    """Calcula la determinante de una matriz M"""
    m = np.copy(M)

    # convierte `m` en una matriz escalonada
    # luego extrae el numero de permutaciones usadas
    numero_permutaciones = MatrizEscalonada(m)[0]

    return (-1)**numero_permutaciones * np.prod( np.diag(m) )


def Resolver(M, b):
    """Encuentra el vector `x` que resuelve `M.x = b`, donde `M` es una
matriz y `b` es un vector o matriz, ambos conocidos."""

    # crea matriz aumentada; si b es vector fila, lo convierte en vector columna 
    if b.ndim==1:
        m = np.hstack( (M, b[:,None]) )
    else:
        m = np.hstack( (M, b) )

    # reescribe `m` como una matriz escalonada
    MatrizEscalonada( m )

    rows, cols = M.shape
    c = m[:,cols:]
    m = m[:,:cols]

    x = np.zeros(b.shape)    

    # resuelve m.x=c donde m es matriz escalonada
    x[-1] = c[-1] / m[-1,-1]
    
    for i in range(rows-2, -1, -1):
        x[i] = ( c[i] - m[i, i+1:] @ x[i+1:] ) / m[i,i]

    return x
    

def Inversa(M):
    """Encuentra la matriz `x` que resuelve `M.x = I`, donde `M` es una
matriz cuadrada e `I` es la matriz identidad. `x` corresponde a la
inversa de `M`."""

    return Resolver(M, np.identity(len(M)) )
