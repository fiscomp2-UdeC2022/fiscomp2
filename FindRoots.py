__all__ = ["Secante", "NewtonRaphson"]


import numpy as np
import matplotlib.pyplot as plt
# import sys, h5py, json


def Secante(f, x0, tol=1e-5, *args, **kwargs):
    """Resuelve f(x)==0 con el metodo de la secante, dados dos valores
    iniciales ingresados en x0=[a,b], hasta llegar a cierta tolerancia
    tol. Si esta rutina se llama con otros argumentos, estos son pasados
    a la funcion f como parametros.

    [argumentos]

    f : Funcion f(x, *args, **kwargs) que depende de un escalar `x` y
        retorna un escalar. *args y **kwargs son posibles parametros
        de la funcion.

    x0: Lista, tupla o array de dos elementos [a,b] que son usados como
        semilla del metodo de la secante.

    tol: Se asume que el metodo converge cuando |f(x)|<tol o |b-a|<tol.

    *args, **kwargs : Otros parametros que son entregados a la funcion `f`.

    [retorna]

    xsol : Numero escalar que es solucion de `f(x)=0`.

    """
    if(len(x0) != 2):
        raise SystemExit(f"argumento `x0` debe un objeto de dos elementos, pero `x0={x0}` fue entregado.")

    a, b = x0

    fa = f(a, *args, **kwargs)
    fb = f(b, *args, **kwargs)

    stop = False
    # itera metodo de secante hasta que |f(b)|<tol o |b-a|<tol.
    while(not stop):
        a, b = b, b - fb * (a-b) / (fa - fb)

        fa, fb = fb, f(b, *args, **kwargs)
        
        stop = np.abs(fb)<tol and np.abs(b-a)<tol

    return b


def NewtonRaphson(f, df, x0, tol=1e-5, *args, **kwargs):
    """Resuelve f(x)==0 con el metodo de Newton-Raphson, usando como
    semilla un valor escalar `x0`, hasta llegar a cierta tolerancia
    tol. Si esta rutina se llama con otros argumentos, estos son
    pasados a la funcion f como parametros.

    [argumentos]

    f : Funcion f(x, *args, **kwargs) que depende de un escalar `x` y
        retorna un escalar. *args y **kwargs son posibles parametros
        de la funcion.

    df : Derivada de la funcion `f`, que acepta los mismos parametros que `f`.

    x0: Numero escalar (real o complejo) que es usado como
        semilla del metodo de Newton-Raphson.

    tol: Se asume que el metodo converge cuando |f(x)|<tol o |f(x)/f'(x)|<tol.

    *args, **kwargs : Otros parametros que son entregados a las funciones `f` y `df`.

    [retorna]

    xsol : Numero escalar (del mismo tipo que x0) que es solucion de `f(x)=0`.

    """
    stop = False

    x = x0
    
    while(not stop):
        val = f(x, *args, **kwargs)
        error = val/df(x, *args, **kwargs)
        x = x - error

        stop = abs(error)<tol and abs(val)<tol
        
    return x


def Biseccion(f, a, b, N=50, tol=1e-5, *args, **kwargs):
    """Retorna todos los valores de `x` que resuelven `f(x)=0` en el
    intervalo a<x<b. Primero busca los ceros usando el método de
    biseccion, dividiendo el intervalo a<x<b en N celdas. En cada una
    de las celdas donde se detecta un cambio de signo de f, se usa el
    metodo de la secante para refinar la busqueda hasta una tolerancia
    tol.

    [argumentos]
    
    f : Funcion f(x, *args, **kwargs) que depende de un escalar `x` y
        retorna un escalar. *args y **kwargs son posibles parametros
        de la funcion.

    a, b : Escalares que representan el intervalo donde buscaremos los ceros de f.

    N : Numero de celdas en que es dividido el intervalo a<b. 

    tol: Se asume que el metodo converge cuando |f(x)|<tol o |f(x)/f'(x)|<tol.

    *args, **kwargs : Otros parametros que son entregados a las funciones `f` y `df`.

    [retorna]

    xsol : Lista de valores que resuelven `f(xsol)=0`. 

    """
    # crea una lista de N valores entre a y b
    xval = np.linspace(a,b,N)

    # evalua f en cada valor de x
    fval = f(xval, *args, **kwargs)

    # busca celdas donde f(x) cambia de signo
    i = np.where(fval[1:] * fval[:-1] < 0)[0]

    zeros = [Secante(f, x0=[a,b], tol=tol, *args, **kwargs) for a,b in zip(xval[i],xval[i+1])]

    return zeros
    
