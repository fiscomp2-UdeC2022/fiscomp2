# Módulos del curso Física Computacional II (UdeC)

En este repositorio, encontrará varios módulos escritos para el curso de Física Computacional II.

- [FindRoots](FindRoots.py): Rutinas para encontrar ceros de funciones de una variable.
  - `Secante(f, x0)` encuentra una solución de `f(x)=0` cerca de `x=x0` usando el método de la secante. 
  - `NewtonRaphson(f, df, x0)` usa el método de Newton-Raphson para encontrar una solución de `f(x)=0` a partir de `x=x0`. Es necesario conocer la derivada analítica `df(x)` de la función `f(x)`.
  - `Biseccion(f, a, b)` encuentra todos los ceros de `f(x)` en el intervalo `a<x<b` usando el método de la bisección. Las soluciones son refinadas con el método de la secante.
  
- [GaussianElimination](GaussianElimination.py): Implementa el método de eliminación Gaussiana con pivoteo para algunas operaciones matriciales elementales.
  - `Det(M)` calcula la determinante de una matriz `M`.
  - `Inversa(M)` encuentra la matriz inversa `M^{-1}`.
  - `Resolver(M,b)` encuentra el vector `x` que resuelve `M.x=b`, con `M` una matriz cuadrada y `b` un vector conocido.

- [PolyClass](PolyClass.py): Define la clase `Poly`, el cual corresponde a polinomios de variable real o compleja. La clase define algunas operaciones elementales entre polinomios, como la suma, resta, multiplicación y derivada de polinomios.

- [PolyInterpolation](PolyInterpolation.py): Define tres métodos para calcular el polinomio interpolante de una serie de punto `yi=f(xi)`. Estos tres métodos son programados de modo que heredan las propiedades de la clase `Poly` (ver [PolyClass](PolyClass.py). 
  - `Vandermonde(xi, yi)` encuentra los coeficientes del polinomio más
    simple que interpola los puntos `yi=f(xi)` usando inversión de la
    matriz de Vandermonde generada a partir de los puntos `xi`.
  - `Lagrange(xi, yi)` construye polinomios de Lagrange usando la definición de forma directa.
  - `Neville(xi, yi)` construye polinomios de Lagrange usando el método iterativo de Neville.
  
- [Quadratures](Quadratures.py): Define algunos métodos de integración de una función de una variable.
  - `Rectangle(f, a, b, N)` usa el método del rectángulo para estimar la integral de `f(x)` entre `a<x<b` usando `N` nodos.
  - `MidPoint(f, a, b, N)` usa el método del punto medio.
  - `Trapezoide(f, a, b, N)` usa el método del trapezoide.
  - `Romberg(f, a, b, N)` usa el método de Romberg, el cual es un
    método iterativo basado en la regla trapezoidal y en extrapolación
    de Richardson.


## Forma de uso
Para descargar este repositorio, simplemente escriba:
```git 
git clone git@github.com:fiscomp2-UdeC2022/fiscomp2.git
```
	
Luego, debe indicarle a python cómo encontrar estos módulos. Si estos se descargaron en la carpeta `/home/user/fiscomp2`, entonces:
- En **linux**, debe agregar la siguiente línea al archivo `~/.bashrc`:
  ```bash
  export PYTHONPATH=$PYTHONPATH:/home/user/fiscomp2
  ```
  
- En **Mac**, debe agregar la siguiente línea al archivo `~/.bash_profile`:
  ```bash
  export PYTHONPATH="${PYTHONPATH}:/home/user/fiscomp2"
  ```
	
- En **Windows**, en el menú busque `Computador` y use el botón secudario de su mouse en él. 
  - Seleccione `Propiedades` y luego busque la pestaña `Configuraciones avanzadas de sistema`. 
  - Luego, busque el botón `variables de entorno`. 
  - Finalmente, edite (o cree si es encesario) la variable `PYTHONPATH` de modo que incluya `/home/user/fiscomp2;` 
  
  No tengo windows para probar estos pasos :yum: Si hay algún paso que falte o es impreciso, ¡agradezco sus correcciones! Me basé principalmente en instrucciones de [aquí](https://stackoverflow.com/questions/3701646/how-to-add-to-the-pythonpath-in-windows-so-it-finds-my-modules-packages) y [aquí](https://bic-berkeley.github.io/psych-214-fall-2016/using_pythonpath.html). 


Ahora debería poder importar los módulos normalmente. Por ejemplo, para usar la función `Quadratures.Romberg`, puede escribir:
```python
from Quadratures import Romberg
f = lambda x : x**2

Romberg(f, a=0, b=5, show=True)
```
lo cual debería retornar:

```python
iter 0: 5	 62.5
iter 1: 2.5	 [46.875, 41.666666666666664]
iter 2: 1.25	 [42.96875, 41.666666666666664, 41.666666666666664]
Out[]: 41.666666666666664

```
