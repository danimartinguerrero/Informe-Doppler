# Plantillas de gráficas con ajuste ponderado (Práctica 8, efecto Doppler)

Para hacer gráficas con **el mismo estilo que las del informe** (colores, tipografía,
barras de error, panel de residuos y caja azulada con el resultado).

Instala lo necesario una vez: `pip install numpy scipy matplotlib`

## `mi_grafica.py`: para meter tus propios datos

1. Abre `mi_grafica.py` y rellena la sección **TUS DATOS**:
   - los puntos: `X`, `UX`, `Y`, `UY`, escritos a mano o calculados a partir de las
     lecturas brutas con `media(lecturas, resolucion)`;
   - si quieres, `GRUPOS` (0 = círculos, 1 = cuadrados) y sus nombres en la leyenda;
   - `MODELO`: `"libre"` (y = a x + b), `"origen"` (y = a x) o `None` (sin ajuste);
   - una curva teórica opcional, la referencia con la que comparar y los textos
     (ejes, pie de figura, título de la caja).
   Viene rellenado con un ejemplo para que veas cómo se escribe todo.
2. Ejecuta `python mi_grafica.py`.
3. En `salida/` aparecen `NOMBRE.tex` (la figura), `NOMBRE_resultado.tex` (la caja con
   el resultado) y `NOMBRE.png` (vista previa). Por pantalla salen todos los números.
4. Sube los `.tex` a Overleaf, junto a `main.tex`, y escribe donde quieras la figura:
   `\input{NOMBRE}` y `\input{NOMBRE_resultado}`.

Los cálculos están en `herramientas.py`, comentados; no hace falta tocarlo.

## `ajuste_ponderado.py`: los casos del informe ya preparados

Si solo quieres rehacer uno de los ajustes del informe (casos 1–4) con los datos de las
tablas, basta con elegir `CASO` y `MODELO` al principio y ejecutarlo.
