r"""
Plantilla para hacer TU gráfica con ajuste ponderado y el estilo del informe
=============================================================================

Rellena la sección «TUS DATOS» y ejecuta:

    python mi_grafica.py

En la carpeta `salida/` te aparecerán:
  * <NOMBRE>.tex           -> la figura, lista para \input{<NOMBRE>} en main.tex
  * <NOMBRE>_resultado.tex -> la caja azulada con el resultado redondeado
  * <NOMBRE>.png           -> una vista previa rápida (la buena es la del .tex)

y por pantalla todos los números sin redondear para que puedas comentarlos.

Lo he dejado rellenado con un ejemplo (el ajuste f_R frente a 1 + v/ṽ_s con el
receptor en movimiento) para que veas cómo se escribe cada cosa. Cámbialo por lo
tuyo. Los cálculos están en `herramientas.py`, explicados; no hace falta tocarlo.
"""

from pathlib import Path

import numpy as np

from herramientas import (media, ajuste_ponderado, resumen, figura_tex,
                          caja_resultado, vista_previa)

# =============================================================================
# TUS DATOS
# =============================================================================

# Nombre de los archivos que se generan (sin extensión). Usa uno distinto para
# cada gráfica, así no se pisan.
NOMBRE = "ajuste_receptor"

# -----------------------------------------------------------------------------
# 1) Los puntos.
#
# Necesito, para cada punto: x, su incertidumbre ux, y, y su incertidumbre uy.
# Tienes dos formas de dármelos:
#
#   a) Si ya los tienes calculados, escríbelos directamente en las listas X, UX,
#      Y, UY de más abajo (todas con la misma longitud).
#   b) Si tienes las lecturas brutas, usa media(lecturas, resolucion): te
#      devuelve la media y su incertidumbre (u^2 = s^2/n + Δ^2/12). Es lo que hago
#      en el ejemplo.
#
# Si alguna x no tiene incertidumbre, pon 0.
# -----------------------------------------------------------------------------

# --- Ejemplo con lecturas brutas (receptor en movimiento) --------------------
FE = 40607.667             # Hz, frecuencia del emisor medida en reposo
VS_REF = 344.85            # m/s, velocidad del sonido de referencia
U_VS_REF = 0.236           # m/s, su incertidumbre

# (velocidades en m/s, frecuencias en Hz, ¿se acercan?)
series = [
    ([0.350, 0.351, 0.349, 0.347, 0.342, 0.344], [40799, 40720, 40733, 40800, 40676, 40708], True),
    ([0.378, 0.381, 0.374, 0.368, 0.365, 0.369], [40429, 40520, 40438, 40500, 40557, 40489], False),
    ([0.234, 0.228, 0.229, 0.226, 0.227, 0.226], [40655, 40663, 40660, 40662, 40654, 40670], True),
    ([0.250, 0.251, 0.251, 0.253, 0.251, 0.252], [40582, 40584, 40588, 40588, 40583, 40583], False),
    ([0.136, 0.140, 0.139, 0.140, 0.139, 0.140], [40637, 40649, 40652, 40646, 40637, 40642], True),
    ([0.147, 0.146, 0.145, 0.145, 0.146, 0.150], [40596, 40611, 40599, 40598, 40598, 40610], False),
]

X, UX, Y, UY, GRUPOS = [], [], [], [], []
for velocidades, frecuencias, se_acercan in series:
    v, uv = media(velocidades, resolucion=0.001)
    f, uf = media(frecuencias, resolucion=1)
    v = v if se_acercan else -v                  # velocidad relativa con signo
    # La abscisa es x = 1 + v/ṽ_s. Como depende de v y de ṽ_s, propago las dos:
    # u^2(x) = (uv/ṽ_s)^2 + (v·u(ṽ_s)/ṽ_s^2)^2
    X.append(1 + v / VS_REF)
    UX.append(np.hypot(uv / VS_REF, v * U_VS_REF / VS_REF**2))
    Y.append(f)
    UY.append(uf)
    GRUPOS.append(0 if se_acercan else 1)        # 0 = círculos, 1 = cuadrados

# --- Si prefieres meter los puntos a mano, borra el ejemplo y usa esto: -------
# X  = [ ... ]
# UX = [ ... ]
# Y  = [ ... ]
# UY = [ ... ]
# GRUPOS = None      # o una lista de 0 y 1 para separar dos tipos de punto

# Nombres de los dos tipos de punto en la leyenda ("" para no ponerlos)
NOMBRES_GRUPOS = ("acercamiento", "alejamiento")

# Puntos que quieres dejar fuera del ajuste (índices empezando en 0). Siguen
# saliendo en la tabla de números, pero no en la gráfica ni en la recta.
EXCLUIR = []

# -----------------------------------------------------------------------------
# 2) El ajuste.
#
#   "libre"  -> y = a x + b  (dos parámetros; te dice si b es compatible con 0)
#   "origen" -> y = a x      (si la teoría dice que la recta pasa por el origen)
#   None     -> no ajusta nada, solo pinta los puntos (y la curva teórica)
# -----------------------------------------------------------------------------
MODELO = "origen"

# Valor con el que comparar la pendiente (o None si no hay)
REFERENCIA = FE

# Curva teórica a dibujar en ámbar discontinuo (o None). Es una función de la x
# real (no de la escalada). Ejemplo: la recta que saldría con f_E en reposo.
TEORICA = lambda x: FE * x
NOMBRE_TEORICA = r"$\fE$ en reposo"

# -----------------------------------------------------------------------------
# 3) Cómo se ve el eje x.
#
# Si tus x son números muy parecidos (aquí 0,999... y 1,001...), en el eje se
# pinta escala·(x − origen) para que se lea bien. Si no te hace falta, deja
# ESCALA_X = 1 y ORIGEN_X = 0.
# -----------------------------------------------------------------------------
ESCALA_X = 1e3
ORIGEN_X = 1.0

# -----------------------------------------------------------------------------
# 4) Textos. En los .tex puedes usar las macros del informe (\fE, \fR, \vst, \si...).
# -----------------------------------------------------------------------------
ETIQUETA_X = r"$10^{3}\,(v/\vst)$"
ETIQUETA_Y = r"$\fR$ (\si{\hertz})"
PIE = (r"Frecuencia recibida frente a $1+v/\vst$ con el receptor en movimiento. "
       r"Círculos: acercamiento; cuadrados: alejamiento; barras: incertidumbre típica. "
       r"Recta continua: ajuste ponderado; discontinua: predicción con la $\fE$ en reposo. "
       r"Debajo, los residuos normalizados $r_i/u_i$.")
ETIQUETA_LATEX = "fig:ajuste-receptor"          # para \cref{...}

# Caja de resultado: símbolo de la pendiente y unidades (formato siunitx)
SIMBOLO = r"\fE"
UNIDAD = r"\hertz"
UNIDAD_B = r"\hertz"                            # solo si MODELO = "libre"
TITULO_CAJA = "Ajuste 3 --- receptor móvil"

# Etiquetas sencillas para la vista previa en PNG (sin macros de LaTeX)
ETIQUETA_X_PNG = "10³ (v / ṽs)"
ETIQUETA_Y_PNG = "fR (Hz)"

# =============================================================================
# A PARTIR DE AQUÍ NO HACE FALTA TOCAR NADA
# =============================================================================

salida = Path(__file__).with_name("salida")
salida.mkdir(exist_ok=True)

X, UX, Y, UY = (np.asarray(t, float) for t in (X, UX, Y, UY))
GRUPOS = np.zeros(len(X), int) if GRUPOS is None else np.asarray(GRUPOS)

print(f"\n=== {NOMBRE}: {len(X)} puntos ===")
for i, (xi, uxi, yi, uyi) in enumerate(zip(X, UX, Y, UY)):
    marca = "  (excluido)" if i in EXCLUIR else ""
    print(f"  [{i}] x = {xi:.6g} ± {uxi:.3g}    y = {yi:.6g} ± {uyi:.3g}{marca}")

usar = np.array([i not in EXCLUIR for i in range(len(X))])
X, UX, Y, UY, GRUPOS = X[usar], UX[usar], Y[usar], UY[usar], GRUPOS[usar]

ajuste = None
if MODELO is not None:
    ajuste = ajuste_ponderado(X, UX, Y, UY, MODELO)
    resumen(ajuste, REFERENCIA)

xs, uxs = figura_tex(salida / f"{NOMBRE}.tex", X, UX, Y, UY, grupos=GRUPOS,
                     nombres_grupos=NOMBRES_GRUPOS, ajuste=ajuste, teorica=TEORICA,
                     nombre_teorica=NOMBRE_TEORICA, escala_x=ESCALA_X, origen_x=ORIGEN_X,
                     etiqueta_x=ETIQUETA_X, etiqueta_y=ETIQUETA_Y, pie=PIE,
                     etiqueta=ETIQUETA_LATEX)
print(f"\n  Escrito: {salida / (NOMBRE + '.tex')}")
if ajuste is not None:
    caja_resultado(salida / f"{NOMBRE}_resultado.tex", ajuste, SIMBOLO, UNIDAD, UNIDAD_B, TITULO_CAJA)
    print(f"  Escrito: {salida / (NOMBRE + '_resultado.tex')}")
vista_previa(salida / f"{NOMBRE}.png", xs, uxs, Y, UY, grupos=GRUPOS, nombres_grupos=NOMBRES_GRUPOS,
             ajuste=ajuste, teorica=TEORICA, escala_x=ESCALA_X, origen_x=ORIGEN_X,
             etiqueta_x=ETIQUETA_X_PNG, etiqueta_y=ETIQUETA_Y_PNG)
print(f"  Escrito: {salida / (NOMBRE + '.png')}")
