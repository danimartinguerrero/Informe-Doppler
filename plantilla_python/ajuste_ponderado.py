r"""
Plantilla de ajustes ponderados — Práctica 8: Efecto Doppler con ultrasonidos
==============================================================================

Genera los ajustes lineales del informe con el MISMO método y el MISMO estilo
gráfico que la figura del primer ajuste del documento LaTeX:

  * medias de cada serie con incertidumbre u^2 = s^2/n + Δ^2/12 (tipo A + B);
  * mínimos cuadrados ponderados con varianza efectiva u_i^2 = u^2(y_i) + a^2 u^2(x_i);
  * reescalado de incertidumbres por sqrt(chi^2/nu) cuando chi^2/nu > 1;
  * incertidumbre ampliada al 95 % con la t de Student;
  * redondeo: una cifra significativa en la incertidumbre (dos si empieza por 1).

Qué produce (en la carpeta `salida/`):
  * figura_casoN.tex    -> figura pgfplots lista para \input{} en main.tex
                           (usa los estilos `datos`, `puntoA`, `puntoB`, `ajuste`
                           y `teorica` que ya están definidos en el preámbulo);
  * resultado_casoN.tex -> caja «destacado» con el resultado redondeado;
  * vista_casoN.png     -> vista previa aproximada con matplotlib (opcional);
  * un resumen por pantalla con todos los números sin redondear.

Uso:
  1. Edita la sección «1. CONFIGURACIÓN» (normalmente basta con elegir CASO).
  2. Ejecuta:  python ajuste_ponderado.py
  3. Copia los .tex de `salida/` junto a main.tex y escribe donde corresponda:
        \input{figura_caso3}
        \input{resultado_caso3}

Casos (véase la tabla de linealizaciones del informe):
  1 -> emisor móvil:    y = v,   x = 1 - fE/fR,        pendiente = v_s
  2 -> emisor móvil:    y = fR,  x = (1 - v/ṽs)^-1,    pendiente = fE
  3 -> receptor móvil:  y = fR,  x = 1 + v/ṽs,         pendiente = fE
  4 -> receptor móvil:  Δf = fR - fE(1 - v/ṽs)^-1 frente a v/ṽs, SIN ajuste
       (término cuadrático), con la curva teórica -fE (v/vs)^2/(1 - v/vs).

Requisitos: numpy y scipy (y matplotlib solo para la vista previa).

Cómo está pensado
-----------------
Te he dejado comentarios por todo el código explicando qué hago en cada paso y,
sobre todo, por qué. No hace falta entenderlo todo para usarlo: con tocar la
sección 1 basta. Pero si algún número te sale raro, leyendo los comentarios
deberías poder seguir el cálculo de principio a fin.
"""

import math
from pathlib import Path

import numpy as np
from scipy import stats

# =============================================================================
# 1. CONFIGURACIÓN (EDITAR AQUÍ)
# =============================================================================

# Lo primero es decirle al script qué recta quieres. Para tus ajustes serán el
# caso 2 (emisor móvil) y el 3 (receptor móvil); el 4 es la gráfica del término
# cuadrático, que no lleva ajuste. El 1 es el mío, por si quieres comprobar que te
# salen los mismos números que en el informe (v_s = 244,7 m/s).
CASO = 3            # 1, 2, 3 o 4 (véase arriba)

# Modelo de la recta:
#   "origen" -> y = a x      (lo que predice la teoría; recomendado en los casos 2 y 3,
#                             porque la abscisa apenas se separa de 1 y una ordenada
#                             libre queda completamente correlacionada con la pendiente)
#   "libre"  -> y = a x + b  (permite comprobar si b es compatible con cero; caso 1)
MODELO = "origen"

# Series que se excluyen del ajuste (índices de la lista de series, empezando en 0).
# Ejemplo: EXCLUIR = [0, 1] descarta las dos series a velocidad máxima.
EXCLUIR = []

# Estos son los valores con los que comparamos los resultados. f_E es la frecuencia
# que medimos con todo en reposo (media de las seis lecturas) y ṽ_s la velocidad del
# sonido que sale de la temperatura del laboratorio. No los cambies salvo que se
# corrija algo en el informe.
# Magnitudes de referencia (valores del informe)
FE, U_FE = 40607.667, 0.357        # Hz, frecuencia de emisión en reposo (tipo A + B)
VS_REF, U_VS_REF = 344.85, 0.236   # m/s, velocidad del sonido de referencia

# Resoluciones de los instrumentos
RES_V = 0.001   # m/s
RES_F = 1.0     # Hz

# Aquí van los datos tal cual los apuntamos. Ya están metidos los de las tablas del
# informe, así que en principio no tienes que escribir nada. Si quieres cambiar
# algo, respeta el formato: seis velocidades y seis frecuencias por serie.
# Datos brutos: (selector, sentido, [6 velocidades en m/s], [6 frecuencias en Hz]).
# «sentido» se indica como "acercamiento" o "alejamiento": el script asigna el
# signo de la velocidad relativa v (positiva al acercarse emisor y receptor).
# Las velocidades se escriben en valor absoluto, tal como las da el temporizador.
DATOS_EMISOR = [   # emisor en movimiento, receptor en reposo (casos 1 y 2)
    ("Máxima", "acercamiento", [0.341, 0.358, 0.354, 0.358, 0.354, 0.354],
                               [40699, 40708, 40705, 40703, 40677, 40686]),
    ("Máxima", "alejamiento",  [0.381, 0.388, 0.382, 0.382, 0.380, 0.382],
                               [40452, 40346, 40382, 40293, 40428, 40265]),
    ("Media",  "acercamiento", [0.227, 0.222, 0.226, 0.222, 0.223, 0.225],
                               [40655, 40664, 40657, 40661, 40652, 40656]),
    ("Media",  "alejamiento",  [0.241, 0.240, 0.240, 0.244, 0.244, 0.243],
                               [40601, 40587, 40586, 40586, 40595, 40590]),
    ("Mínima", "acercamiento", [0.127, 0.128, 0.127, 0.127, 0.128, 0.126],
                               [40635, 40642, 40652, 40627, 40629, 40627]),
    ("Mínima", "alejamiento",  [0.127, 0.129, 0.131, 0.127, 0.128, 0.128],
                               [40608, 40607, 40605, 40601, 40603, 40601]),
]
DATOS_RECEPTOR = [  # emisor en reposo, receptor en movimiento (casos 3 y 4)
    # En las tablas del informe, «Negativo» = el receptor se acerca al emisor.
    ("Máxima", "acercamiento", [0.350, 0.351, 0.349, 0.347, 0.342, 0.344],
                               [40799, 40720, 40733, 40800, 40676, 40708]),
    ("Máxima", "alejamiento",  [0.378, 0.381, 0.374, 0.368, 0.365, 0.369],
                               [40429, 40520, 40438, 40500, 40557, 40489]),
    ("Media",  "acercamiento", [0.234, 0.228, 0.229, 0.226, 0.227, 0.226],
                               [40655, 40663, 40660, 40662, 40654, 40670]),
    ("Media",  "alejamiento",  [0.250, 0.251, 0.251, 0.253, 0.251, 0.252],
                               [40582, 40584, 40588, 40588, 40583, 40583]),
    ("Mínima", "acercamiento", [0.136, 0.140, 0.139, 0.140, 0.139, 0.140],
                               [40637, 40649, 40652, 40646, 40637, 40642]),
    ("Mínima", "alejamiento",  [0.147, 0.146, 0.145, 0.145, 0.146, 0.150],
                               [40596, 40611, 40599, 40598, 40598, 40610]),
]

GENERAR_PNG = True                 # vista previa con matplotlib
CARPETA_SALIDA = Path(__file__).with_name("salida")

# =============================================================================
# 2. ESTADÍSTICA BÁSICA (no hace falta tocar nada a partir de aquí)
# =============================================================================


# -----------------------------------------------------------------------------
# ¿Por qué u^2 = s^2/n + Δ^2/12?
#
# Cada punto de la gráfica es la media de seis lecturas, y su incertidumbre
# tiene dos partes que no tienen nada que ver entre sí:
#
#   1) La dispersión de las lecturas (tipo A). Si repites la medida, no sale
#      siempre lo mismo. La desviación típica s nos dice cuánto bailan las
#      lecturas sueltas, pero a nosotros nos interesa la MEDIA, que es más fiable
#      que una lectura sola: su incertidumbre es s/sqrt(n). Por eso promediar
#      seis medidas ayuda.
#
#   2) La resolución del aparato (tipo B). El contador marca de 1 Hz en 1 Hz y el
#      temporizador de 0,001 en 0,001 m/s; dentro de ese escalón no sabemos nada.
#      Lo tratamos como si el valor real pudiera estar en cualquier punto del
#      escalón con la misma probabilidad (distribución uniforme), y eso da una
#      varianza de Δ^2/12. Es lo mismo que hicimos con el termómetro.
#
# Como las dos cosas son independientes, se suman las VARIANZAS (los cuadrados),
# no las incertidumbres. Ojo, el segundo término no es un adorno: si las seis
# lecturas fueran iguales, s valdría 0 y, sin él, diríamos que la incertidumbre
# es nula, que es mentira. Es justo lo que nos pasó con f_E en reposo.
# -----------------------------------------------------------------------------
def media_con_incertidumbre(valores, resolucion):
    """Media, desviación típica y u^2 = s^2/n + Δ^2/12 (notas de incertidumbres)."""
    x = np.asarray(valores, float)
    n = len(x)
    # ddof=1 hace que divida entre n-1 y no entre n: es la desviación típica
    # «muestral», la que viene en las notas.
    s = x.std(ddof=1)
    u = math.sqrt(s**2 / n + resolucion**2 / 12)
    return x.mean(), s, u


# Para presentar los resultados: la incertidumbre con una cifra significativa,
# salvo que empiece por 1, que entonces dejamos dos (0,0118 -> 0,012, pero
# 88,6 -> 9·10^1). El valor se redondea a la misma posición decimal.
def redondea(valor, incertidumbre):
    """Redondeo del informe: 1 cifra significativa en u (2 si la primera es un 1)."""
    exp = math.floor(math.log10(abs(incertidumbre)))
    primera = int(abs(incertidumbre) / 10**exp)
    decimales = -(exp - 1) if primera == 1 else -exp
    v, u = round(valor, decimales), round(incertidumbre, decimales)
    d = max(decimales, 0)
    return f"{v:.{d}f}".replace(".", ","), f"{u:.{d}f}".replace(".", ",")


def coma(x, d):
    return f"{x:.{d}f}".replace(".", ",")


# =============================================================================
# 3. PREPARACIÓN DE LOS PUNTOS SEGÚN EL CASO
# =============================================================================

# Aquí convierto las seis lecturas de cada serie en un único punto (media ±
# incertidumbre) y le pongo signo a la velocidad. En el informe usamos siempre la
# velocidad RELATIVA v: positiva cuando emisor y receptor se acercan y negativa
# cuando se alejan. Así la misma fórmula vale para los dos montajes.
datos = DATOS_EMISOR if CASO in (1, 2) else DATOS_RECEPTOR
series = [s for i, s in enumerate(datos) if i not in EXCLUIR]

etiquetas, acerc, v, uv, f, uf = [], [], [], [], [], []
for selector, sentido, vel, frec in series:
    signo = 1 if sentido == "acercamiento" else -1
    mv, _, umv = media_con_incertidumbre(vel, RES_V)
    mf, _, umf = media_con_incertidumbre(frec, RES_F)
    etiquetas.append(f"{selector}, {sentido}")
    acerc.append(signo > 0)
    v.append(signo * mv)
    uv.append(umv)
    f.append(mf)
    uf.append(umf)
v, uv, f, uf, acerc = map(np.array, (v, uv, f, uf, acerc))

# Cada caso es una de las rectas de la tabla de linealizaciones del informe. Lo que
# cambia es qué pongo en cada eje. La abscisa x no es un dato directo: se calcula
# con magnitudes medidas, así que también tiene incertidumbre. La saco propagando
# (derivada parcial al cuadrado por la incertidumbre al cuadrado, y se suman).
#
# Un detalle de la gráfica: en los casos 2 y 3 la x vale 0,999... o 1,001...,
# y si la pintara tal cual todos los puntos caerían casi en el mismo sitio. Por eso
# en el eje represento 10^3 (x - 1), que es lo mismo pero se lee mucho mejor.
if CASO == 1:
    y, uy = v, uv
    x = 1 - FE / f
    ux = np.sqrt(U_FE**2 / f**2 + FE**2 * uf**2 / f**4)
    escala_x = lambda t: 1e3 * t                        # eje: 10^3 x
    etiqueta_x = r"$10^{3}\,x=10^{3}\,(1-\fE/\fR)$"
    etiqueta_y = r"$v$ (\si{\meter\per\second})"
    teorica = lambda t: VS_REF * t                      # v = ṽs x
    etiqueta_teo = r"$v=\vst\,x$ (valor teórico)"
    magnitud, unidad = r"\vs", r"\meter\per\second"
elif CASO == 2:
    y, uy = f, uf
    q = v / VS_REF
    x = 1 / (1 - q)
    ux = x**2 * np.sqrt((uv / VS_REF) ** 2 + (v * U_VS_REF / VS_REF**2) ** 2)
    escala_x = lambda t: 1e3 * (t - 1)                  # eje: 10^3 (x-1)
    etiqueta_x = r"$10^{3}\,\bigl[(1-v/\vst)^{-1}-1\bigr]$"
    etiqueta_y = r"$\fR$ (\si{\hertz})"
    teorica = lambda t: FE * t
    etiqueta_teo = r"$\fE$ en reposo"
    magnitud, unidad = r"\fE", r"\hertz"
elif CASO == 3:
    y, uy = f, uf
    x = 1 + v / VS_REF
    ux = np.sqrt((uv / VS_REF) ** 2 + (v * U_VS_REF / VS_REF**2) ** 2)
    escala_x = lambda t: 1e3 * (t - 1)
    etiqueta_x = r"$10^{3}\,(v/\vst)$"
    etiqueta_y = r"$\fR$ (\si{\hertz})"
    teorica = lambda t: FE * t
    etiqueta_teo = r"$\fE$ en reposo"
    magnitud, unidad = r"\fE", r"\hertz"
elif CASO == 4:
    q = v / VS_REF
    x, ux = q, np.zeros_like(q)
    y = f - FE / (1 - q)
    uy = np.sqrt(uf**2 + U_FE**2 / (1 - q) ** 2
                 + FE**2 / (VS_REF**2 * (1 - q) ** 4) * (uv**2 + q**2 * U_VS_REF**2))
    escala_x = lambda t: 1e3 * t
    etiqueta_x = r"$10^{3}\,(v/\vst)$"
    etiqueta_y = r"$\Delta f$ (\si{\hertz})"
    teorica = lambda t: -FE * t**2 / (1 - t)
    etiqueta_teo = "término cuadrático predicho"
else:
    raise SystemExit("CASO debe ser 1, 2, 3 o 4")

# =============================================================================
# 4. AJUSTE PONDERADO CON VARIANZA EFECTIVA
# =============================================================================


# -----------------------------------------------------------------------------
# El ajuste ponderado, paso a paso.
#
# No todos los puntos son igual de fiables: hay series con 1 Hz de incertidumbre
# y otras con 30 Hz. Si las tratáramos igual, un punto malo nos podría arrastrar
# la recta entera (eso me pasaba en mi ajuste antes de ponderar). Por eso cada
# punto pesa w = 1/u^2: uno con el doble de incertidumbre pesa cuatro veces menos.
#
# Además, la x también tiene error. Si mueves un punto un poco en horizontal, sobre
# una recta de pendiente a eso equivale a moverlo a·δx en vertical. Así que la
# incertidumbre que uso para el peso es la «varianza efectiva»:
#       u_i^2 = u^2(y_i) + a^2 u^2(x_i)
# El lío es que para calcularla necesito la pendiente a, que es justo lo que
# busco. Se resuelve iterando: empiezo con una pendiente aproximada, calculo los
# pesos, recalculo la pendiente, y repito hasta que deja de cambiar (en la
# práctica converge en unas pocas vueltas).
#
# Dos modelos posibles:
#   - "libre":  y = a x + b. Se ajustan dos parámetros, así que quedan n - 2
#               grados de libertad. Sirve para ver si b sale compatible con 0.
#   - "origen": y = a x, que es lo que dice la teoría. Un solo parámetro, n - 1
#               grados de libertad. En tus casos (2 y 3) es lo recomendable,
#               porque como x está pegada a 1, dejar b libre obliga a extrapolar
#               muchísimo y a y b se vuelven inseparables.
# -----------------------------------------------------------------------------
def ajuste_ponderado(x, ux, y, uy, modelo):
    """Devuelve a, u(a), b, u(b), chi2, nu e incertidumbres efectivas u_i."""
    n = len(x)
    a = np.polyfit(x, y, 1)[0] if modelo == "libre" else (x @ y) / (x @ x)
    for _ in range(500):                                  # iteración de la varianza efectiva
        ue = np.sqrt(uy**2 + a**2 * ux**2)
        w = 1 / ue**2
        if modelo == "libre":
            W = w.sum()
            m = lambda z: (w * z).sum() / W
            a_nueva = (m(x * y) - m(x) * m(y)) / (m(x * x) - m(x) ** 2)
        else:
            a_nueva = (w * x * y).sum() / (w * x * x).sum()
        if abs(a_nueva - a) <= 1e-12 * max(1.0, abs(a)):
            a = a_nueva
            break
        a = a_nueva
    ue = np.sqrt(uy**2 + a**2 * ux**2)
    w = 1 / ue**2
    if modelo == "libre":
        W = w.sum()
        m = lambda z: (w * z).sum() / W
        b = m(y) - a * m(x)
        ua = 1 / math.sqrt(W) / math.sqrt(m(x * x) - m(x) ** 2)
        ub = ua * math.sqrt(m(x * x))
        nu = n - 2
    else:
        b, ub = 0.0, 0.0
        ua = 1 / math.sqrt((w * x * x).sum())
        nu = n - 1
    chi2 = float((w * (y - a * x - b) ** 2).sum())
    return a, ua, b, ub, chi2, nu, ue


CARPETA_SALIDA.mkdir(exist_ok=True)
xs = escala_x(x)
uxs = escala_x(x + ux) - escala_x(x)

# Límites de los ejes con un margen del 12 %
def limites(centro, error):
    lo, hi = float(np.min(centro - error)), float(np.max(centro + error))
    m = 0.12 * (hi - lo if hi > lo else 1)
    return lo - m, hi + m

xmin, xmax = limites(xs, uxs)

print(f"\n=== Caso {CASO} — {len(series)} series (excluidas: {EXCLUIR or 'ninguna'}) ===")
for e, xi, uxi, yi, uyi in zip(etiquetas, xs, uxs, y, uy):
    print(f"  {e:<22s}  x = {xi:+.5f} ± {uxi:.5f}   y = {yi:.5f} ± {uyi:.5f}")

if CASO in (1, 2, 3):
    a, ua, b, ub, chi2, nu, ue = ajuste_ponderado(x, ux, y, uy, MODELO)
    # ¿Se ajusta bien la recta? Si las barras de error están bien estimadas,
    # chi2 debería salir parecido a nu. Si sale mucho mayor (a nosotros nos sale
    # 17-20 veces más), los puntos se separan de la recta más de lo que dicen sus
    # barras, y la incertidumbre de la pendiente calculada solo con los pesos se
    # queda corta. Lo honesto es inflarla multiplicando por sqrt(chi2/nu). Si chi2
    # es menor que nu, no toco nada.
    factor = max(1.0, math.sqrt(chi2 / nu)) if nu > 0 else 1.0
    ua_f, ub_f = ua * factor, ub * factor
    # Para dar el intervalo al 95 % multiplico por la t de Student, no por 2:
    # con tan pocos puntos la t es bastante mayor (2,57 con nu = 5, 4,30 con nu = 2).
    t95 = stats.t.ppf(0.975, nu)
    Ua, Ub = t95 * ua_f, t95 * ub_f
    # Residuos normalizados: cuántas «barras de error» se aparta cada punto de la
    # recta. Entre -2 y 2 es normal; si alguno se va a 5 o 6, ese punto es sospechoso.
    residuos = (y - a * x - b) / ue
    # Comparo con la referencia: si la diferencia es menor que la incertidumbre
    # ampliada, decimos que son compatibles al 95 %.
    ref = VS_REF if CASO == 1 else FE
    print(f"\n  Modelo: {'y = a x + b' if MODELO == 'libre' else 'y = a x'}   nu = {nu}")
    print(f"  a = {a:.4f}   u(a) = {ua:.4f}   chi2 = {chi2:.3f}   chi2/nu = {chi2/nu:.3f}")
    print(f"  factor sqrt(chi2/nu) aplicado = {factor:.4f}  ->  u(a) = {ua_f:.4f}")
    print(f"  t95,{nu} = {t95:.3f}  ->  U(a) = {Ua:.4f}")
    print(f"  Comparación con la referencia {ref:.3f}: |a - ref| = {abs(a-ref):.3f}  "
          f"{'<=' if abs(a-ref) <= Ua else '>'} U(a) = {Ua:.3f}  "
          f"-> {'compatibles' if abs(a-ref) <= Ua else 'NO compatibles'} al 95 %")
    if MODELO == "libre":
        print(f"  b = {b:.5f}   U(b) = {Ub:.5f}  -> "
              f"{'compatible' if abs(b) <= Ub else 'NO compatible'} con cero")
    print("  residuos normalizados r_i/u_i:", np.round(residuos, 2))
else:
    print("\n  Caso 4: sin ajuste. Valores teóricos del término cuadrático (Hz):",
          np.round(teorica(x), 3))

# =============================================================================
# 5. FIGURA PGFPLOTS (mismo estilo que el informe)
# =============================================================================


# -----------------------------------------------------------------------------
# A partir de aquí no hay física: solo escribo la figura en código pgfplots para
# que, al compilarla dentro del informe, salga con exactamente el mismo estilo que
# la mía (colores, tipografía, barras de error, panel de residuos...). Uso los
# estilos que ya están definidos en el preámbulo de main.tex, así que no hay que
# configurar nada en LaTeX: basta con \input{figura_casoN}.
# -----------------------------------------------------------------------------
def coords(px, py, ex, ey):
    return " ".join(f"({a_:.6f},{b_:.6f}) +- ({c_:.6f},{d_:.6f})"
                    for a_, b_, c_, d_ in zip(px, py, ex, ey))


ymin, ymax = limites(y, uy)


def marcas(lo, hi, n=5):
    """Marcas «redondas» (1, 2 o 5 por potencia de 10) y decimales necesarios."""
    paso = (hi - lo) / n
    e = math.floor(math.log10(paso))
    paso = min((m * 10**e for m in (1, 2, 5, 10) if m * 10**e >= paso))
    ini = math.ceil(lo / paso) * paso
    valores = np.arange(ini, hi + 1e-9 * paso, paso)
    return valores, max(0, -math.floor(math.log10(paso)))


yt, ydec = marcas(ymin, ymax)
ytick = "{" + ",".join(f"{t:.{ydec}f}" for t in yt) + "}"
dom = f"{xmin:.4f}:{xmax:.4f}"
es_frecuencia = CASO in (2, 3)
x_real = (lambda s: 1 + s / 1e3) if CASO in (2, 3) else (lambda s: s / 1e3)

# Recta teórica y ajustada, escritas en la variable del eje (pgfplots usa «x»)
if CASO == 1:
    f_teo = f"{VS_REF/1e3:.7f}*x"
elif CASO in (2, 3):
    f_teo = f"{FE:.4f}*(1+x/1000)"
else:
    f_teo = f"-{FE:.4f}*(x/1000)^2/(1-x/1000)"

lineas = [r"\begin{figure}[ht]", r"\centering", r"\begin{tikzpicture}"]
opciones_y = (rf", ytick={ytick}, yticklabel style={{/pgf/number format/fixed, "
              rf"/pgf/number format/precision={ydec}, /pgf/number format/1000 sep={{\,}}}}")
con_residuos = CASO in (1, 2, 3)
lineas.append(
    rf"  \begin{{axis}}[datos, name=principal{CASO}, width=0.86\textwidth, height=6.2cm,"
    + "\n" + rf"      xmin={xmin:.4f}, xmax={xmax:.4f}, ymin={ymin:.6g}, ymax={ymax:.6g}, clip=true,"
    + "\n" + (r"      xlabel={}, xticklabels={}," if con_residuos else rf"      xlabel={{{etiqueta_x}}},")
    + rf" ylabel={{{etiqueta_y}}}, legend pos=north west{opciones_y}]")
lineas.append(rf"    \addplot[teorica, domain={dom}, samples={2 if CASO != 4 else 60}] {{{f_teo}}};")
lineas.append(rf"    \addlegendentry{{{etiqueta_teo}}}")
if con_residuos:
    if CASO == 1:
        f_aj = f"{a/1e3:.8f}*x+({b:.8f})"
    else:
        f_aj = f"{a:.5f}*(1+x/1000)+({b:.5f})"
    lineas.append(rf"    \addplot[ajuste, domain={dom}, samples=2] {{{f_aj}}};")
    lineas.append(r"    \addlegendentry{ajuste ponderado}")
if acerc.any():
    lineas.append(r"    \addplot[puntoA] coordinates {" + coords(xs[acerc], y[acerc], uxs[acerc], uy[acerc]) + "};")
    lineas.append(r"    \addlegendentry{acercamiento}")
if (~acerc).any():
    lineas.append(r"    \addplot[puntoB] coordinates {" + coords(xs[~acerc], y[~acerc], uxs[~acerc], uy[~acerc]) + "};")
    lineas.append(r"    \addlegendentry{alejamiento}")
lineas.append(r"  \end{axis}")
if con_residuos:
    rmax = max(3.2, 1.25 * float(np.max(np.abs(residuos) + 1)))
    lineas += [
        rf"  \begin{{axis}}[datos, at={{(principal{CASO}.south west)}}, anchor=north west, yshift=-4pt,",
        rf"      width=0.86\textwidth, height=3.0cm, xmin={xmin:.4f}, xmax={xmax:.4f},",
        rf"      ymin={-rmax:.2f}, ymax={rmax:.2f}, xlabel={{{etiqueta_x}}}, ylabel={{$r_i/u_i$}}]",
        rf"    \addplot[grispizarra, line width=0.5pt, densely dashed, domain={dom}, samples=2] {{0}};",
        r"    \addplot[only marks, mark=*, mark size=1.7pt, azulprof,",
        r"      error bars/.cd, y dir=both, y explicit,",
        r"      error bar style={line width=0.5pt, azulprof}, error mark options={rotate=90, mark size=1.8pt, azulprof}]",
        "      coordinates {" + " ".join(f"({p:.6f},{r:.4f}) +- (0,1)" for p, r in zip(xs, residuos)) + "};",
        r"  \end{axis}",
    ]
lineas.append(r"\end{tikzpicture}")
pie = {
    1: r"Velocidad relativa del emisor frente a $1-\fE/\fR$.",
    2: r"Frecuencia recibida frente a $(1-v/\vst)^{-1}$ con el emisor en movimiento.",
    3: r"Frecuencia recibida frente a $1+v/\vst$ con el receptor en movimiento.",
    4: r"Diferencias $\Delta f=\fR-\fE\,(1-v/\vst)^{-1}$ con el receptor en movimiento, sin ajuste.",
}[CASO]
pie += (r" Círculos: acercamiento; cuadrados: alejamiento; barras: incertidumbre típica.")
if con_residuos:
    pie += (r" Recta continua: ajuste ponderado; discontinua: predicción teórica."
            r" Debajo, los residuos normalizados $r_i/u_i$.")
else:
    pie += r" Curva discontinua: término cuadrático predicho."
lineas += [rf"\caption{{{pie}}}", rf"\label{{fig:caso{CASO}}}", r"\end{figure}", ""]
(CARPETA_SALIDA / f"figura_caso{CASO}.tex").write_text("\n".join(lineas), encoding="utf-8")

# Caja de resultado (entorno «destacado» del informe)
if con_residuos:
    va, Ua_s = redondea(a, Ua)
    caja = [rf"\begin{{destacado}}[Ajuste del caso {CASO}]", r"\vspace{-4pt}", r"\begin{equation}",
            rf"  {magnitud}=(\num{{{va}}}\pm\num{{{Ua_s}}})\,\si{{{unidad}}}"]
    if MODELO == "libre":
        vb, Ub_s = redondea(b, Ub)
        uni_b = r"\meter\per\second" if CASO == 1 else r"\hertz"
        caja[-1] += rf",\qquad b=(\num{{{vb}}}\pm\num{{{Ub_s}}})\,\si{{{uni_b}}}"
    caja[-1] += "."
    caja += [r"\end{equation}",
             rf"{{\small Incertidumbres ampliadas al \SI{{95}}{{\percent}}: $\nu={nu}$, "
             rf"$t_{{95,{nu}}}=\num{{{coma(t95, 2)}}}$, $\chi^{{2}}/\nu=\num{{{coma(chi2/nu, 1)}}}$.}}",
             r"\end{destacado}", ""]
    (CARPETA_SALIDA / f"resultado_caso{CASO}.tex").write_text("\n".join(caja), encoding="utf-8")

print(f"\n  Escrito: {CARPETA_SALIDA / f'figura_caso{CASO}.tex'}")
if con_residuos:
    print(f"  Escrito: {CARPETA_SALIDA / f'resultado_caso{CASO}.tex'}")

# =============================================================================
# 6. VISTA PREVIA CON MATPLOTLIB (misma paleta; tipografía aproximada)
# =============================================================================

# La vista previa en PNG es solo para que veas rápido cómo queda sin compilar el
# informe. La que va al documento es la del .tex; esta se parece pero no es idéntica.
if GENERAR_PNG:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.ticker import FuncFormatter
    except ImportError:
        print("  (matplotlib no está instalado: se omite la vista previa)")
    else:
        AZULPROF, AZULACERO, AZULCLARO = "#1B365D", "#3B6EA8", "#B7CBE3"
        GRIS, GRISCLARO, ACENTO = "#5B6770", "#E4E8EC", "#B8800A"
        plt.rcParams.update({"font.family": "serif",
                             "font.serif": ["Palatino", "TeX Gyre Pagella", "Palatino Linotype",
                                            "URW Palladio L", "DejaVu Serif"],
                             "mathtext.fontset": "dejavuserif", "font.size": 9,
                             "axes.edgecolor": GRIS, "axes.labelcolor": "black",
                             "xtick.color": GRIS, "ytick.color": GRIS})
        fmt = FuncFormatter(lambda t, _: f"{t:g}".replace(".", ","))
        filas = 2 if con_residuos else 1
        fig, ejes = plt.subplots(filas, 1, figsize=(6.2, 4.3 if con_residuos else 3.4),
                                 sharex=True, gridspec_kw={"height_ratios": [3, 1.2]} if con_residuos else None)
        ax = ejes[0] if con_residuos else ejes
        malla = np.linspace(xmin, xmax, 200)
        ax.plot(malla, teorica(x_real(malla)), "--", color=ACENTO, lw=1.0, label="teórica")
        if con_residuos:
            ax.plot(malla, a * x_real(malla) + b, color=AZULPROF, lw=1.3, label="ajuste ponderado")
        for mascara, marcador, color, relleno, nombre in ((acerc, "o", AZULPROF, AZULPROF, "acercamiento"),
                                                          (~acerc, "s", AZULACERO, AZULCLARO, "alejamiento")):
            if mascara.any():
                ax.errorbar(xs[mascara], y[mascara], xerr=uxs[mascara], yerr=uy[mascara], fmt=marcador,
                            ms=4.5, color=color, mfc=relleno, mec=color, ecolor=color, elinewidth=0.7,
                            capsize=2.2, label=nombre)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.grid(color=GRISCLARO, lw=0.5)
        ax.yaxis.set_major_formatter(fmt)
        ax.legend(frameon=False, fontsize=7, loc="upper left")
        for lado in ("top", "right"):
            ax.spines[lado].set_visible(False)
        if con_residuos:
            ax2 = ejes[1]
            ax2.axhline(0, color=GRIS, lw=0.6, ls="--")
            ax2.errorbar(xs, residuos, yerr=1, fmt="o", ms=3.5, color=AZULPROF, capsize=2, elinewidth=0.6)
            ax2.grid(color=GRISCLARO, lw=0.5)
            ax2.set_ylabel(r"$r_i/u_i$")
            for lado in ("top", "right"):
                ax2.spines[lado].set_visible(False)
            ax2.xaxis.set_major_formatter(fmt)
            ax_x = ax2
        else:
            ax.xaxis.set_major_formatter(fmt)
            ax_x = ax
        limpia = lambda t: (t.replace(r"\fE", "f_E").replace(r"\fR", "f_R").replace(r"\vst", r"\tilde{v}_s")
                            .replace(r"\,", " ").replace(r"\bigl", "").replace(r"\bigr", "")
                            .replace(r"(\si{\meter\per\second})", "(m/s)").replace(r"(\si{\hertz})", "(Hz)"))
        ax_x.set_xlabel(limpia(etiqueta_x))
        ax.set_ylabel(limpia(etiqueta_y))
        fig.tight_layout()
        ruta = CARPETA_SALIDA / f"vista_caso{CASO}.png"
        fig.savefig(ruta, dpi=200)
        print(f"  Escrito: {ruta}")
