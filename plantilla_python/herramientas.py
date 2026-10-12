r"""
Herramientas comunes para hacer ajustes y gráficas con el estilo del informe.

No hace falta tocar nada de este archivo: lo usa `mi_grafica.py`. Aquí están:

  * media(...)            -> media de unas lecturas y su incertidumbre (tipo A + B)
  * ajuste_ponderado(...) -> recta por mínimos cuadrados ponderados
  * redondea(...)         -> redondeo con el criterio del informe
  * figura_tex(...)       -> escribe la figura en pgfplots, con nuestro estilo
  * caja_resultado(...)   -> escribe la caja azulada con el resultado
  * vista_previa(...)     -> PNG rápido con matplotlib, para ver cómo queda
"""

import math
from pathlib import Path

import numpy as np
from scipy import stats


# -----------------------------------------------------------------------------
# Media de unas lecturas repetidas
#
# Cada punto de la gráfica suele ser la media de varias lecturas. Su incertidumbre
# tiene dos partes independientes:
#   - lo que «bailan» las lecturas (tipo A): s/sqrt(n), la incertidumbre de la media;
#   - el escalón del aparato (tipo B): si marca de Δ en Δ, el valor real puede
#     estar en cualquier punto del escalón, lo que da una varianza Δ^2/12.
# Como son independientes, se suman los cuadrados: u^2 = s^2/n + Δ^2/12.
# El segundo término importa: si todas las lecturas fueran iguales, s = 0, pero
# eso no quiere decir que la medida sea perfecta.
# -----------------------------------------------------------------------------
def media(lecturas, resolucion):
    """Devuelve (media, incertidumbre típica) de una lista de lecturas."""
    x = np.asarray(lecturas, float)
    n = len(x)
    s = x.std(ddof=1) if n > 1 else 0.0          # ddof=1: dividimos entre n-1
    return float(x.mean()), math.sqrt(s**2 / n + resolucion**2 / 12)


# -----------------------------------------------------------------------------
# Redondeo del informe: la incertidumbre con una cifra significativa, salvo que
# empiece por 1 (entonces dos). El valor, a la misma posición decimal.
# -----------------------------------------------------------------------------
def redondea(valor, incertidumbre):
    exp = math.floor(math.log10(abs(incertidumbre)))
    primera = int(abs(incertidumbre) / 10**exp)
    decimales = -(exp - 1) if primera == 1 else -exp
    v, u = round(valor, decimales), round(incertidumbre, decimales)
    d = max(decimales, 0)
    return f"{v:.{d}f}".replace(".", ","), f"{u:.{d}f}".replace(".", ",")


def _coma(x, d):
    return f"{x:.{d}f}".replace(".", ",")


# -----------------------------------------------------------------------------
# Ajuste ponderado con varianza efectiva
#
# Cada punto pesa w = 1/u^2, así los más precisos mandan más en la recta.
# Si la x también tiene error, un desplazamiento horizontal δx equivale a uno
# vertical a·δx, así que el peso se calcula con u^2 = u^2(y) + a^2 u^2(x).
# Como a sale del propio ajuste, se itera hasta que deja de cambiar.
#
# Después miramos chi^2: si es mucho mayor que los grados de libertad, los puntos
# se separan de la recta más de lo que dicen sus barras, y entonces inflamos la
# incertidumbre multiplicando por sqrt(chi^2/nu). Por último, el intervalo al 95 %
# se da con la t de Student (con pocos puntos es bastante mayor que 2).
# -----------------------------------------------------------------------------
def ajuste_ponderado(x, ux, y, uy, modelo="libre"):
    """Ajusta y = a x + b ("libre") o y = a x ("origen"). Devuelve un diccionario."""
    x, ux, y, uy = (np.asarray(t, float) for t in (x, ux, y, uy))
    n = len(x)
    a = np.polyfit(x, y, 1)[0] if modelo == "libre" else (x @ y) / (x @ x)
    for _ in range(500):
        w = 1 / (uy**2 + a**2 * ux**2)
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
        b, ub, nu = 0.0, 0.0, n - 1
        ua = 1 / math.sqrt((w * x * x).sum())
    chi2 = float((w * (y - a * x - b) ** 2).sum())
    factor = max(1.0, math.sqrt(chi2 / nu)) if nu > 0 else 1.0
    t95 = float(stats.t.ppf(0.975, nu)) if nu > 0 else float("nan")
    return dict(a=a, b=b, ua=ua * factor, ub=ub * factor, Ua=t95 * ua * factor,
                Ub=t95 * ub * factor, chi2=chi2, nu=nu, factor=factor, t95=t95,
                ue=ue, residuos=(y - a * x - b) / ue, modelo=modelo)


def resumen(r, referencia=None, nombre="a"):
    """Imprime los números del ajuste sin redondear."""
    print(f"\n  Modelo: {'y = a x + b' if r['modelo'] == 'libre' else 'y = a x'}   nu = {r['nu']}")
    print(f"  {nombre} = {r['a']:.6g}   u = {r['ua']:.4g}   U(95 %) = {r['Ua']:.4g}"
          f"   (t = {r['t95']:.3f}, factor sqrt(chi2/nu) = {r['factor']:.3f})")
    print(f"  chi2 = {r['chi2']:.3f}   chi2/nu = {r['chi2'] / r['nu']:.3f}")
    if r["modelo"] == "libre":
        ok = abs(r["b"]) <= r["Ub"]
        print(f"  b = {r['b']:.6g}   U(95 %) = {r['Ub']:.4g}   -> "
              f"{'compatible' if ok else 'NO compatible'} con cero")
    if referencia is not None:
        ok = abs(r["a"] - referencia) <= r["Ua"]
        print(f"  Referencia {referencia:.6g}: diferencia {abs(r['a'] - referencia):.4g} -> "
              f"{'compatibles' if ok else 'NO compatibles'} al 95 %")
    print("  Residuos normalizados:", np.round(r["residuos"], 2))


# -----------------------------------------------------------------------------
# Figura en pgfplots con el estilo del informe
# -----------------------------------------------------------------------------
def _marcas(lo, hi, n=5):
    paso = (hi - lo) / n
    e = math.floor(math.log10(paso))
    paso = min(m * 10**e for m in (1, 2, 5, 10) if m * 10**e >= paso)
    valores = np.arange(math.ceil(lo / paso) * paso, hi + 1e-9 * paso, paso)
    return valores, max(0, -math.floor(math.log10(paso)))


def _limites(c, e):
    lo, hi = float(np.min(c - e)), float(np.max(c + e))
    m = 0.12 * (hi - lo if hi > lo else 1)
    return lo - m, hi + m


def figura_tex(ruta, x, ux, y, uy, *, grupos=None, nombres_grupos=("", ""),
               ajuste=None, teorica=None, nombre_teorica="valor teórico",
               escala_x=1.0, origen_x=0.0, etiqueta_x="$x$", etiqueta_y="$y$",
               pie="", etiqueta="fig:mi-grafica", residuos=True):
    """Escribe la figura en `ruta` (.tex) usando los estilos de main.tex."""
    x, ux, y, uy = (np.asarray(t, float) for t in (x, ux, y, uy))
    grupos = np.zeros(len(x), int) if grupos is None else np.asarray(grupos)
    xs = escala_x * (x - origen_x)                 # lo que se dibuja en el eje
    uxs = abs(escala_x) * ux
    xmin, xmax = _limites(xs, uxs)
    ymin, ymax = _limites(y, uy)
    yt, ydec = _marcas(ymin, ymax)
    ytick = ",".join(f"{t:.{ydec}f}" for t in yt)
    real = lambda s: s / escala_x + origen_x       # del eje a la x verdadera
    malla = np.linspace(xmin, xmax, 80)
    curva = lambda f: " ".join(f"({s:.6f},{f(real(s)):.8g})" for s in malla)
    puntos = lambda m: " ".join(f"({a:.6f},{b:.8g}) +- ({c:.6f},{d:.8g})"
                                for a, b, c, d in zip(xs[m], y[m], uxs[m], uy[m]))
    con_res = residuos and ajuste is not None
    L = [r"\begin{figure}[ht]", r"\centering", r"\begin{tikzpicture}",
         rf"  \begin{{axis}}[datos, name=principal, width=0.86\textwidth, height=6.2cm,",
         rf"      xmin={xmin:.6g}, xmax={xmax:.6g}, ymin={ymin:.8g}, ymax={ymax:.8g}, clip=true,",
         (r"      xlabel={}, xticklabels={}," if con_res else rf"      xlabel={{{etiqueta_x}}},")
         + rf" ylabel={{{etiqueta_y}}}, legend pos=north west, ytick={{{ytick}}},",
         rf"      yticklabel style={{/pgf/number format/fixed, /pgf/number format/precision={ydec}}}]"]
    if teorica is not None:
        L += [rf"    \addplot[teorica] coordinates {{{curva(teorica)}}};",
              rf"    \addlegendentry{{{nombre_teorica}}}"]
    if ajuste is not None:
        a, b = ajuste["a"], ajuste["b"]
        L += [rf"    \addplot[ajuste] coordinates {{{curva(lambda t: a * t + b)}}};",
              r"    \addlegendentry{ajuste ponderado}"]
    for g, estilo in ((0, "puntoA"), (1, "puntoB")):
        m = grupos == g
        if m.any():
            L.append(rf"    \addplot[{estilo}] coordinates {{{puntos(m)}}};")
            if nombres_grupos[g]:
                L.append(rf"    \addlegendentry{{{nombres_grupos[g]}}}")
    L.append(r"  \end{axis}")
    if con_res:
        r = ajuste["residuos"]
        rmax = max(3.2, 1.25 * float(np.max(np.abs(r) + 1)))
        L += [r"  \begin{axis}[datos, at={(principal.south west)}, anchor=north west, yshift=-4pt,",
              rf"      width=0.86\textwidth, height=3.0cm, xmin={xmin:.6g}, xmax={xmax:.6g},",
              rf"      ymin={-rmax:.2f}, ymax={rmax:.2f}, xlabel={{{etiqueta_x}}}, ylabel={{$r_i/u_i$}}]",
              rf"    \addplot[grispizarra, line width=0.5pt, densely dashed] coordinates {{({xmin:.6g},0) ({xmax:.6g},0)}};",
              r"    \addplot[only marks, mark=*, mark size=1.7pt, azulprof,",
              r"      error bars/.cd, y dir=both, y explicit,",
              r"      error bar style={line width=0.5pt, azulprof}, error mark options={rotate=90, mark size=1.8pt, azulprof}]",
              "      coordinates {" + " ".join(f"({s:.6f},{q:.4f}) +- (0,1)" for s, q in zip(xs, r)) + "};",
              r"  \end{axis}"]
    L += [r"\end{tikzpicture}", rf"\caption{{{pie}}}", rf"\label{{{etiqueta}}}", r"\end{figure}", ""]
    Path(ruta).write_text("\n".join(L), encoding="utf-8")
    return xs, uxs


def caja_resultado(ruta, ajuste, simbolo="a", unidad="", unidad_b="", titulo="Resultado"):
    """Escribe la caja azulada (entorno `destacado`) con el resultado redondeado."""
    va, Ua = redondea(ajuste["a"], ajuste["Ua"])
    linea = rf"  {simbolo}=(\num{{{va}}}\pm\num{{{Ua}}})\,\si{{{unidad}}}"
    if ajuste["modelo"] == "libre":
        vb, Ub = redondea(ajuste["b"], ajuste["Ub"])
        linea += rf",\qquad b=(\num{{{vb}}}\pm\num{{{Ub}}})\,\si{{{unidad_b}}}"
    nu = ajuste["nu"]
    L = [rf"\begin{{destacado}}[{titulo}]", r"\vspace{-4pt}", r"\begin{equation}", linea + ".",
         r"\end{equation}",
         rf"{{\small Incertidumbres ampliadas al \SI{{95}}{{\percent}}: $\nu={nu}$, "
         rf"$t_{{95,{nu}}}=\num{{{_coma(ajuste['t95'], 2)}}}$, "
         rf"$\chi^{{2}}/\nu=\num{{{_coma(ajuste['chi2'] / nu, 1)}}}$.}}",
         r"\end{destacado}", ""]
    Path(ruta).write_text("\n".join(L), encoding="utf-8")


# -----------------------------------------------------------------------------
# Vista previa con matplotlib (misma paleta; la buena es la del .tex)
# -----------------------------------------------------------------------------
def vista_previa(ruta, xs, uxs, y, uy, *, grupos=None, nombres_grupos=("", ""),
                 ajuste=None, teorica=None, escala_x=1.0, origen_x=0.0,
                 etiqueta_x="x", etiqueta_y="y"):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.ticker import FuncFormatter
    except ImportError:
        print("  (matplotlib no está instalado: no hago la vista previa)")
        return
    AZULPROF, AZULACERO, AZULCLARO = "#1B365D", "#3B6EA8", "#B7CBE3"
    GRIS, GRISCLARO, ACENTO = "#5B6770", "#E4E8EC", "#B8800A"
    plt.rcParams.update({"font.family": "serif",
                         "font.serif": ["Palatino", "TeX Gyre Pagella", "URW Palladio L", "DejaVu Serif"],
                         "font.size": 9, "axes.edgecolor": GRIS, "xtick.color": GRIS, "ytick.color": GRIS})
    xs, uxs, y, uy = (np.asarray(t, float) for t in (xs, uxs, y, uy))
    grupos = np.zeros(len(xs), int) if grupos is None else np.asarray(grupos)
    con_res = ajuste is not None
    fig, ejes = plt.subplots(2 if con_res else 1, 1, figsize=(6.2, 4.3 if con_res else 3.4), sharex=True,
                             gridspec_kw={"height_ratios": [3, 1.2]} if con_res else None)
    ax = ejes[0] if con_res else ejes
    lo, hi = _limites(xs, uxs)
    malla = np.linspace(lo, hi, 200)
    real = malla / escala_x + origen_x
    if teorica is not None:
        ax.plot(malla, [teorica(t) for t in real], "--", color=ACENTO, lw=1.0, label="teórica")
    if con_res:
        ax.plot(malla, ajuste["a"] * real + ajuste["b"], color=AZULPROF, lw=1.3, label="ajuste ponderado")
    for g, (mk, col, rel) in enumerate((("o", AZULPROF, AZULPROF), ("s", AZULACERO, AZULCLARO))):
        m = grupos == g
        if m.any():
            ax.errorbar(xs[m], y[m], xerr=uxs[m], yerr=uy[m], fmt=mk, ms=4.5, color=col, mfc=rel,
                        mec=col, ecolor=col, elinewidth=0.7, capsize=2.2, label=nombres_grupos[g] or None)
    fmt = FuncFormatter(lambda t, _: f"{t:g}".replace(".", ","))
    ax.yaxis.set_major_formatter(fmt)
    ax.grid(color=GRISCLARO, lw=0.5)
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    ax.set_ylabel(etiqueta_y)
    for lado in ("top", "right"):
        ax.spines[lado].set_visible(False)
    eje_x = ax
    if con_res:
        ax2 = ejes[1]
        ax2.axhline(0, color=GRIS, lw=0.6, ls="--")
        ax2.errorbar(xs, ajuste["residuos"], yerr=1, fmt="o", ms=3.5, color=AZULPROF, capsize=2, elinewidth=0.6)
        ax2.grid(color=GRISCLARO, lw=0.5)
        ax2.set_ylabel("r/u")
        for lado in ("top", "right"):
            ax2.spines[lado].set_visible(False)
        eje_x = ax2
    eje_x.xaxis.set_major_formatter(fmt)
    eje_x.set_xlabel(etiqueta_x)
    fig.tight_layout()
    fig.savefig(ruta, dpi=200)
    plt.close(fig)
