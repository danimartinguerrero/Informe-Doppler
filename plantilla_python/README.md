# Plantilla de ajustes ponderados (Práctica 8, efecto Doppler)

`ajuste_ponderado.py` hace los ajustes lineales del informe con el mismo método
que el primer ajuste y escribe las figuras en pgfplots con el mismo estilo.

## Uso

1. Instala lo necesario: `pip install numpy scipy matplotlib`
2. Abre `ajuste_ponderado.py` y, en la sección **1. CONFIGURACIÓN**, elige:
   - `CASO`: 2 (emisor móvil, frecuencia), 3 (receptor móvil) o 4 (término cuadrático);
   - `MODELO`: `"origen"` (recta por el origen, recomendado en los casos 2 y 3) o `"libre"`;
   - `EXCLUIR`: series que quieras descartar, por ejemplo `[0, 1]` para las de velocidad máxima.
   Los datos de las tablas del informe ya están cargados; solo hay que cambiarlos si se corrige alguno.
3. Ejecuta `python ajuste_ponderado.py`.
4. En la carpeta `salida/` aparecen:
   - `figura_casoN.tex`: la figura, lista para `\input{figura_casoN}` en `main.tex`;
   - `resultado_casoN.tex`: la caja azul con el resultado redondeado;
   - `vista_casoN.png`: una vista previa rápida.
   Por pantalla se imprimen todos los números sin redondear (pendiente, χ², factor
   √(χ²/ν), incertidumbre ampliada y compatibilidad con la referencia).

Copia los `.tex` junto a `main.tex` (en Overleaf, súbelos al proyecto) y añade las
líneas `\input{...}` donde vayan los resultados. Los estilos `datos`, `puntoA`,
`puntoB`, `ajuste` y `teorica` ya están definidos en el preámbulo del informe.
