#!/usr/bin/env python3
"""
visualizar.py - Visualizaciones del Laboratorio 8 (Descubrimiento de Motifs).

Lee los CSV que genera el core en C++ (output/matriz_*.csv) y produce, para cada
conjunto de secuencias:
  - un heatmap de la matriz de frecuencias (posicion x {A,C,G,T}),
  - un grafico de barras del porcentaje de conservacion por posicion, y
  - un informe HTML autocontenido (output/reporte.html) con las imagenes
    embebidas en base64 (sin enlaces externos) y las tablas de resultados.

Uso:
    python3 visualizar.py
"""

import base64
import csv
import os

import matplotlib
matplotlib.use("Agg")  # backend sin ventana: solo escribe archivos PNG
import matplotlib.pyplot as plt
import numpy as np

# Carpeta donde el core deja los CSV y donde dejaremos las imagenes/HTML.
DIR_SALIDA = "output"

# Orden de bases usado en toda la matriz (coincide con el core en C++).
BASES = ["A", "C", "G", "T"]

# Conjuntos a visualizar: (id, titulo legible, ruta del CSV).
CONJUNTOS = [
    ("ejercicio1", "Ejercicio principal (10 secuencias)",
     os.path.join(DIR_SALIDA, "matriz_ejercicio1.csv")),
    ("anexo1", "Anexo 1 (8 secuencias)",
     os.path.join(DIR_SALIDA, "matriz_anexo1.csv")),
]


def leer_matriz(ruta):
    """Lee un CSV de matriz de frecuencias y devuelve una lista de filas (dict).

    Cada fila tiene: posicion (int), A/C/G/T (int), consenso (str),
    conservacion (float).
    """
    filas = []
    with open(ruta, newline="") as f:
        for fila in csv.DictReader(f):
            filas.append({
                "posicion": int(fila["posicion"]),
                "A": int(fila["A"]),
                "C": int(fila["C"]),
                "G": int(fila["G"]),
                "T": int(fila["T"]),
                "consenso": fila["consenso"],
                "conservacion": float(fila["conservacion"]),
            })
    return filas


def consenso_completo(filas):
    """Concatena los tokens de consenso por posicion (p. ej. TACGATG[ACT]C)."""
    return "".join(fila["consenso"] for fila in filas)


def heatmap_frecuencias(filas, titulo, ruta_png):
    """Genera el heatmap de la matriz de frecuencias y lo guarda como PNG."""
    posiciones = [fila["posicion"] for fila in filas]
    # Matriz 4 x L: filas = bases, columnas = posiciones.
    matriz = np.array([[fila[b] for fila in filas] for b in BASES])

    fig, ax = plt.subplots(figsize=(1.1 * len(posiciones) + 2, 3.2))
    im = ax.imshow(matriz, cmap="Blues", aspect="auto")

    # Etiquetas de ejes.
    ax.set_xticks(range(len(posiciones)))
    etiquetas_x = [f"{p}\n{fila['consenso']}" for p, fila in zip(posiciones, filas)]
    ax.set_xticklabels(etiquetas_x)
    ax.set_yticks(range(len(BASES)))
    ax.set_yticklabels(BASES)
    ax.set_xlabel("Posicion / consenso")
    ax.set_ylabel("Base")
    ax.set_title(f"Matriz de frecuencias - {titulo}")

    # Anotar el conteo en cada celda.
    maximo = matriz.max() if matriz.size else 1
    for i in range(matriz.shape[0]):
        for j in range(matriz.shape[1]):
            valor = matriz[i, j]
            color = "white" if valor > maximo / 2 else "black"
            ax.text(j, i, str(valor), ha="center", va="center", color=color)

    fig.colorbar(im, ax=ax, label="Frecuencia")
    fig.tight_layout()
    fig.savefig(ruta_png, dpi=120)
    plt.close(fig)


def grafico_conservacion(filas, titulo, ruta_png):
    """Genera el grafico de barras de conservacion por posicion (PNG)."""
    posiciones = [fila["posicion"] for fila in filas]
    valores = [fila["conservacion"] for fila in filas]
    # Verde para posiciones totalmente conservadas, naranja para variables.
    colores = ["#2ca02c" if v >= 100.0 else "#ff7f0e" for v in valores]

    fig, ax = plt.subplots(figsize=(1.1 * len(posiciones) + 2, 3.2))
    barras = ax.bar(range(len(posiciones)), valores, color=colores)
    ax.set_xticks(range(len(posiciones)))
    ax.set_xticklabels([f"{p}\n{fila['consenso']}"
                        for p, fila in zip(posiciones, filas)])
    ax.set_ylim(0, 109)
    ax.set_ylabel("Conservacion (%)")
    ax.set_xlabel("Posicion / consenso")
    ax.set_title(f"Conservacion por posicion - {titulo}")

    # Linea de referencia del 100% y etiquetas de valor.
    ax.axhline(100, color="gray", linestyle="--", linewidth=0.8)
    for rect, v in zip(barras, valores):
        ax.text(rect.get_x() + rect.get_width() / 2, v + 1,
                f"{v:.1f}", ha="center", va="bottom", fontsize=8)

    fig.tight_layout()
    fig.savefig(ruta_png, dpi=120)
    plt.close(fig)


def png_a_base64(ruta_png):
    """Codifica un PNG en base64 para embeberlo directamente en el HTML."""
    with open(ruta_png, "rb") as f:
        return base64.b64encode(f.read()).decode("ascii")


def tabla_html(filas):
    """Construye la tabla HTML de la matriz de frecuencias y conservacion."""
    cabecera = ("<tr><th>Pos</th><th>A</th><th>C</th><th>G</th><th>T</th>"
                "<th>Consenso</th><th>Conservacion</th></tr>")
    cuerpo = ""
    for fila in filas:
        resaltado = "" if fila["conservacion"] >= 100.0 else ' class="var"'
        cuerpo += (f"<tr{resaltado}><td>{fila['posicion']}</td>"
                   f"<td>{fila['A']}</td><td>{fila['C']}</td>"
                   f"<td>{fila['G']}</td><td>{fila['T']}</td>"
                   f"<td>{fila['consenso']}</td>"
                   f"<td>{fila['conservacion']:.1f}%</td></tr>")
    return f"<table>{cabecera}{cuerpo}</table>"


def main():
    secciones = []
    for id_conj, titulo, ruta_csv in CONJUNTOS:
        if not os.path.exists(ruta_csv):
            print(f"Aviso: no existe {ruta_csv}. Ejecuta primero ./motif")
            continue

        filas = leer_matriz(ruta_csv)
        ruta_heat = os.path.join(DIR_SALIDA, f"heatmap_{id_conj}.png")
        ruta_cons = os.path.join(DIR_SALIDA, f"conservacion_{id_conj}.png")
        heatmap_frecuencias(filas, titulo, ruta_heat)
        grafico_conservacion(filas, titulo, ruta_cons)

        consenso = consenso_completo(filas)
        global_pct = sum(f["conservacion"] for f in filas) / len(filas)
        n_conservadas = sum(1 for f in filas if f["conservacion"] >= 100.0)
        print(f"{titulo}: consenso={consenso}  global={global_pct:.1f}%  "
              f"posiciones 100%={n_conservadas}/{len(filas)}")

        secciones.append(f"""
        <section>
          <h2>{titulo}</h2>
          <p><b>Secuencia consenso:</b> <code>{consenso}</code> &nbsp;|&nbsp;
             <b>Longitud:</b> {len(filas)} &nbsp;|&nbsp;
             <b>Conservacion global:</b> {global_pct:.1f}% &nbsp;|&nbsp;
             <b>Posiciones 100%:</b> {n_conservadas}/{len(filas)}</p>
          <img src="data:image/png;base64,{png_a_base64(ruta_heat)}" alt="heatmap">
          <img src="data:image/png;base64,{png_a_base64(ruta_cons)}" alt="conservacion">
          {tabla_html(filas)}
        </section>""")

    html = f"""<!DOCTYPE html>
<html lang="es">
<head>
  <meta charset="utf-8">
  <title>Laboratorio 8 - Descubrimiento de Motifs</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 2rem; color: #222; }}
    h1 {{ color: #1f4e79; }}
    h2 {{ color: #2e6da4; border-bottom: 2px solid #ddd; padding-bottom: 4px; }}
    section {{ margin-bottom: 3rem; }}
    img {{ max-width: 100%; height: auto; display: block; margin: 1rem 0; }}
    table {{ border-collapse: collapse; margin-top: 1rem; }}
    th, td {{ border: 1px solid #ccc; padding: 4px 10px; text-align: center; }}
    th {{ background: #1f4e79; color: #fff; }}
    tr.var {{ background: #fff3e0; }}
    code {{ background: #f0f0f0; padding: 2px 6px; border-radius: 4px; }}
  </style>
</head>
<body>
  <h1>Laboratorio 8 - Descubrimiento de Motifs</h1>
  <p>Motif conservado identificado en cada conjunto de secuencias. Las filas
     resaltadas indican posiciones con variacion (mutaciones puntuales).</p>
  {''.join(secciones)}
</body>
</html>"""

    ruta_html = os.path.join(DIR_SALIDA, "reporte.html")
    with open(ruta_html, "w") as f:
        f.write(html)
    print(f"Informe HTML generado en: {ruta_html}")


if __name__ == "__main__":
    main()
