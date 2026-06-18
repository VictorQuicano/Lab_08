# Laboratorio 8 — Descubrimiento de Motifs

Curso de Bioinformatica (EPCC - UNSA). Analisis de secuencias de ADN para
identificar un motif conservado: seleccion de k-mers candidatos, localizacion de
ocurrencias, extraccion de la region comun, alineamiento multiple, matriz de
frecuencias, secuencia consenso, evaluacion de conservacion y reporte del motif.

## Estructura

```
lab_08/
├── motif_discovery.cpp   # core: pipeline completo de descubrimiento de motifs (C++17)
├── visualizar.py         # heatmap, grafico de conservacion e informe HTML
├── informe_lab08.md      # informe del laboratorio (entregable)
├── data/
│   ├── ejercicio1.fasta  # 10 secuencias del ejercicio principal
│   └── anexo1.fasta      # 8 secuencias del Anexo 1
└── output/               # resultados generados (CSV, PNG, HTML)
```

## Uso

```bash
# 1. Compilar el core
g++ -std=c++17 -O2 -Wall motif_discovery.cpp -o motif

# 2. Ejecutar indicando el archivo FASTA a analizar (un archivo por ejecucion).
#    El CSV de salida se nombra a partir del archivo de entrada
#    (output/matriz_<nombre>.csv).
./motif data/ejercicio1.fasta
./motif data/anexo1.fasta

# 3. Generar las visualizaciones y el informe HTML a partir de los CSV
python3 visualizar.py

# 4. Abrir output/reporte.html en el navegador
```

El motif conservado esperado es `TACGATG[ACT]C` (longitud 9): las primeras siete
posiciones (`TACGATG`) estan totalmente conservadas y sirven de ancla; las dos
ultimas presentan variaciones puntuales entre secuencias.
