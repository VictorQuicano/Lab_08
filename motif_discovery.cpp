// ============================================================================
// Laboratorio 8 - Descubrimiento de Motifs
// Curso de Bioinformatica (EPCC - UNSA)
//
// Core del laboratorio. Implementa el pipeline completo de descubrimiento de
// motifs sobre conjuntos de secuencias de ADN:
//   1. Lectura de secuencias en formato FASTA.
//   2. Seleccion de k-mers candidatos (k = 7, 8, 9).
//   3. Localizacion de las ocurrencias de cada candidato.
//   4. Extraccion de la region conservada.
//   5. Alineamiento multiple de las regiones extraidas.
//   6. Construccion de la matriz de frecuencias.
//   7. Obtencion de la secuencia consenso.
//   8. Evaluacion del grado de conservacion.
//   9. Reporte final del motif.
//
// Este archivo se construye de forma incremental; por ahora cubre el paso 1.
// Compilar con:  g++ -std=c++17 -O2 -Wall motif_discovery.cpp -o motif
// ============================================================================

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

// Bases validas del alfabeto de ADN estandar.
const string BASES = "ACGT";

// Una secuencia leida de un archivo FASTA: su cabecera y sus bases.
struct Secuencia {
    string nombre;  // cabecera sin el caracter '>'
    string bases;   // secuencia en mayusculas, solo {A,C,G,T}
};

// Devuelve true si la base pertenece al alfabeto de ADN {A,C,G,T}.
bool esBaseValida(char c) {
    return BASES.find(c) != string::npos;
}

// Elimina de la secuencia cualquier caracter que no sea {A,C,G,T} y la pasa a
// mayusculas. Equivale a filter_valid_bases del lab_04: asi los separadores
// 'N' o caracteres ambiguos no generan k-mers artificiales mas adelante.
string filtrarBases(const string& seq) {
    string limpia;
    limpia.reserve(seq.size());
    for (char c : seq) {
        char mayus = toupper(static_cast<unsigned char>(c));
        if (esBaseValida(mayus)) {
            limpia.push_back(mayus);
        }
    }
    return limpia;
}

// Lee un archivo FASTA y devuelve todas las secuencias en orden de aparicion.
// Las lineas que empiezan con '>' son cabeceras; las lineas siguientes hasta
// la proxima '>' se concatenan como la secuencia de ese registro. Es el
// equivalente en C++ de parse_fasta del lab_04.
vector<Secuencia> leerFasta(const string& ruta) {
    ifstream archivo(ruta);
    if (!archivo.is_open()) {
        cerr << "Error: no se pudo abrir el archivo " << ruta << "\n";
        return {};
    }

    vector<Secuencia> secuencias;
    string linea;
    Secuencia actual;
    bool hayActual = false;  // indica si ya empezamos a leer un registro

    while (getline(archivo, linea)) {
        // Quitar posibles retornos de carro de archivos con saltos de Windows.
        if (!linea.empty() && linea.back() == '\r') {
            linea.pop_back();
        }
        if (linea.empty()) {
            continue;  // ignorar lineas en blanco
        }

        if (linea[0] == '>') {
            // Guardar el registro anterior antes de empezar uno nuevo.
            if (hayActual) {
                secuencias.push_back(actual);
            }
            actual = Secuencia{};
            actual.nombre = linea.substr(1);  // quitar el '>'
            hayActual = true;
        } else {
            // Acumular las bases del registro actual (ya filtradas).
            actual.bases += filtrarBases(linea);
        }
    }
    // Guardar el ultimo registro leido.
    if (hayActual) {
        secuencias.push_back(actual);
    }

    return secuencias;
}

int main() {
    // Demo del paso 1: leer ambos conjuntos e imprimir un resumen.
    vector<string> rutas = {"data/ejercicio1.fasta", "data/anexo1.fasta"};

    for (const string& ruta : rutas) {
        vector<Secuencia> secuencias = leerFasta(ruta);
        cout << "Archivo: " << ruta << "\n";
        cout << "Secuencias leidas: " << secuencias.size() << "\n";
        for (const Secuencia& s : secuencias) {
            cout << "  " << s.nombre << " (len " << s.bases.size() << "): "
                 << s.bases << "\n";
        }
        cout << "\n";
    }

    return 0;
}
