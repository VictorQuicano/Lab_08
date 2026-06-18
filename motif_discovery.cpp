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
// Este archivo se construye de forma incremental; por ahora cubre los pasos 1-2.
// Compilar con:  g++ -std=c++17 -O2 -Wall motif_discovery.cpp -o motif
// ============================================================================

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
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

// ----------------------------------------------------------------------------
// Paso 2 del laboratorio: identificacion de k-mers candidatos.
// ----------------------------------------------------------------------------

// Estadisticas de un k-mer dentro del conjunto de secuencias:
//   total     = numero de apariciones sumando todas las secuencias.
//   cobertura = numero de secuencias distintas en las que aparece (1..N).
struct EstadKmer {
    string kmer;
    int total = 0;
    int cobertura = 0;
};

// Cuenta todos los k-mers del conjunto con ventana deslizante (porta
// count_kmers del lab_04) y, a la vez, calcula en cuantas secuencias aparece
// cada uno. Devuelve un vector con la estadistica de cada k-mer distinto.
vector<EstadKmer> estadisticasKmers(const vector<Secuencia>& secuencias, int k) {
    map<string, int> total;
    map<string, int> cobertura;

    for (const Secuencia& s : secuencias) {
        if (k > static_cast<int>(s.bases.size())) {
            continue;  // la secuencia es mas corta que k
        }
        set<string> vistosEnEsta;  // k-mers distintos dentro de esta secuencia
        for (size_t i = 0; i + k <= s.bases.size(); ++i) {
            string km = s.bases.substr(i, k);
            total[km]++;
            vistosEnEsta.insert(km);
        }
        // Cada k-mer suma 1 a su cobertura por cada secuencia donde aparece.
        for (const string& km : vistosEnEsta) {
            cobertura[km]++;
        }
    }

    vector<EstadKmer> estad;
    estad.reserve(total.size());
    for (const auto& par : total) {
        estad.push_back(EstadKmer{par.first, par.second, cobertura[par.first]});
    }
    return estad;
}

// Ranking ingenuo: k-mers ordenados por frecuencia total (de mayor a menor).
// En secuencias con mucho ruido repetido (Anexo 1) este ranking lo dominan los
// repetidos de fondo, no el motif; por eso necesitamos tambien el discriminante.
vector<EstadKmer> rankingPorFrecuencia(vector<EstadKmer> estad) {
    sort(estad.begin(), estad.end(), [](const EstadKmer& a, const EstadKmer& b) {
        if (a.total != b.total) return a.total > b.total;
        return a.kmer < b.kmer;
    });
    return estad;
}

// Ranking discriminante para motifs. Un motif conservado tiende a aparecer
// ~una vez por secuencia en (casi) todas las secuencias, mientras que un
// repetido de fondo aparece muchas veces dentro de cada secuencia. Por eso:
//   1) nos quedamos solo con k-mers cuya media de apariciones por secuencia
//      donde aparecen (total/cobertura) sea <= 1.5, y
//   2) los ordenamos por cobertura (mas secuencias primero) y, a igual
//      cobertura, por menor total (mas "limpio", una vez por secuencia).
vector<EstadKmer> rankingDiscriminante(vector<EstadKmer> estad) {
    vector<EstadKmer> filtrados;
    for (const EstadKmer& e : estad) {
        double mediaPorSeq = static_cast<double>(e.total) / e.cobertura;
        if (mediaPorSeq <= 1.5) {
            filtrados.push_back(e);
        }
    }
    sort(filtrados.begin(), filtrados.end(),
         [](const EstadKmer& a, const EstadKmer& b) {
             if (a.cobertura != b.cobertura) return a.cobertura > b.cobertura;
             if (a.total != b.total) return a.total < b.total;
             return a.kmer < b.kmer;
         });
    return filtrados;
}

// Selecciona el k-mer candidato a motif: el mejor del ranking discriminante.
EstadKmer seleccionarCandidato(const vector<Secuencia>& secuencias, int k) {
    vector<EstadKmer> disc = rankingDiscriminante(estadisticasKmers(secuencias, k));
    if (disc.empty()) {
        return EstadKmer{};
    }
    return disc.front();
}

// Imprime, para un valor de k, los dos rankings y el candidato elegido.
void mostrarCandidatos(const vector<Secuencia>& secuencias, int k, int topN) {
    vector<EstadKmer> estad = estadisticasKmers(secuencias, k);
    vector<EstadKmer> porFrec = rankingPorFrecuencia(estad);
    vector<EstadKmer> disc = rankingDiscriminante(estad);
    int n = static_cast<int>(secuencias.size());

    cout << "  --- k = " << k << " ---\n";
    cout << "  Top por frecuencia total (ranking ingenuo):\n";
    for (int i = 0; i < topN && i < static_cast<int>(porFrec.size()); ++i) {
        cout << "    " << porFrec[i].kmer << "  total=" << porFrec[i].total
             << "  cobertura=" << porFrec[i].cobertura << "/" << n << "\n";
    }
    cout << "  Top discriminante (~1 vez por secuencia, mayor cobertura):\n";
    for (int i = 0; i < topN && i < static_cast<int>(disc.size()); ++i) {
        cout << "    " << disc[i].kmer << "  cobertura=" << disc[i].cobertura
             << "/" << n << "  total=" << disc[i].total << "\n";
    }
    if (!disc.empty()) {
        cout << "  -> Candidato a motif: " << disc.front().kmer
             << " (en " << disc.front().cobertura << "/" << n << " secuencias)\n";
    }
}

int main() {
    // Demo de los pasos 1-2: leer cada conjunto y listar k-mers candidatos.
    vector<string> rutas = {"data/ejercicio1.fasta", "data/anexo1.fasta"};
    vector<int> kValores = {7, 8, 9};

    for (const string& ruta : rutas) {
        vector<Secuencia> secuencias = leerFasta(ruta);
        cout << "==================================================\n";
        cout << "Archivo: " << ruta
             << "  (" << secuencias.size() << " secuencias)\n";
        cout << "==================================================\n";
        for (int k : kValores) {
            mostrarCandidatos(secuencias, k, 5);
            cout << "\n";
        }
    }

    return 0;
}
