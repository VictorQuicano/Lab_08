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

// ----------------------------------------------------------------------------
// Paso 3 del laboratorio: localizacion de las ocurrencias.
// ----------------------------------------------------------------------------

// Devuelve todas las posiciones (indice base 0) donde aparece 'patron' dentro
// de 'secuencia'. Una lista vacia significa que no aparece.
vector<int> todasPosiciones(const string& patron, const string& secuencia) {
    vector<int> posiciones;
    if (patron.empty()) return posiciones;
    size_t desde = 0;
    while (true) {
        size_t p = secuencia.find(patron, desde);
        if (p == string::npos) break;
        posiciones.push_back(static_cast<int>(p));
        desde = p + 1;  // permitir solapamientos
    }
    return posiciones;
}

// Devuelve la primera posicion de 'patron' en 'secuencia', o -1 si no aparece.
int primeraPosicion(const string& patron, const string& secuencia) {
    size_t p = secuencia.find(patron);
    return (p == string::npos) ? -1 : static_cast<int>(p);
}

// Selecciona el k-mer ancla (k = 7) que sirve de nucleo conservado para
// localizar el motif. Criterio: maxima cobertura (presente en mas secuencias)
// y, a igualdad de cobertura, el que aparece mas a la izquierda en promedio,
// porque ese nucleo marca el inicio del motif y la region se extiende hacia la
// derecha. En ambos conjuntos esto selecciona TACGATG.
string seleccionarAncla(const vector<Secuencia>& secuencias) {
    vector<EstadKmer> disc = rankingDiscriminante(estadisticasKmers(secuencias, 7));
    if (disc.empty()) return "";

    int maxCobertura = disc.front().cobertura;  // ya viene ordenado por cobertura
    string ancla;
    double mejorPosMedia = 1e18;

    for (const EstadKmer& e : disc) {
        if (e.cobertura != maxCobertura) continue;
        // Posicion media de la primera ocurrencia en las secuencias donde aparece.
        long sumaPos = 0;
        int conteo = 0;
        for (const Secuencia& s : secuencias) {
            int p = primeraPosicion(e.kmer, s.bases);
            if (p >= 0) { sumaPos += p; ++conteo; }
        }
        double media = (conteo > 0) ? static_cast<double>(sumaPos) / conteo : 1e18;
        if (media < mejorPosMedia) {
            mejorPosMedia = media;
            ancla = e.kmer;
        }
    }
    return ancla;
}

// Imprime, para el ancla elegida, su posicion en cada secuencia y comprueba si
// las ocurrencias caen en regiones similares.
void mostrarOcurrencias(const vector<Secuencia>& secuencias, const string& ancla) {
    cout << "  Ancla conservada (nucleo del motif): " << ancla << "\n";
    cout << "  Posicion (base 0) de la ocurrencia en cada secuencia:\n";
    int presentes = 0, minPos = 1e9, maxPos = -1;
    for (const Secuencia& s : secuencias) {
        vector<int> pos = todasPosiciones(ancla, s.bases);
        cout << "    " << s.nombre << ": ";
        if (pos.empty()) {
            cout << "(no aparece)";
        } else {
            for (size_t i = 0; i < pos.size(); ++i) {
                cout << pos[i] << (i + 1 < pos.size() ? ", " : "");
            }
            ++presentes;
            minPos = min(minPos, pos.front());
            maxPos = max(maxPos, pos.front());
        }
        cout << "\n";
    }
    cout << "  Aparece en " << presentes << "/" << secuencias.size()
         << " secuencias. Rango de la primera posicion: [" << minPos << ", "
         << maxPos << "]";
    if (maxPos - minPos <= 5) {
        cout << " -> ocurrencias en regiones muy similares.\n";
    } else {
        cout << " -> posiciones dispersas (motif embebido en distinto contexto).\n";
    }
}

// ----------------------------------------------------------------------------
// Pasos 4 y 5 del laboratorio: extraccion de la region conservada y
// alineamiento multiple.
// ----------------------------------------------------------------------------

// Una region del motif extraida de una secuencia concreta.
struct Region {
    string nombre;  // secuencia de origen
    int posicion;   // posicion (base 0) donde empieza la region
    string bases;   // subsecuencia extraida (longitud = longitud del motif)
};

// Determina la longitud del motif extendiendo el ancla hacia la derecha. Se
// agrega una base mas mientras exista un k-mer que tenga el ancla como prefijo
// y siga apareciendo en al menos 2 secuencias (sigue conservado). Se limita a
// k = 9 porque el laboratorio prueba k = 7, 8 y 9. En ambos conjuntos da 9.
int determinarLongitudMotif(const vector<Secuencia>& secuencias,
                            const string& ancla, int maxK = 9) {
    int L = static_cast<int>(ancla.size());
    for (int k = L + 1; k <= maxK; ++k) {
        vector<EstadKmer> estad = estadisticasKmers(secuencias, k);
        int mejorCobertura = 0;
        for (const EstadKmer& e : estad) {
            // El k-mer debe empezar exactamente con el ancla (extension a la derecha).
            if (e.kmer.compare(0, ancla.size(), ancla) == 0) {
                mejorCobertura = max(mejorCobertura, e.cobertura);
            }
        }
        if (mejorCobertura >= 2) {
            L = k;  // el prefijo conservado sigue presente, extendemos
        } else {
            break;  // ya no esta conservado, detenemos la extension
        }
    }
    return L;
}

// Extrae de cada secuencia la region del motif: ubica el ancla y toma una
// ventana de longitud 'longitud' a partir de ahi. Se conservan las variantes
// (no se corrigen las mutaciones puntuales). Las secuencias donde el ancla no
// aparece, o donde la ventana se saldria del limite, se omiten.
vector<Region> extraerRegiones(const vector<Secuencia>& secuencias,
                               const string& ancla, int longitud) {
    vector<Region> regiones;
    for (const Secuencia& s : secuencias) {
        int p = primeraPosicion(ancla, s.bases);
        if (p < 0) continue;
        if (p + longitud > static_cast<int>(s.bases.size())) continue;
        regiones.push_back(Region{s.nombre, p, s.bases.substr(p, longitud)});
    }
    return regiones;
}

// Imprime las regiones extraidas apiladas (alineamiento multiple por columnas)
// con una linea de marcas: '*' si la columna esta totalmente conservada y ' '
// si presenta variacion entre las secuencias.
void mostrarAlineamiento(const vector<Region>& regiones) {
    if (regiones.empty()) return;
    int L = static_cast<int>(regiones.front().bases.size());

    for (const Region& r : regiones) {
        cout << "    " << r.nombre << "\t" << r.bases << "\n";
    }
    // Linea de conservacion por columna.
    cout << "    conserv.\t";
    for (int col = 0; col < L; ++col) {
        char base = regiones.front().bases[col];
        bool conservada = true;
        for (const Region& r : regiones) {
            if (r.bases[col] != base) { conservada = false; break; }
        }
        cout << (conservada ? '*' : ' ');
    }
    cout << "\n";
}

// ----------------------------------------------------------------------------
// Pasos 6 y 7 del laboratorio: matriz de frecuencias y secuencia consenso.
// ----------------------------------------------------------------------------

// Matriz de frecuencias por posicion. conteo[pos] guarda el numero de
// apariciones de A, C, G y T (en ese orden) en esa columna del alineamiento.
struct MatrizFrecuencias {
    int L = 0;                      // longitud del motif (numero de columnas)
    int numSecuencias = 0;          // numero de regiones alineadas (filas)
    vector<vector<int>> conteo;     // conteo[L][4], orden de bases = "ACGT"
};

// Indice de una base dentro del orden "ACGT" (0..3), o -1 si no es valida.
int indiceBase(char c) {
    size_t p = BASES.find(c);
    return (p == string::npos) ? -1 : static_cast<int>(p);
}

// Construye la matriz de frecuencias contando, por cada posicion del
// alineamiento, cuantas veces aparece cada nucleotido.
MatrizFrecuencias matrizFrecuencias(const vector<Region>& regiones) {
    MatrizFrecuencias m;
    if (regiones.empty()) return m;
    m.L = static_cast<int>(regiones.front().bases.size());
    m.numSecuencias = static_cast<int>(regiones.size());
    m.conteo.assign(m.L, vector<int>(4, 0));

    for (const Region& r : regiones) {
        for (int col = 0; col < m.L; ++col) {
            int idx = indiceBase(r.bases[col]);
            if (idx >= 0) m.conteo[col][idx]++;
        }
    }
    return m;
}

// Imprime la matriz de frecuencias como una tabla posicion x {A,C,G,T}.
void mostrarMatriz(const MatrizFrecuencias& m) {
    cout << "    Pos\tA\tC\tG\tT\n";
    for (int col = 0; col < m.L; ++col) {
        cout << "    " << (col + 1);
        for (int b = 0; b < 4; ++b) cout << "\t" << m.conteo[col][b];
        cout << "\n";
    }
}

// Obtiene la secuencia consenso. Por cada posicion se elige la base mas
// frecuente. Cuando la columna presenta 3 o mas bases distintas (variabilidad
// real) se usa notacion degenerada entre corchetes con todas las bases
// observadas, p. ej. [ACT]; con una sola base se reporta esa base y con dos
// (mutacion puntual) se reporta la mas frecuente. Esta convencion reproduce el
// consenso TACGATG[ACT]C del marco teorico.
string consenso(const MatrizFrecuencias& m) {
    string resultado;
    for (int col = 0; col < m.L; ++col) {
        int distintas = 0, mejorBase = 0, mejorConteo = -1;
        for (int b = 0; b < 4; ++b) {
            if (m.conteo[col][b] > 0) distintas++;
            if (m.conteo[col][b] > mejorConteo) {
                mejorConteo = m.conteo[col][b];
                mejorBase = b;
            }
        }
        if (distintas >= 3) {
            string grupo = "[";
            for (int b = 0; b < 4; ++b) {
                if (m.conteo[col][b] > 0) grupo += BASES[b];
            }
            grupo += "]";
            resultado += grupo;
        } else {
            resultado += BASES[mejorBase];
        }
    }
    return resultado;
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

        string ancla = seleccionarAncla(secuencias);
        cout << "  [Paso 3] Localizacion de ocurrencias\n";
        mostrarOcurrencias(secuencias, ancla);
        cout << "\n";

        int L = determinarLongitudMotif(secuencias, ancla);
        vector<Region> regiones = extraerRegiones(secuencias, ancla, L);
        cout << "  [Paso 4] Extraccion de la region conservada (longitud "
             << L << ")\n";
        for (const Region& r : regiones) {
            cout << "    " << r.nombre << " (pos " << r.posicion << "): "
                 << r.bases << "\n";
        }
        cout << "  [Paso 5] Alineamiento multiple de las regiones\n";
        mostrarAlineamiento(regiones);
        cout << "\n";

        MatrizFrecuencias m = matrizFrecuencias(regiones);
        cout << "  [Paso 6] Matriz de frecuencias\n";
        mostrarMatriz(m);
        cout << "  [Paso 7] Secuencia consenso: " << consenso(m) << "\n";
        cout << "\n";
    }

    return 0;
}
