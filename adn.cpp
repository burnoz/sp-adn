#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <map>

using namespace std;

// Mapa de aminoacidos y sus codones
// map<string, string> aminoacids = {
//     {"Ala", "GCT, GCC, GCA, GCG"},
//     {"Arg", "CGT, CGC, CGA, CGG; AGA, AGG"},
//     {"Asn", "AAT, AAC"},
//     {"Asp", "GAT, GAC"},
//     {"Cys", "TGT, TGC"},
//     {"Gln", "CAA, CAG"},
//     {"Glu", "GAA, GAG"},
//     {"Gly", "GGT, GGC, GGA, GGG"},
//     {"His", "CAT, CAC"},
//     {"Ile", "ATT, ATC, ATA"},
//     {"Leu", "CTT, CTC, CTA, CTG; TTA, TTG"},
//     {"Lys", "AAA, AAG"},
//     {"Met", "ATG"},
//     {"Phe", "TTT, TTC"},
//     {"Pro", "CCT, CCC, CCA, CCG"},
//     {"Ser", "TCT, TCC, TCA, TCG; AGT, AGC"},
//     {"Thr", "ACT, ACC, ACA, ACG"},
//     {"Trp", "TGG"},
//     {"Tyr", "TAT, TAC"},
//     {"Val", "GTT, GTC, GTA, GTG"},
//     {"START", "ATG"},
//     {"STOP", "TAA, TGA, TAG"}
// };

map<string, string> aminoacids = {
    {"A", "GCT,GCC,GCA,GCG"},
    {"R", "CGT,CGC,CGA,CGG,AGA,AGG"},
    {"N", "AAT,AAC"},
    {"D", "GAT,GAC"},
    {"C", "TGT,TGC"},
    {"Q", "CAA,CAG"},
    {"E", "GAA,GAG"},
    {"G", "GGT,GGC,GGA,GGG"},
    {"H", "CAT,CAC"},
    {"I", "ATT,ATC,ATA"},
    {"L", "CTT,CTC,CTA,CTG,TTA,TTG"},
    {"K", "AAA,AAG"},
    {"M", "ATG"},
    {"F", "TTT,TTC"},
    {"P", "CCT,CCC,CCA,CCG"},
    {"S", "TCT,TCC,TCA,TCG,AGT,AGC"},
    {"T", "ACT,ACC,ACA,ACG"},
    {"W", "TGG"},
    {"Y", "TAT,TAC"},
    {"V", "GTT,GTC,GTA,GTG"},
    {"X", "TAA,TGA,TAG"} // STOP
};

// Algoritmo de Manacher
// Complejidad: O(n)
vector<int> manacher(string palabra){
    // Texto modificado: @t#ex#t#o#$
    string modified = "@";

    for(int i = 0; i < palabra.size(); i++){
        modified += "#" + palabra.substr(i, 1);
    }

    modified += "#$";

    int centro = 0;
    int limite = 0;
    int gap = 0;
    vector<int> P(modified.size(), 0);

    for(int i = 1; i < modified.size() - 1; i++){
        // Si i esta dentro del limite, usa la simetria
        if(i < limite){
            int simetrico = 2 * centro - i;
            P[i] = min(limite - i, P[simetrico]);
        }

        // Intenta expandir el palindromo alrededor de i
        gap = P[i] + 1;

        // Mientras los caracteres a ambos lados coincidan, expande
        while(modified[i + gap] == modified[i - gap]){
            P[i]++;
            gap++;
        }

        // Si el palindromo se extiende mas alla del limite, actualiza centro y limite
        if(i + P[i] > limite){
            limite = i + P[i];
            centro = i;
        }
    }

    return P;
}

// Obtiene el prefijo mas largo que es sufijo
// Complejidad: O(n)
vector<int> getLPS(string text){
    int n = text.size();

    vector<int> lps(n, 0);

    int len_prefix = 0; // Longitud del prefijo mas largo que es sufijo
    int i = 1;

    while(i < n){
        // Si hay coincidencia, extiende el prefijo
        if(text[i] == text[len_prefix]){
            len_prefix++;
            lps[i] = len_prefix;
            i++;
        }

        else{
            // Si no hay coincidencia y len_prefix no es 0, retrocede
            if(len_prefix != 0){
                len_prefix = lps[len_prefix - 1];
            }

            // Si len_prefix es 0, no hay prefijo, avanza
            else{
                lps[i] = 0;
                i++;
            }
        }
    }

    return lps;
}

// KMP
// Complejidad: O(n + m), donde n es el tamano del texto y m el del patron
vector<int> kmp(string text, string pattern){
    int n = text.size();
    int m = pattern.size();

    vector<int> lps = getLPS(pattern);
    vector<int> occurrences;

    int len_prefix = 0;

    int i = 0;
    int j = 0;

    while(i < n){
        // Si hay coincidencia, avanza en ambos
        if(text[i] == pattern[j]){
            i++;
            j++;
        }

        // Si se encontro el patron, guarda la posicion
        if(j == m){
            occurrences.push_back(i - j);
            j = lps[j - 1];
        }

        // Si no hay coincidencia
        else if(i < n && text[i] != pattern[j]){
            // Si j no es 0, retrocede en el patron
            if(j > 0){
                j = lps[j - 1];
            }

            // Si j es 0, avanza en el texto
            else{
                i++;
            }
        }
    }

    return occurrences;
}

string getSequence(string filename) {
    ifstream file(filename);
    string line;
    string sequence; 
    string header;
    
    bool first = true;
    
    while(getline(file, line)){
        if(first){
            first = false;
            header = line;
            continue; 
        }

        sequence += line;
    }

    // cout << "Header: " << header << endl;
    // cout << "Sequence: " << sequence << endl;

    return sequence;
}


vector<pair<string, string>> getProteinsFromFile(string filename){
    ifstream file(filename);
    string line;

    vector<pair<string, string>> proteins; // Par de nombre y secuencia

    string currentName = "";
    string currentSequence = "";

    while(getline(file, line)){
        if(line.empty()){
            continue; // Salta lineas vacias
        } 

        if(line[0] == '>'){
            // Si ya hay una secuencia en proceso, guardala
            if(!currentName.empty()){
                proteins.push_back({currentName, currentSequence});
                currentSequence = "";
            }

            currentName = line.substr(1); // Quita el '>'
        }

        else{
            currentSequence += line; // Agrega a la secuencia actual
        }
    }

    // Agrega la ultima secuencia si existe
    if(!currentName.empty()){
        proteins.push_back({currentName, currentSequence});
    }

    return proteins;
}

string getProteins(string sequence) {
    string protein = "";
    
    for(int i = 0; i < sequence.size(); i += 3){
        string codon = sequence.substr(i, 3);
        
        for(auto [aminoacid, codons] : aminoacids){
            if(kmp(codons, codon).size() > 0){
                protein += aminoacid;
                break;
            }
        }
    }
    
    return protein;
}

// Traduce un codon a su aminoacido
char codonToAminoacid(string codon) {
    for(auto [aminoacid, codons] : aminoacids) {
        if(kmp(codons, codon).size() > 0) {
            return aminoacid[0];
        }
    }

    return '?';
}

int main(){
    // Archivos
    string genMFile = "gen-M";
    string genSFile = "gen-S";
    string genORF1ABFile = "gen-ORF1AB";
    string genomeSeqFile = "SARS-COV-2-MN908947.3";

    vector<string> genFiles = {genMFile, genSFile, genORF1ABFile};

    vector<string> genSequences;
    for(int i = 0; i < genFiles.size(); i++){
        genSequences.push_back(getSequence(genFiles[i] + ".txt"));
    }

    string genomeSequence = getSequence(genomeSeqFile + ".txt");

    // Parte 1, indices de aparición de cada uno de los tres genes
    for(int i = 0; i < genSequences.size(); i++){
        vector<int> occurrences = kmp(genomeSequence, genSequences[i]);
        cout << genFiles[i] << " ";

        if(occurrences.size() > 0){
            cout << "encontrado en la posicion: " << occurrences[0] << endl;
        }

        else{
            cout << "no encontrado" << endl;
        }

        // Primeros 12 caracteres del gen
        cout << "Primeros 12 caracteres del gen: " << genSequences[i].substr(0, 12) << endl;
        cout << endl;
    }

    // Parte 2, palíndromo mas largo en cada uno de los tres genes
    for(int i = 0; i < genSequences.size(); i++){
        vector<int> P = manacher(genSequences[i]);

        int max_index = 0;
        int max_length = 0;

        for(int j = 0; j < P.size(); j++){
            if(P[j] > max_length){
                max_length = P[j];
                max_index = j;
            }
        }

        int inicio = (max_index - max_length) / 2;

        cout << "Palindromo mas largo en " << genFiles[i] << " desde " << inicio << " hasta " << inicio + max_length - 1 << ":\n";
        cout << "Longitud: " << max_length << endl;
        cout << "Palindromo: " << genSequences[i].substr(inicio, max_length) << endl;
        cout << endl;

        // Guarda el palindromo en un archivo
        ofstream palFile;
        palFile.open(genFiles[i] + "-palindrome.txt");
        palFile << genSequences[i].substr(inicio, max_length);
        palFile.close();
    }

    // Parte 3, dado el archivo de las proteinas, el cual contiene sus secuencias de aminoacidos, encuentra en cuales secciones del virus se produce cada proteina
    vector<pair<string, string>> proteins = getProteinsFromFile("seq-proteins.txt");
    vector<int> slipperyPos;

    // Convierte considerando los 3 reading frames y los guarda en un vector y archivo
    vector<string> genomeProteins(3);
    for(int frame = 0; frame < 3; frame++){
        string frameSeq = genomeSequence.substr(frame);
        slipperyPos = kmp(genomeSequence.substr(frame), "TTTAAAC");

        string proteinsInFrame = getProteins(frameSeq);
        genomeProteins[frame] = proteinsInFrame;

        // Guarda en archivo
        ofstream frameFile;
        frameFile.open("genome-proteins-frame-" + to_string(frame) + ".txt");
        frameFile << proteinsInFrame;
        frameFile.close();
    }

    for(int i = 0; i < proteins.size(); i++){
        string proteinName = proteins[i].first;
        string proteinSeq = proteins[i].second;

        bool found = false;

        for(int frame = 0; frame < 3; frame++){
            vector<int> occurrences = kmp(genomeProteins[frame], proteinSeq);

            if(occurrences.size() > 0){
                found = true;
                cout << "Proteina " << proteinName << " encontrada en reading frame " << frame << " en las posiciones: ";
                
                for(int j = 0; j < occurrences.size(); j++){
                    cout << occurrences[j] * 3 << " ";
                }

                cout << endl;
            }
        }

        if(!found){
            cout << "Proteina " << proteinName << " no encontrada en ningun reading frame" << endl;
        }
    }

    // Parte 4, Compara las versiones del genoma del virus de Wuhan, 2019 vs Texas, 2020. Determina donde  difieren, y si tales diferencias resultan en aminoácidos diferentes
    // Comparación Wuhan vs Texas
    string wuhan = getSequence("SARS-COV-2-MN908947.3.txt");
    string texas = getSequence("SARS-COV-2-MT106054.1.txt");
    int minLen = min(wuhan.size(), texas.size());
    
    cout << "Diferencias:" << endl;
    
    for(int i = 0; i + 2 < minLen; i += 3) {
        string codonWuhan = wuhan.substr(i, 3);
        string codonTexas = texas.substr(i, 3);

        if(codonWuhan != codonTexas) {
            char aminoWuhan = codonToAminoacid(codonWuhan);
            char aminoTexas = codonToAminoacid(codonTexas);
            
            cout << " Indice: " << i << " | Wuhan: " << codonWuhan << " (" << aminoWuhan << ")" << " | Texas: " << codonTexas << " (" << aminoTexas << ")";
            
            if(aminoWuhan != aminoTexas){
                cout << " <-- Cambio de aminoacido";
            }

            cout << endl;
        }
    }
}
