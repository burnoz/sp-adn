#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <map>

using namespace std;

// Mapa de aminoacidos y sus codones
map<string, string> aminoacids = {
    {"Ala", "GCT, GCC, GCA, GCG"},
    {"Arg", "CGT, CGC, CGA, CGG; AGA, AGG"},
    {"Asn", "AAT, AAC"},
    {"Asp", "GAT, GAC"},
    {"Cys", "TGT, TGC"},
    {"Gln", "CAA, CAG"},
    {"Glu", "GAA, GAG"},
    {"Gly", "GGT, GGC, GGA, GGG"},
    {"His", "CAT, CAC"},
    {"Ile", "ATT, ATC, ATA"},
    {"Leu", "CTT, CTC, CTA, CTG; TTA, TTG"},
    {"Lys", "AAA, AAG"},
    {"Met", "ATG"},
    {"Phe", "TTT, TTC"},
    {"Pro", "CCT, CCC, CCA, CCG"},
    {"Ser", "TCT, TCC, TCA, TCG; AGT, AGC"},
    {"Thr", "ACT, ACC, ACA, ACG"},
    {"Trp", "TGG"},
    {"Tyr", "TAT, TAC"},
    {"Val", "GTT, GTC, GTA, GTG"},
    {"START", "ATG"},
    {"STOP", "TAA, TGA, TAG"}
};

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

    cout << "Header: " << header << endl;
    cout << "Sequence: " << sequence << endl;

    return sequence;
}

string getProtein(string sequence) {
    string protein = "";
    
    for(int i = 0; i < sequence.size(); i += 3){
        string codon = sequence.substr(i, 3);
        
        for(auto const& [aminoacid, codons] : aminoacids){
            if(kmp(codons, codon).size() > 0){
                if(aminoacid == "STOP"){
                    return protein;
                }
                
                protein += aminoacid + "-";
                break;
            }
        }
    }
    
    if(!protein.empty()){
        protein.pop_back();
    }
    
    return protein;
}

int main(){
    string sequence = getSequence("SARS-COV-2-MN908947.3.txt");
    string pattern = getSequence("gen-M.txt");

    vector<int> occurrences = kmp(sequence, pattern);

    cout << "Gen encontrado en las posiciones: ";
    
    for(int i = 0; i < occurrences.size(); i++){
        cout << occurrences[i] << " ";
    }

    string protein = getProtein(pattern);
    cout << "\nProteina traducida: " << protein << endl;
}
