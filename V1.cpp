#include <complex>
#include <iterator>
#include <string>
#include <vector>
#include <iostream>
#include <random>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <iostream>
#include <fstream>

namespace fs = std::filesystem;
using namespace std;
constexpr double PI = 3.14159265359; 

typedef struct{ 
  float x; 
  float y; 
  float z;
} vec3;

float norm(vec3 v){
  return sqrt(pow(v.x, 2) + pow(v.y, 2) + pow(v.z, 2));
}

void neutronvect(int n, vector<vec3> &positions, vector<vec3> &vitesses , vector<float> &temps , vec3 NouvellePosition){
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<float> dis(-1.0f, 1.0f);

  for (int i = 0; i < n; i++){
    positions.push_back(NouvellePosition);
    vec3 v;

    v.x = dis(gen);
    v.y = dis(gen);
    v.z = dis(gen);
    float Vnorm = norm(v);
    v.x = v.x/Vnorm * 10; 
    v.y = v.y/Vnorm * 10;
    v.z = v.z/Vnorm * 10;
    vitesses.push_back(v);
    temps.push_back(0.0);
  }
}

void deplacement(vector<vec3> &positions, vector<vec3> &vitesses, vector<float> &temps, float t){
  for (int i = 0; i < positions.size(); i++) {
    positions[i].x += vitesses[i].x * t;
    positions[i].y += vitesses[i].y * t;
    positions[i].z += vitesses[i].z * t;
    temps[i] += t;
  }
}
void moderation(vector<vec3> &vitesses, float t){
  for (int i = 0; i < vitesses.size(); i++) {
    vitesses[i].x *= pow(0.9, t); 
    vitesses[i].y *= pow(0.9, t); 
    vitesses[i].z *= pow(0.9, t); 
  }
}
void combustible(vector<vec3> &comb, float rayon, float hauteur, float nbNoyaudtau, float enrichissement, int n){
  std::random_device rd;
  std::mt19937 gen(rd());
  vec3 v;
  v.x = 0;
  v.y = 0;
  v.z = hauteur/2;
  comb.push_back(v);
  int k = 1;
  while (k < n) {
    int i = 0;
    float angle = 2 * PI / (nbNoyaudtau * k);
    while (i < nbNoyaudtau * k +1) {
      float x = cos(angle*i) * rayon * k; 
      float y = sin(angle*i) * rayon * k;
      if(i == nbNoyaudtau && k == 1){
        x = 0;
        y = 0;
      }
      float j = 0;
      while (j < hauteur){
        int m = 1;
        while (m < enrichissement) {
          vec3 nv;
          std::normal_distribution<double> normalDistx(x, 1);
          nv.x = normalDistx(gen);
          std::normal_distribution<double> normalDisty(y, 1);
          nv.y = normalDisty(gen);
          nv.z = j;
          comb.push_back(nv);
        m++;
        }
      j+=0.1;
      }
    i++;
    }
  k++;
  }
}
void distance(vector<float> &dist, vector<vec3> &coordonne, vec3 point){
    for (int j = 0; j < coordonne.size(); j++ ) {
      float s = 0;
      s += pow(coordonne[j].x - point.x, 2.f);
      s += pow(coordonne[j].y - point.y, 2.f);
      s += pow(coordonne[j].z - point.z, 2.f);
      dist.push_back(sqrt(s));
    }
}


void DesintegrationAux(vector <vec3> &positions, vector<float> &temps){
  int dmin = 3;
  int bulle_max = 30;
  for (int i = 0; i < positions.size(); i++) {
    int counter = 0; 
    for (int j = 0; j < positions.size(); j++ ) {
      float s = 0;
      s += pow(positions[j].x - positions[i].x, 2.f);
      s += pow(positions[j].y - positions[i].y, 2.f);
      s += pow(positions[j].z - positions[i].z, 2.f);
      s = sqrt(s);
      if (s > 0 && s <= dmin){
        counter++;
        if (counter > bulle_max){
          temps[j] = std::numeric_limits<float>::max();
        }
      }
    }
  }
}

void Desintegration(vector<vec3> &positionsNew, vector<vec3> &vitessesNew, vector<float> &tempsNew, vector<vec3> &positions, vector<vec3> &vitesses, vector<float> &temps, float dinf, float hauteur){
    float AgeFermi = 2.f;
    
    DesintegrationAux(positions, temps);

    

    for (int i = 0; i < positions.size(); i++) {
        float dist = pow(positions[i].x, 2.f) +
                     pow(positions[i].y, 2.f) +
                     pow(positions[i].z - hauteur / 2, 2.f);

        if (temps[i] < AgeFermi && sqrt(dist) < dinf) {
            positionsNew.push_back(positions[i]);
            vitessesNew.push_back(vitesses[i]);
            tempsNew.push_back(temps[i]);
        }
    }
    if(positionsNew.size() == 0) {
        std::cerr << "Erreur : Plus de neutrons" << std::endl;
    }

}


void fission(vector<vec3> &positions, vector<vec3> &vitesses, vector<float> &temps, vector<vec3> &comb, float rayon){
  std::random_device rd;
  std::mt19937 gen(rd());
  float v0 = 1.8f;
  float m = 0.2;
  float max = positions.size();
  for (int i = 0; i < max; i++) {
    vector<float> nearbyindices;
    
    float noyau = 0;
    float min = std::numeric_limits<float>::max();

    for (int j = 0; j < comb.size(); j++) {
      float s = 0;
      s += pow(comb[j].x - positions[i].x, 2.f);
      s += pow(comb[j].y - positions[i].y, 2.f);
      s += pow(comb[j].z - positions[i].z, 2.f);
      s = sqrt(s);
      if(s < rayon){
        nearbyindices.push_back(s);
      }
      if(s < min){
        noyau = j;
        min = s;
      }
    }
    if (nearbyindices.size() != 0) {
    
      float Ec = pow(vitesses[i].x, 2.0f) +
                 pow(vitesses[i].y, 2.0f) +
                 pow(vitesses[i].z, 2.0f);
      Ec = Ec * m / 2;
      float R;
      if (Ec < v0){
        R = 0;
      } else {
        R = pow((1 - sqrt(Ec-v0))/(1 + sqrt(Ec-v0)), 2.f);
      }
      std::normal_distribution<double> normalDist1(0, 1);
      if(normalDist1(gen) > R) {
        neutronvect(3, positions, vitesses, temps, comb[noyau]);
        temps[i] = std::numeric_limits<float>::max();
      } else{
        std::normal_distribution<double> normalDist2(-1, 1);
        vitesses[i].x *= normalDist2(gen);
        vitesses[i].y *= normalDist2(gen);
        vitesses[i].z *= normalDist2(gen);
        temps[i] = 0.f;
      }
    }
  }
}

void ecrireDansFichier(const std::string& dossier, const std::string& name, std::vector<vec3>& vect) {
    // Définir le chemin du fichier
    std::string cheminFichier = "./" + dossier + "/" + name + ".csv";

    // Créer le dossier s'il n'existe pas
    if (!fs::exists(dossier)) {
        fs::create_directory(dossier);
    }

    // Ouvrir le fichier en mode écriture
    std::ofstream fichier(cheminFichier);

    // Vérifier si le fichier est ouvert avec succès
    if (fichier.is_open()) {
        // Écrire les en-têtes
        fichier << "x,y,z" << std::endl;

        // Écrire les données dans le fichier
        for (const auto& ligne : vect) {
            // Écrire chaque composant de vec3 séparément
            fichier << ligne.x << "," << ligne.y << "," << ligne.z << std::endl;
        }

        // Fermer le fichier
        fichier.close();
    } else {
        std::cerr << "Erreur : Impossible d'ouvrir le fichier." << std::endl;
    }
}

void cycle(int n, float t, float rayon, vector<vec3> &positions, vector<vec3> &vitesses, vector<float> &temps, vector<vec3> &comb, float dinf, float hauteur, string dossier){
  ecrireDansFichier(dossier, "comb", comb);
  vector<vec3> positionsNew;
  vector<vec3> vitessesNew;
  vector<float> tempsNew;
  float sommeK = 1;
  float tmp = positions.size();
  for (int i = 0; i < n; i++) {
    cout << i << "/" << n << "->" << positions.size() << endl;
    cout << "keff = " << sommeK/(i+1) << endl;
    deplacement(positions, vitesses, temps, t);
    fission(positions, vitesses, temps, comb, rayon);
    moderation(vitesses, t);

    positionsNew.clear();
    vitessesNew.clear();
    tempsNew.clear();

    Desintegration(positionsNew, vitessesNew, tempsNew, positions, vitesses, temps, dinf, hauteur);

    positions.resize(positionsNew.size());
    vitesses.resize(vitessesNew.size());
    temps.resize(tempsNew.size());
    
    positions = positionsNew;
    vitesses = vitessesNew;
    temps = tempsNew;

    sommeK += positions.size()/tmp;
    tmp = positions.size();

    ecrireDansFichier(dossier, to_string(i), positions);
  }
}


int main(int argc, char* argv[]) {
    if (argc != 9) {
        std::cerr << "Usage: " << argv[0] << " <n> <t> <r> <rep> <uma> <h> <b> <e>" << std::endl;
        return 1;
    }

    int N = 100;
    double n = std::stod(argv[1]);
    double t = std::stod(argv[2]);
    double rayon = std::stod(argv[3]);
    double rep = std::stod(argv[4]);
    double uma = std::stod(argv[5]);
    double h = std::stod(argv[6]);
    double b = std::stod(argv[7]);
    double e = std::stod(argv[8]);

    std::ostringstream oss;
    oss << n << t << rayon << rep << uma << h << b << e;
    std::string dossier = oss.str();

    vec3 p0;
    p0.x = 0.f;
    p0.y = 0.f;
    p0.z = h/2;

    vector<vec3> positions;
    vector<vec3> vitesses;
    vector<float> temps;
    vector<vec3> comb; 

    neutronvect(N, positions, vitesses, temps, p0);
    combustible(comb, rayon, h, b, e, rep);

    cycle(n, t, uma, positions, vitesses, temps, comb, rayon*rep, h, dossier);
    return 0;
}
