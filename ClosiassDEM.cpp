#include <complex>
#include <iterator>
#include <string>
#include <vector>
#include <iostream>
#include <random>
#include <cmath>
#include <algorithm>

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
  int k = 0;
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
  vector <float> dist;
  for (int i = 0; i < positions.size(); i++) {
    distance(dist, positions, positions[i]);
    int counter = 0; 
    for (int j = 0; j < dist.size(); j++) {
      if (dist[j] > 0 && dist[j] <= dmin){
        counter++;
        if (counter > bulle_max){
          temps[j] = std::numeric_limits<float>::max();
        }
      }
    }
  }
}

void Desintegration(vector<vec3> &positionsNew, vector<vec3> &vitessesNew, vector<float> &tempsNew, const vector<vec3> &positions, const vector<vec3> &vitesses, const vector<float> &temps, float dinf, float hauteur){
    float AgeFermi = 400.f;

    for (int i = 0; i < positions.size(); i++) {
        vec3 PointCentre;
        PointCentre.x = 0;
        PointCentre.y = 0;
        PointCentre.z = hauteur / 2;
        
        float dist = pow(positions[i].x - PointCentre.x, 2.f) +
                     pow(positions[i].y - PointCentre.y, 2.f) +
                     pow(positions[i].z - PointCentre.z, 2.f);

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
  float m = 0.1;
  for (int i = 0; i < positions.size(); i++) {
    vector<float> dist;
    vector<float> nearbyindices;
    
    distance(dist, comb, positions[i]);
    
    for (int j = 0; j < dist.size(); i++) {
      if (dist[i] <= rayon) {
        nearbyindices.push_back(dist[i]);  
      }
    }
    if (nearbyindices.size() != 0) {
      auto min = std::min_element(nearbyindices.begin(), nearbyindices.end());
      int noyau = std::distance(nearbyindices.begin(), min);
    
      float Ec = pow(vitesses[i].x, 2.0f) +
                 pow(vitesses[i].y, 2.0f) +
                 pow(vitesses[i].z, 2.0f);
      Ec = Ec * m;
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
        positions[i].x *= normalDist2(gen);
        positions[i].y *= normalDist2(gen);
        positions[i].z *= normalDist2(gen);
        temps[i] = 0.f;
      }
    }
  }
}
void cycle(int n, float t, float rayon, vector<vec3> &positions, vector<vec3> &vitesses, vector<float> &temps, vector<vec3> &comb, float dinf, float hauteur, string dossier);

/*
void Desintegration(vector<vec3> &positionsNew, vector<vec3> &vitessesNew, vector<float> &tempsNew, vector<vec3> &positions, vector<vec3> &vitesses, vector<float> &temps, float dinf, float hauteur){
  float AgeFermi = 2.f;
  DesintegrationAux(positions, temps);
  
  vec3 PointCentre;
  
  PointCentre.x = 0;
  PointCentre.y = 0;
  PointCentre.z = hauteur/2;
  
  vector<float> dist;
  distance(dist, positions, PointCentre);
  
  for (int i = 0; i < positions.size(); i++) {  
    if (temps[i] < AgeFermi && dist[i] < dinf){
      positionsNew.push_back(positions[i]);
      vitessesNew.push_back(vitesses[i]);
      tempsNew.push_back(temps[i]);
    }
  }
}
*/


int main(int argc, char* argv[]) {
    vec3 p0;
    int N = 10;
    p0.x = 0.f;
    p0.y = 0.f;
    p0.z = 1.f;
    vector<vec3> positions;
    vector<vec3> vitesses;
    vector<float> temps;
    vector<vec3> comb; 
    neutronvect(N, positions, vitesses, temps, p0);
    combustible(comb, 1.f, 2.f,6.f , 2.f, 3);
    vector<vec3> positionsNew;
    vector<vec3> vitessesNew;
    vector<float> tempsNew;
  

    for (int i = 0; i < 12 ; i++){

    // Affichage des résultats pour les premières particules
    /*
    for (int i = 0; i < 10; ++i) {
        cout << "Particule " << i << " - Position : (" << positions[i].x << ", " << positions[i].y << ", " << positions[i].z << ")" << endl;
        cout << "Particule " << i << " - Vitesse : (" << vitesses[i].x << ", " << vitesses[i].y << ", " << vitesses[i].z << ")" << endl;
        cout << "Particule " << i << " - Temps : " << temps[i] << endl;
    }
    */
    moderation(vitesses, 2.1f);
    /*
    for (size_t i = 0; i < positions.size(); ++i) {
        std::cout << "Nouvelle position " << i << ": (" << positions[i].x << ", " << positions[i].y << ", " << positions[i].z << ")" << std::endl;
        std::cout << "Nouveau temps " << i << ": " << temps[i] << std::endl;
        std::cout << "Nouvelle vitesse " << i << ": (" << vitesses[i].x << ", " << vitesses[i].y << ", " << vitesses[i].z << ")" << std::endl;
    }
    */
    //deplacement(positions, vitesses, temps, 2.1f);
    /*
    for (size_t i = 0; i < positions.size(); ++i) {
        std::cout << "Nouvelle position " << i << ": (" << positions[i].x << ", " << positions[i].y << ", " << positions[i].z << ")" << std::endl;
        std::cout << "Nouveau temps " << i << ": " << temps[i] << std::endl;
    }
    */

    
    }
    for (int i = 0; i < 1; i++) {
        Desintegration(positionsNew, vitessesNew, tempsNew, positions, vitesses, temps, 1000.f, 1.f);

        // For debugging purposes, print the size of positionsNew after each iteration
        //cout << "Iteration " << i << " - PositionsNew size: " << positionsNew.size() << endl;

        // Further processing or output if needed
    }
    return 0;
}



