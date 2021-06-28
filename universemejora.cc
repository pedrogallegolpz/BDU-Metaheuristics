extern "C" {
#include "cec17.h"
}
#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <random>


////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
//  
//                  FUNCIONES AUXILIARES
//
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
/*
    Calcula la norma 2 de un vector
*/
double l2_norm(std::vector<double> &v){
    double norm = 0;

    for(int i=0; i<v.size(); i++){
        norm += v[i]*v[i];
    }
    norm = sqrt(norm);

    return norm;
}


/*
    Calcula la distancia euclidea entre dos puntos.
*/
double euclidean_distance(std::vector<double> &p1, std::vector<double> &p2){
    double dist = 0;

    if(p1.size()!=p2.size()){
        std::cerr <<"Warning at euclidean_distance: p1 has different size than p2.\n";
    }

    for(int i=0; i<p1.size(); i++){
        dist += (p1[i]-p2[i])*(p1[i]-p2[i]);
    }
    dist = sqrt(dist);

    return dist;
}


void clip(std::vector<double> &sol, int lower, int upper) {
  for (auto &val : sol) {
    if (val < lower) {
      val = lower;
    }
    else if (val > upper) {
      val = upper;
    }
  }
}

void increm_bias(std::vector<double> &bias, std::vector<double> dif) {
  for (unsigned i = 0; i < bias.size(); i++) {
    bias[i] = 0.2*bias[i]+0.4*(dif[i]+bias[i]);
  }
}

void decrement_bias(std::vector<double> &bias, std::vector<double> dif) {
  for (unsigned i = 0; i < bias.size(); i++) {
    bias[i] = bias[i]-0.4*(dif[i]+bias[i]);
  }
}


////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
//  
//                  OPERADORES
//
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
/*
    Perturbamos todas las coordenadas ligeramente
*/
void create_nghbr(std::vector<double> & pivot, std::vector<double> & neighbour, double dist, int seed){
    // Iniciamos distribución para valores aleatorios
    std::uniform_real_distribution<> dis(-dist, dist);
    std::mt19937 gen(seed); // Inicio semilla

    // Modificamos únicamente un tercio del planeta;
    for(int i=0; i<pivot.size(); i++){
        neighbour[i]=pivot[i]+dis(gen);
    }

}


/*
    Orbitar
*/
void orbit(std::vector<double> &planet, std::vector<double> &star, int seed){
    std::vector<double> movement(planet.size());
    double distancia_objetivo = 0.5*euclidean_distance(planet, star);
    double dist_obj_cuadrado_restante = distancia_objetivo*distancia_objetivo;


    // Iniciamos distribución para valores aleatorios
    // set some values:
    std::vector<int> indices;
    for (int i=0; i<planet.size(); ++i) indices.push_back(i); 

    // using built-in random generator:
    std::random_shuffle ( indices.begin(), indices.end() );    
    
    std::uniform_real_distribution<> dis_valor(-1,1);
    std::mt19937 gen(seed); // Inicio semilla

    for(int ii=0; ii<planet.size(); ii++){
        int i = indices[ii];
        double valor = dis_valor(gen);
        planet[i] = star[i]+valor*sqrt(dist_obj_cuadrado_restante);

        dist_obj_cuadrado_restante -= valor*valor*dist_obj_cuadrado_restante;
    }


}



/**
 * Aplica el Solis Wets
 *
 * @param  sol solucion a mejorar.
 * @param fitness fitness de la solución.
 */
template <class Random>
void soliswets(std::vector<double> &sol, double &fitness, double delta, int maxevals, int lower, int upper, Random &random) {
  const size_t dim = sol.size();
  std::vector<double> bias (dim), dif (dim), newsol (dim);
  double newfit;
  size_t i;

  int evals = 0;
  int num_success = 0;
  int num_failed = 0;

  while (evals < maxevals) {
    std::uniform_real_distribution<double> distribution(0.0, delta);

    for (i = 0; i < dim; i++) {
      dif[i] = distribution(random);
      newsol[i] = sol[i] + dif[i] + bias[i];
    }

    clip(newsol, lower, upper);
    newfit = cec17_fitness(&newsol[0]);
    evals += 1;

    if (newfit < fitness) {
      sol = newsol;
      fitness = newfit;
      increm_bias(bias, dif);
      num_success += 1;
      num_failed = 0;
    }
    else if (evals < maxevals) {

      for (i = 0; i < dim; i++) {
        newsol[i] = sol[i] - dif[i] - bias[i];
      }

      clip(newsol, lower, upper);
      newfit = cec17_fitness(&newsol[0]);
      evals += 1;

      if (newfit < fitness) {
        sol = newsol;
        fitness = newfit;
        decrement_bias(bias, dif);
        num_success += 1;
        num_failed = 0;
      }
      else {
        for (i = 0; i < dim; i++) {
          bias[i] /= 2;
        }

        num_success = 0;
        num_failed += 1;
      }
    }

    if (num_success >= 5) {
      num_success = 0;
      delta *= 2;
    }
    else if (num_failed >= 3) {
      num_failed = 0;
      delta /= 2;
    }
  }

}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
//  
//                  ACCIONES
//
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
/*

*/
void universe_init_(std::vector<std::vector<std::vector<double>>> & universe, int dim, int n_planetas, int n_galaxias, int seed){
    // Creamos vector auxiliar de planetas, vector auxiliar de galaxias y el universo que va a guardar la información de todo
    std::vector<double> planeta(dim);
    std::vector<std::vector<double>> galaxia(n_planetas);

    // Iniciamos distribución para valores aleatorios
    std::uniform_real_distribution<> dis(-100.0, 100.0);
    std::mt19937 gen(seed); // Inicio semilla
    
    // Generamos una solución aleatoria por galaxia
    for(int i=0; i<n_galaxias; i++){
        for(int k=0; k<dim; k++){
            planeta[k] = dis(gen);
        }
        galaxia[0] = planeta;
        universe[i] = galaxia;
    }


    
    // Creamos planetas cercanos a los planetas iniciales de cada galaxia para obtener el resto de planetas
    for(int i=0; i<n_galaxias; i++){
        for(int j=1; j<n_planetas; j++){
            // Creamos los nuevos planetas vecinos en cada galaxia
            create_nghbr(universe[i][0], planeta, j, seed+i+j*n_planetas);

            // Lo guardamos como planeta de la galaxia
            universe[i][j]=planeta;
        }
    }

}


/*
    Comienza el proceso de órbita de los planetas
        - galaxy:   Conjunto de planetas de una determinada galaxia
        - star:     Planeta con mejor fitness. Planeta eje sobre el que giran los demás
*/
void run_planets_orbit(std::vector<std::vector<double>> &galaxy, int &star, std::vector<double> &fitness_galaxy, int max_tries, int &evals, int seed){

    for(int i=0; i<galaxy.size(); i++){
        if(i!=star){
            double distancia = euclidean_distance(galaxy[i], galaxy[star]);
            double fitness;

            std::vector<double> old = galaxy[i];

            if(distancia>=1){
                int tries=0;
                bool ha_mejorado=false;
                bool seguir;
                do{
                    // Orbitamos
                    orbit(galaxy[i], galaxy[star], seed);
                    
                    // Actualizamos resultados
                    // Evaluamos la inicialización
                    fitness = cec17_fitness(&galaxy[i][0]);
                    evals++;
                    tries++;
                    seed+=i*10+1;

                    if(fitness<fitness_galaxy[star] && ha_mejorado==false){
                        ha_mejorado=true;
                        star = i;
                        fitness_galaxy[i]=fitness;
                    }

                    
                    //seguir = ((fitness_galaxy[star]<=fitness && ha_mejorado==false) || (fitness_galaxy[star]>fitness && ha_mejorado==true));
                    seguir = ha_mejorado==false;
                }while( seguir==true && tries<max_tries );

                if(ha_mejorado==false){
                    // Recuperamos la info
                    galaxy[i]=old;
                }

            }else{
                // Si tienen mismas coordenadas les separamos haciendo un vecino nuevo
                create_nghbr(galaxy[star], galaxy[i], galaxy[i].size(), seed*13*17);
                
                // Actualizamos resultados
                // Evaluamos la inicialización
                fitness_galaxy[i]  = cec17_fitness(&galaxy[i][0]);
                evals++;
            }

            seed++;
            
            // Si es la mejor de todas se guarda
            if(fitness_galaxy[i]<fitness_galaxy[star]){
                star = i;
            }
        }else{
            // Búsqueda local para codificación real: método Solis Wets
            std::mt19937 gen(seed*10+1); // Inicio semilla
            
            double old_fitness = fitness_galaxy[i];

            soliswets(galaxy[i], fitness_galaxy[i], 3, 5*max_tries,-100,100, gen);

            if(fitness_galaxy[i]>=old_fitness){
                // Si ha mejorado, actualizamos la info
                fitness_galaxy[i] = old_fitness;
            }
        }
    }
}


bool atraer(std::vector<std::vector<double>> &galaxy, std::vector<double> &pto_atraccion, int &star,  std::vector<double> &fitness_galaxy, int &evals, double velocidad){
    // set some values:
    std::vector<int> indices;
    for (int i=0; i<galaxy[0].size(); ++i) indices.push_back(i); 

    int n_cambios = velocidad * galaxy[0].size();

    double old_fit = fitness_galaxy[star];
    // using built-in random generator:
    std::random_shuffle ( indices.begin(), indices.end() );
    for(int ii=0; ii<galaxy[0].size(); ii++){
        int i = indices[ii];
        for(int j=0; j<galaxy.size(); j++){
            if(ii<n_cambios){
                galaxy[j][i]=pto_atraccion[i];
            }
        }
    }

    //Actualizamos fitness
    for(int i=0; i<galaxy.size(); i++){
        fitness_galaxy[i]=cec17_fitness(&galaxy[i][0]);
        evals++;
        if(i==0 || fitness_galaxy[i]<fitness_galaxy[star]){
            star=i;
        }
    }
    bool mejora=false;
    if(old_fit>fitness_galaxy[star]){
        mejora = true;
    }

    return mejora;
    
}



////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
//  
//                  MAIN
//
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
int main() {

    std::vector<int> dims={10,30};

    int n_ejecuciones=10;
    std::vector<int> seeds(n_ejecuciones);
    for(int i=0; i<n_ejecuciones; i++){
        seeds[i]=i+1;
    }
    int n_galaxias = 50;
    int n_planetas = 7;

    bool verbose = true;

    std::uniform_real_distribution<> dis(-100.0, 100.0);

    for(int d=0; d<dims.size();d++){
        int dim=dims[d];
        int epoca=0;
        for (int funcid = 1; funcid <= 30; funcid++) {
            epoca++;

            if(verbose){
                std::cout << "Ejecutando función: " << funcid << std::endl;
            }
            for(int s=0; s<n_ejecuciones; s++){       
                int seed = seeds[s];

                std::vector<double> bestsol(dim);

                // Creamos el universo:
                //          - La idea es que el mejor planeta de la galaxia i, sea el planeta j==best_planet_galaxys(i),
                //            el cual tiene un score o fitness: score==fitnes_galaxys(i)
                int galaxias_activas = n_galaxias;
                std::vector<int> star_galaxys(n_galaxias);            // Guardamos el índice del mejor planeta de cada Galaxia
                std::vector<std::vector<std::vector<double>>> universo(n_galaxias);    
                std::vector<std::vector<double>> fit_univ(n_galaxias);
                std::vector<double> learning_rate(n_galaxias);
                for(int i=0; i<n_galaxias;i++){
                    fit_univ[i].resize(n_planetas);
                }   
                

                double fitness;
                double best = -1;
                double best_galaxy = -1;

                // Inicializamos los datos para nuestro algoritmo
                cec17_init("universemejora", funcid, dim);


                //std::cerr <<"Warning: output by console, if you want to create the output file you have to comment cec17_print_output()\n";
                //cec17_print_output(); // Comment to generate the output file

                

                int evals = 0;
                ////////////////////////////////////////////////////////////
                // Inicializamos el universo
                universe_init_(universo, dim, n_planetas, n_galaxias, seed);

                // Evaluamos la inicialización
                for(int i=0; i<universo.size(); i++){
                    for(int j=0; j<universo[0].size(); j++){
                        // Evalúo el planeta j-ésimo de la i-ésima galaxia
                        fit_univ[i][j] = cec17_fitness(&universo[i][j][0]);
                        evals++;
                        
                        // Si es la primera evaluación de la galaxia, o es la mejor de todas se guarda
                        if(j==0 || fit_univ[i][j] < fit_univ[i][star_galaxys[i]]){
                            star_galaxys[i]=j;

                            // Si es la primera evaluación del universo o es la mejor de todo el universo, se guarda
                            if (evals == 1 || fit_univ[i][j] < best) {
                                best = fit_univ[i][j];
                                bestsol = universo[i][j];
                                best_galaxy = i;
                            }
                        }
                    }
                    learning_rate[i]=0.5;
                }

                

                // Ejecutamos las épocas
                int max_evaluaciones = 10000*dim;
                while(evals < max_evaluaciones ) {
                    for(int i=0; i<universo.size() && evals<max_evaluaciones; i++){
                        // Hacemos orbitar los planetas
                        int intentos = 10;
                        run_planets_orbit(universo[i], star_galaxys[i], fit_univ[i], intentos,evals, seed);

                        // Actualizamos resultados en el universo
                        if(fit_univ[i][star_galaxys[i]] < best){
                            best = fit_univ[i][star_galaxys[i]];
                            bestsol = universo[i][star_galaxys[i]];
                            best_galaxy = i;
                        }

                        seed++;
                    }                    
                    

                    // Efectuamos la fuerza gravitacional del agujero negro
                    for(int i=0; i<universo.size(); i++){
                        double velocidad=learning_rate[i];

                        bool mejora = atraer(universo[i], bestsol, star_galaxys[i], fit_univ[i], evals, velocidad);

                        // Si mejora actualizamos el learning rate
                        if(mejora){
                            learning_rate[i]+=(1-learning_rate[i])*0.25;
                        }else{
                            learning_rate[i]-=(learning_rate[i])*0.25;
                        }

                        // Actualizamos resultados en el universo
                        if(fit_univ[i][star_galaxys[i]] < best){
                            best = fit_univ[i][star_galaxys[i]];
                            bestsol = universo[i][star_galaxys[i]];
                            best_galaxy = i;
                        }

                        double distancia = euclidean_distance(universo[i][star_galaxys[i]], bestsol);

                        // Comprobamos si algún planeta se ha salido del espacio
                        for(int j=0; j<universo[i][star_galaxys[i]].size(); j++){
                            if(universo[i][star_galaxys[i]][j]>=100 || universo[i][star_galaxys[i]][j]<=-100 || best_galaxy==i || distancia<0.01){
                                // Damos por perdida la galaxia y creamos una nueva galaxia centrada
                                // Iniciamos distribución para valores aleatorios
                                std::uniform_real_distribution<> dis(-100.0, 100.0);
                                std::mt19937 gen(seed+evals); // Inicio semilla
                                
                                // Generamos una solución aleatoria por galaxia
                                for(int k=0; k<dim; k++){
                                    universo[i][0][k] = dis(gen); //bestsol[k]+dis(gen)*(1-evals/max_evaluaciones);
                                    if(universo[i][0][k]>100){
                                        universo[i][0][k]=100;
                                    }else if(universo[i][0][k]<-100){
                                        universo[i][0][k]=-100;
                                    }
                                }
                                fit_univ[i][0] = cec17_fitness(&universo[i][0][0]);
                                
                                // Creamos planetas cercanos a los planetas iniciales de cada galaxia para obtener el resto de planetas
                                for(int j=1; j<n_planetas; j++){
                                    // Creamos los nuevos planetas vecinos en cada galaxia
                                    create_nghbr(universo[i][0], universo[i][j], j, seed+evals);
                                    fit_univ[i][j] = cec17_fitness(&universo[i][j][0]);
                                }

                                learning_rate[i]=0.5;

                            }
                        }
                    }
                

                }   
                ////////////////////////////////////////////////////////////


                std::cout <<"Best Universe[F" <<funcid <<"]: " << std::scientific <<cec17_error(best) << std::endl;
            }
        }
    }
}
