// Bacterial coexistence simulation
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <time.h>
#include <vector>
#include "Bacterium.h"
#include "Strain.h"
#include "FoodSource.h"

// =====================================================================

// PROGRAM PARAMETERS

// General parameters
const int grid_size = 100;   // Size of the grid (along each dimension)
const int num_cycles = 300;   // Number of cycles to simulate
const double growth_multiplier = 1.0;   // Global growth multiplier
const bool all_die = true;    // If false, only bacteria that are not at a food source die spontaneously (starvation?)
const double global_death_rate = 0.002;    // Death probability per iteration

// Strains parameters
const int num_strains = 37;   // Number of strains
const char* strains_fname = "strains.dat";   // Filename for strains data (including name, growth_chance, sender degree and receiver degree)

// Antagonism
const bool consider_antagonism = false;   // Simulate for antangonistic interactions?
const char* antagonism_table_fname = "antagonisms.dat" ;   // If simulating antagonism, provide the antangonsm table filename 

// Growth model
const bool use_real_growth = false;   // Indicates whether measured growth rates should be used (if false, will compute growth rates from sender degree)
const double growth_intercept = 0.2582;   // Growth rate for a strain with sender degree = 0
const double growth_slope = -0.0025;   // Growth rate reduction per additional sender degree

// Food sources
const int food_source_rations = 200;   // Number of rations per new food source
const int food_source_radius = 5;   // Radius of food sources
const double food_source_spawn_rate = 0.2;   // Chance of new food source appearing per iteration

// Starting conditions
const int start_population_per_strain = 100;     // Size of starting population, per strain (same for all strains)
const int start_num_food_sources = 20;   // Starting number of food source

// =====================================================================

// Global vars

// Current cycle number
int cycle;

// Counts bacteria that ever lived
int bac_counter;

// Seed of the RNG
unsigned int global_seed;

// Individual strain population size
int strain_counts[num_strains];

// List of Strains
Strain* strains[num_strains];

// List of Bacteria
std::vector<Bacterium*> bacteria_list;

// List of Food Sources
std::vector<FoodSource*> food_list;

// Antagonism matrix
// antagonism[i][j] is true if strain i antagonizes strain j
bool antagonism[num_strains][num_strains];

// =====================================================================

// Populates the strains array by reading the corresponding datafile
void create_strains() {

  std::ifstream infile;
  std::string str_buf, name;
  double growth;
  int snd_deg, rcv_deg, i;

  printf(">> Reading strain data from \"%s\" ...\n", strains_fname);
  infile.open(strains_fname, std::ifstream::in); 
  if (!infile) {
    printf("ERROR: could not open file\n");
    exit(1);
  }
  
  getline(infile, str_buf);     // skip header
  i = 0;
  while(!infile.eof() && i < num_strains) {
    infile >> str_buf >> name >> growth >> rcv_deg >> snd_deg;
    strains[i] = new Strain(i, name, growth, global_death_rate, snd_deg, rcv_deg);
//    printf("%s %f %i %i\n",  strains[i]->name.c_str(), strains[i]->growth_rate, strains[i]->receiver_degree, strains[i]->sender_degree);
    i++;
  }
  if (i != num_strains) {
    printf("ERROR: expected %i strains, but found only %i!\n", num_strains, i);
    exit(1);
  }
  printf("%d strains loaded\n", i);
  
	infile.close();
  
}

// =====================================================================

// Populates the antangonism matrix by reading the corresponding datafile
void load_antagonisms() {

  std::ifstream infile;
  double value;
  int i, j;

  printf(">> Reading antagonism matrix from \"%s\" ...\n", antagonism_table_fname);
  infile.open(antagonism_table_fname, std::ifstream::in); 
  if (!infile) {
    printf("ERROR: could not open file\n");
    exit(1);
  }
  
  i = 0;
  while(!infile.eof() && i < num_strains) {
    for (j = 0; j < num_strains; j++) {
      infile >> value;
      if (value > 1) antagonism[i][j] = true;
      else antagonism[i][j] = false;
    }
    i++;
  }
  
  // Check sender degrees
  printf("Checking sender degrees ... ");
  int deg;
  for (i = 0; i < num_strains; i++) {
    deg = 0;
    for (j = 0; j < num_strains; j++) {
//      printf("%d ", antagonism[i][j]);
      if (antagonism[i][j]) deg++;
    }
    if (deg != strains[i]->sender_degree) {
      printf("\nERROR: stored sender degree (%d) of strain %i (%s) does not match that from antagonism table (%d)\n", strains[i]->sender_degree, strains[i]->ID, strains[i]->name.c_str(), deg);
      exit(1);
    }
  }
  printf("OK\n");

  // Check receiver degrees
  printf("Checking receiver degrees ... ");
  for (i = 0; i < num_strains; i++) {
    deg = 0;
    for (j = 0; j < num_strains; j++) {
//      printf("%d ", antagonism[i][j]);
      if (antagonism[j][i]) deg++;
    }
    if (deg != strains[i]->receiver_degree) {
      printf("ERROR: stored receiver degree (%d) of strain %i (%s) does not match that from antagonism table (%d)\n", strains[i]->receiver_degree, strains[i]->ID, strains[i]->name.c_str(), deg);
      exit(1);
    }
  }
  printf("OK\n");
  
	infile.close();
  
}

// =====================================================================

// Update the growth rates from the linear growth formula
void compute_growth_rates() {
  int i;
  for (i = 0; i < num_strains; i++) {
    strains[i]->growth_rate = strains[i]->sender_degree*growth_slope + growth_intercept;
    if (strains[i]->growth_rate < 0) {
      printf("ERROR: negative growth rate for strain %i (%s)\n", strains[i]->ID, strains[i]->name.c_str());
    }
    printf("%s | %i -> %f\n", strains[i]->name.c_str(), strains[i]->sender_degree, strains[i]->growth_rate);
  }
  
}

// =====================================================================

// Adds a new bacterium of strain at (x,y)
void add_bacterium(Strain* strain, int x, int y) {
  bacteria_list.push_back(new Bacterium(bac_counter, strain, x, y));
  bac_counter++;
}

// Removes the bacterium residing at 'idx' in the bacteria_list,
// freeing its memory use
void remove_bacterium(int idx) {
  delete bacteria_list[idx];
  bacteria_list[idx] = bacteria_list.back();
  bacteria_list.pop_back();
//  bacteria_list.erase(bacteria_list.begin()+idx);
}

// =====================================================================

// Adds a new food source at (x,y)
void add_food_source(int x, int y, double radius, int rations) {
  food_list.push_back(new FoodSource(x, y, radius, rations));
}

// Removes the food source residing at 'idx' in the food_list
void remove_food_source(int idx) {
//  printf("Removing food source %i\n", idx);
  delete food_list[idx];
  food_list[idx] = food_list.back();
  food_list.pop_back();
//  food_list.erase(food_list.begin()+idx);
}

// =====================================================================

// Creates the starting bacteria and populates the list
void create_starting_bacteria() {
  
  int s, i;
  Strain* strain;
  
  for (s = 0; s < num_strains; s++) {
    strain = strains[s];
    for (i = 1; i <= start_population_per_strain; i++) {
      add_bacterium(strain, randint(1,grid_size), randint(1,grid_size));
//      printf("%zu / %zu\n", bacteria_list.size(), bacteria_list.capacity());
    }
  }
//  for (i = 0; i < bacteria_list.size(); i++) {
//    printf("%d : strain %s at %d,%d\n", i, strains[bacteria_list[i]->strain]->name.c_str(), bacteria_list[i]->x, bacteria_list[i]->y);
//  }
}

// =====================================================================

// Creates the starting food sources
void create_starting_food() {
  for (int i = 0; i < start_num_food_sources; i++) {
    add_food_source(randint(1, grid_size), randint(1, grid_size), food_source_radius, food_source_rations);
  }
}

// Spawns 
void spawn_food_sources() {
  if (randreal() <= food_source_spawn_rate) {
    add_food_source(randint(1, grid_size), randint(1, grid_size), food_source_radius, food_source_rations);
//    printf("Spawned food source at %i,%i with %i rations\n", food_list.back()->x, food_list.back()->y, food_list.back()->rations);
  }
}

// =====================================================================

// Applies the grid topology to the coordinate pair (x,y)
void apply_topology(int& x, int& y) {
  // Toroidal topology
  if (x < 1) x = grid_size;
  else if (x > grid_size) x = 1;
  if (y < 1) y = grid_size;
  else if (y > grid_size) y = 1;  
}

// =====================================================================

// Move bacteria not at a food source
void move_bacteria() {

  Bacterium* bac;
  unsigned int i, counter = 0;

  for (i = 0; i < bacteria_list.size(); i++) {
    bac = bacteria_list[i];
    if (!bac->isAtFoodSource) {
      bac->move();
      apply_topology(bac->x, bac->y);
      counter++;
    }
  }
//  printf("Moved %i / %i bacteria\n", counter, (int)bacteria_list.size());
  
}

// =====================================================================

// Kill bacteria according to their death rate
void bacteria_death() {

  Bacterium* bac;
  unsigned int idx, i, size;

  size = bacteria_list.size();
  idx = 0;
  for (i = 0; i < size; i++) {
    bac = bacteria_list[idx];
    if (all_die || !bac->isAtFoodSource) {
      if (randreal() <= bac->strain->death_rate) {
        remove_bacterium(idx);
      } else {
        idx++;
      }
    }
  }
  
}

// =====================================================================

// Check whether bacteria are at a food source, flagging them and adding
// them to the food source's list
void flag_bacteria_at_food_sources() {

  Bacterium* bac;
  FoodSource* fs;
  unsigned int i, j, counter; 

  for (i = 0; j < food_list.size(); j++)
    food_list[j]->clearBacteriaHere();

  counter = 0;
  for (i = 0; i < bacteria_list.size(); i++) {
    bac = bacteria_list[i];
    bac->isAtFoodSource = false;
      
    for (j = 0; j < food_list.size(); j++) {
      fs = food_list[j];
    
      if (fs->contains(bac->x, bac->y)) {
        fs->bacteriaHere.push_back(bac);
        if (!bac->isAtFoodSource) {
          bac->isAtFoodSource = true;
          counter++;
        }
      }
      
    }
    
  }
//  printf("Bacteria at food sources: %i\n", counter);

}

// =====================================================================

// Implements interactions that occur at food sources
// Specifically:
// - bacteria antagonize each other (if enabled)
// - bacteria eat and reproduce
// - food sources are exhausted
void food_sources_interactions() {

  FoodSource* fs = NULL;
  Bacterium* bac = NULL;
  unsigned int idx, size, num_bacs, i, j, k, will_eat;

  bool verbose = false;   // for debug
  
  size = food_list.size();
  idx = 0;
  for (k = 0; k < size; k++) {
    fs = food_list[idx];

    num_bacs = (int)fs->bacteriaHere.size();
    if (verbose) printf("Source %i at %i,%i with %i rations has %i bacteria\n", idx, fs->x, fs->y, fs->rations, num_bacs);      
    if (num_bacs == 0) continue;

    // Antagonisms
    if (consider_antagonism) {
      
    }

    // Check food
    if (fs->rations < num_bacs) {

      // If not enough food for everyone, randomly select who will eat
      will_eat = fs->rations;

      // Partial Fisher-Yates shuffle
      for (i = 0; i < will_eat; i++) {
        j = randint(i,num_bacs-1);
        if (verbose) {
          printf ("i=%i <-> j=%i ", i, j); fflush(stdout);
          printf ("| %i\n", fs->bacteriaHere[j]->strain->ID);
        }
        if (j != i) {
          bac = fs->bacteriaHere[j];
          fs->bacteriaHere[j] = fs->bacteriaHere[i];
          fs->bacteriaHere[i] = bac;
        }
      }
      
    } else {
      
      // If enough food, everyone gets to eat
      will_eat = num_bacs;
      
    }
  
    // Debug
//    for (i = 0; i < will_eat; i++) {
//      printf("bac %i ", i); fflush(stdout);
//      printf("%f | ", fs->bacteriaHere[i]->strain->growth_rate); fflush(stdout);
//    }

    if (verbose) printf("%i bacteria will eat\n", will_eat);
    
    // Now feed the bacteria, reducing the food source size, and let
    // the fed bacteria have a chance to reproduce
    if (verbose) printf("Bacteria that reproduced: ");
    for (i = 0; i < will_eat; i++) {
      fs->rations--;
      bac = fs->bacteriaHere[i];
      if (randreal() <= bac->strain->growth_rate) {
        if (verbose) printf("%i ", bac->ID);
        add_bacterium(bac->strain, randint(1,grid_size), randint(1,grid_size));
      }
    }
    if (verbose) printf("\n");
    
    // If the food source is exhausted, remove it
    if (fs->rations == 0) {
      remove_food_source(idx);
      if (verbose) printf("Removed food source (%i left)\n", (int)food_list.size());
    } else {
      idx++;
    }
    
  }

}

// =====================================================================

// Counts the number of bacteria of each strain
void compute_strain_counts() {
  int i;
  for (i = 0; i < num_strains; i++) {
    strain_counts[i] = 0;
  }
  for (i = 0; i < (int)bacteria_list.size(); i++) {
    strain_counts[bacteria_list[i]->strain->ID] += 1;
  }
}

// =====================================================================

// Reports the current state of the simulation
// If fp is given, writes the report to that file; else, to sdtout
void report(FILE* fp) {
  int i;
  if (fp == NULL) fp = stdout;
  fprintf(fp, "%i %i %i", cycle, (int)bacteria_list.size(), (int)food_list.size());
  for (i = 0; i < num_strains; i++) {
    fprintf(fp, " %i", strain_counts[i]);
  } 
  fprintf(fp, "\n");
}
void report() {
  report(NULL);
}

// =====================================================================

// Sets the seed of the RNG
void set_seed(unsigned int seed) {
  srand(seed);
  global_seed = seed;
}
// Resets the seed of the RNG (to a new random seed)
void reset_seed() {
  srand(time(NULL));
  set_seed(rand());
}

// =====================================================================
// =====================================================================
// =====================================================================

int main (int argc, char** argv) {

  FILE* fp;
  fp = fopen("test.dat","w");
 
  // Initializations
  time_t raw_time = time(NULL);
  printf(">> %s", ctime(&raw_time));
  
  // Set RNG seed
  reset_seed();
  //set_seed(1749352398);  // segfaults on Terminus
  //set_seed(964473454);  // segfaults on Seldon
  printf(">> RNG seed: %u\n", global_seed);
  
  // Create strains using data from file
  create_strains();
  
  // Load antagonism matrix (if applicable)
  if (consider_antagonism) {
    printf(">> Antagonisms ENABLED\n");
    load_antagonisms();
  } else {
    printf(">> Antagonisms DISABLED\n");
  }

  // Compute growth_rates if not using real ones
  if (use_real_growth) {
    printf(">> Using measured growth rates\n");
  } else {
    printf(">> Computing growth rates using linear relation ...\n");
    compute_growth_rates();
  }
  
  // Create starting bacteria
  bac_counter = 0;
  printf(">> Creating bacteria ...\n");
  create_starting_bacteria();
  printf("Population: %i\n", (int)bacteria_list.size());
  compute_strain_counts();
  
  // Create food sources
  printf(">> Creating food sources ...\n");
  create_starting_food();
  printf("Food sources: %i\n", (int)food_list.size());
  flag_bacteria_at_food_sources();

  // Ready
  printf(">> Starting simulation\n");
  report();
//  getchar();
 
  // Main loop over cycles
  for (cycle = 1; cycle <= num_cycles; cycle++) {

    move_bacteria();

    flag_bacteria_at_food_sources();

    // For each food source:
    // - Check antagonisms, killing antagonized bacteria
    // - Have bacteria eat (just the lucky ones if not enough food)
    // - Have bacteria that ate reproduce
    // - Remove food source if emptied
    food_sources_interactions();

    bacteria_death();
    
    spawn_food_sources();

    compute_strain_counts();

//    printf("\n");
    report();
  
  }
  
}
