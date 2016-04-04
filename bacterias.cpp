// Bacterial coexistence simulation
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <string>
#include <time.h>
#include <vector>
#include "Bacterium.h"
#include "Strain.h"
#include "FoodSource.h"

// =====================================================================

// PROGRAM PARAMETERS -- values to be read from file

// General parameters
int    grid_size;           // Size of the grid (along each dimension)
int    num_cycles;          // Number of cycles to simulate
double growth_multiplier;   // Global growth multiplier
bool   all_die;             // If false, only bacteria that are not at a food source die spontaneously (starvation?)
double global_death_rate;   // Death probability per iteration

// Output parameters
bool  report_to_screen;   // Report simulation progress to screen?
char* log_fname;          // Filename for logging
char* output_fname;       // Filename for data output
int   output_interval;    // Number of iterations between data outputs

// Strains parameters
int   num_strains;      // Number of strains
char* strains_fname;    // Filename for strains data (including name, growth_chance, sender degree and receiver degree)

// Antagonism
bool  consider_antagonism;       // Simulate for antangonistic interactions?
char* antagonism_table_fname;    // If simulating antagonism, provide the antangonsm table filename 

// Growth model
bool   use_real_growth;    // Indicates whether measured growth rates should be used (if false, will compute growth rates from sender degree)
double growth_intercept;   // Growth rate for a strain with sender degree = 0
double growth_slope;       // Growth rate reduction per additional sender degree

// Food sources
int    food_source_rations;      // Number of rations per new food source
int    food_source_radius;       // Radius of food sources
double food_source_spawn_rate;   // Chance of new food source appearing per iteration

// Starting conditions
int start_population_per_strain;   // Size of starting population, per strain (same for all strains)
int start_num_food_sources;        // Starting number of food source

// =====================================================================

// Global vars

// Current cycle number
int cycle;

// Counts bacteria that ever lived
int bac_counter;

// Seed of the RNG
unsigned int global_seed;

// Individual strain population size
int* strain_counts;

// List of Strains
Strain** strains;

// List of Bacteria
std::vector<Bacterium*> bacteria_list;

// List of Food Sources
std::vector<FoodSource*> food_list;

// Antagonism matrix
// antagonism[i][j] is true if strain i antagonizes strain j
bool** antagonisms;

// Filenames
char* params_fname;

// Log and data output files
FILE* log_file;
FILE* output_file;

// String buffer
char* strbuf;

// =====================================================================

// Populates the strains array by reading the corresponding datafile
void create_strains() {

  std::ifstream infile;
  std::string str_buf, name;
  double growth;
  int snd_deg, rcv_deg, i;

  fprintf(log_file, ">> Reading strain data from \"%s\" ...\n", strains_fname);
  infile.open(strains_fname, std::ifstream::in); 
  if (!infile) {
    fprintf(log_file, "ERROR: could not open file\n");
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
    fprintf(log_file, "ERROR: expected %i strains, but found only %i!\n", num_strains, i);
    exit(1);
  }
  fprintf(log_file, "%d strains loaded\n", i);
  
	infile.close();
  
}

// =====================================================================

// Populates the antangonism matrix by reading the corresponding datafile
void load_antagonisms() {

  std::ifstream infile;
  double value;
  int i, j;

  fprintf(log_file, ">> Reading antagonism matrix from \"%s\" ...\n", antagonism_table_fname);
  infile.open(antagonism_table_fname, std::ifstream::in); 
  if (!infile) {
    fprintf(log_file, "Error: could not open file\n");
    exit(1);
  }
  
  i = 0;
  while(!infile.eof() && i < num_strains) {
    for (j = 0; j < num_strains; j++) {
      infile >> value;
      if (value > 1) antagonisms[i][j] = true;
      else antagonisms[i][j] = false;
    }
    i++;
  }
  
  // Check sender degrees
  fprintf(log_file, "Checking sender degrees ... ");
  int deg;
  for (i = 0; i < num_strains; i++) {
    deg = 0;
    for (j = 0; j < num_strains; j++) {
//      printf("%d ", antagonism[i][j]);
      if (antagonisms[i][j]) deg++;
    }
    if (deg != strains[i]->sender_degree) {
      fprintf(log_file, "\nError: stored sender degree (%d) of strain %i (%s) does not match that from antagonism table (%d)\n", strains[i]->sender_degree, strains[i]->ID, strains[i]->name.c_str(), deg);
      exit(1);
    }
  }
  printf("OK\n");

  // Check receiver degrees
  fprintf(log_file, "Checking receiver degrees ... ");
  for (i = 0; i < num_strains; i++) {
    deg = 0;
    for (j = 0; j < num_strains; j++) {
//      printf("%d ", antagonisms[i][j]);
      if (antagonisms[j][i]) deg++;
    }
    if (deg != strains[i]->receiver_degree) {
      fprintf(log_file, "Error: stored receiver degree (%d) of strain %i (%s) does not match that from antagonism table (%d)\n", strains[i]->receiver_degree, strains[i]->ID, strains[i]->name.c_str(), deg);
      exit(1);
    }
  }
  fprintf(log_file, "Degrees OK\n");
  
	infile.close();
  
}

// =====================================================================

// Update the growth rates from the linear growth formula
void compute_growth_rates() {
  int i;
  for (i = 0; i < num_strains; i++) {
    strains[i]->growth_rate = strains[i]->sender_degree*growth_slope + growth_intercept;
    if (strains[i]->growth_rate < 0) {
      fprintf(log_file, "ERROR: negative growth rate for strain %i (%s)\n", strains[i]->ID, strains[i]->name.c_str());
    }
    fprintf(log_file, "%s | %i -> %f\n", strains[i]->name.c_str(), strains[i]->sender_degree, strains[i]->growth_rate);
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
    if (verbose) fprintf(log_file, "Source %i at %i,%i with %i rations has %i bacteria\n", idx, fs->x, fs->y, fs->rations, num_bacs);      
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
          fprintf(log_file, "i=%i <-> j=%i ", i, j); fflush(stdout);
          fprintf(log_file, "| %i\n", fs->bacteriaHere[j]->strain->ID);
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

    if (verbose) printf("%i bacteria will eat\n", will_eat);
    
    // Now feed the bacteria, reducing the food source size, and let
    // the fed bacteria have a chance to reproduce
    if (verbose) printf("Bacteria that reproduced: ");
    for (i = 0; i < will_eat; i++) {
      fs->rations--;
      bac = fs->bacteriaHere[i];
      if (randreal() <= bac->strain->growth_rate) {
        if (verbose) fprintf(log_file, "%i ", bac->ID);
        add_bacterium(bac->strain, randint(1,grid_size), randint(1,grid_size));
      }
    }
    if (verbose) printf("\n");
    
    // If the food source is exhausted, remove it
    if (fs->rations == 0) {
      remove_food_source(idx);
      if (verbose) fprintf(log_file, "Removed food source (%i left)\n", (int)food_list.size());
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

// Does dynamic allocation of data arrays
void do_allocations() {
  
  int i;
  
  strain_counts = new int[num_strains];
  strains = new Strain*[num_strains];
  antagonisms = new bool*[num_strains];
  for (i = 0; i < num_strains; i++) {
    antagonisms[i] = new bool[num_strains];
  }
  
}

// Deallocates data arrays
void deallocate() {
  
  int i;
  delete strain_counts;
  delete strains;
  for (i = 0; i < num_strains; i++) {
    delete antagonisms[i];
  }
  delete antagonisms;
  
}

// =====================================================================

// Returns the next data-holding line, skipping comments and blank lines
std::string get_next_line(std::ifstream& file, const char* param_name) {

  char firstchar;
  std::string line;
  int i;

  while (std::getline(file, line)) {

//    std::cout << "--" << line << std::endl;
    
    if (line.length() == 0) continue;
    firstchar = '\0';
    for (i = 0; i < line.length(); i++) {
      if (line[i] != ' ') {
        firstchar = line[i];
        break;
      }
    }
    if (firstchar == '\0' or firstchar == '#') continue;
    else return line;
    
  }
  
  printf("Error! Reached end of file while looking for '%s'\n", param_name);
  exit(EXIT_FAILURE);
  
}

// Load parameters from the given filename
void load_parameters(const char* params_fname) {

  std::ifstream params_file(params_fname);
  std::string line;

  line = get_next_line(params_file, "grid_size");
  grid_size = atoi(line.c_str());
  
  line = get_next_line(params_file, "num_cycles");
  num_cycles = atoi(line.c_str());
  
  line = get_next_line(params_file, "growth_multiplier");
  growth_multiplier = atof(line.c_str());
  
  line = get_next_line(params_file, "global_death_rate");
  global_death_rate = atof(line.c_str());

  line = get_next_line(params_file, "all_die");
  all_die = (atoi(line.c_str()) == 1);
  
  line = get_next_line(params_file, "report_to_screen");
  report_to_screen = (atoi(line.c_str()) == 1);
  
  line = get_next_line(params_file, "log_fname");
  log_fname = new char[strlen(line.c_str())];
  strcpy(log_fname, line.c_str());

  line = get_next_line(params_file, "output_fname");
  output_fname = new char[strlen(line.c_str())];
  strcpy(output_fname, line.c_str());
  
  line = get_next_line(params_file, "output_interval");
  output_interval = atoi(line.c_str());
  
  line = get_next_line(params_file, "num_strains");
  num_strains = atoi(line.c_str());

  line = get_next_line(params_file, "strains_fname");
  strains_fname = new char[strlen(line.c_str())];
  strcpy(strains_fname, line.c_str());
  
  line = get_next_line(params_file, "consider_antagonism");
  consider_antagonism = (atoi(line.c_str()) == 1);
  
  line = get_next_line(params_file, "antagonism_table_fname");
  antagonism_table_fname = new char[strlen(line.c_str())];
  strcpy(antagonism_table_fname, line.c_str());
  
  line = get_next_line(params_file, "use_real_growth");
  use_real_growth = (atoi(line.c_str()) == 1);
  
  line = get_next_line(params_file, "growth_intercept");
  growth_intercept = atof(line.c_str());

  line = get_next_line(params_file, "growth_slope");
  growth_slope = atof(line.c_str());

  line = get_next_line(params_file, "food_source_rations");
  food_source_rations = atoi(line.c_str());

  line = get_next_line(params_file, "food_source_radius");
  food_source_radius = atoi(line.c_str());

  line = get_next_line(params_file, "food_source_spawn_rate");
  food_source_spawn_rate = atof(line.c_str());

  line = get_next_line(params_file, "start_population_per_strain");
  start_population_per_strain = atoi(line.c_str());
  
  line = get_next_line(params_file, "start_num_food_sources");
  start_num_food_sources = atoi(line.c_str());
  
}
// =====================================================================

// Writes the current state of the simulation to the data ouput file
void data_output() {
  int i;
  fprintf(output_file, "%i %i %i", cycle, (int)bacteria_list.size(), (int)food_list.size());
  for (i = 0; i < num_strains; i++) {
    fprintf(output_file, " %i", strain_counts[i]);
  } 
  fprintf(output_file, "\n");
}

// =====================================================================

// Logs a message to the logfile, echoing to screen if enabled
void log(const char* message) {
  
  fputs(message, log_file);
  if (report_to_screen) fputs(message, stdout);
  
}

// =====================================================================
// =====================================================================
// =====================================================================

int main (int argc, char** argv) {

  // Load parameters from the specified file
  if (argc != 2) {
    printf("Error! Must provide a parameters filename as first argument\n");
    exit(EXIT_FAILURE);
  }
  params_fname = new char[strlen(argv[1])];
  strcpy(params_fname, argv[1]);
  printf(">> Reading %s ...\n", params_fname);
  load_parameters(params_fname);
  printf(">> Parameters read successfully\n");
  
  // Open files
  log_file = fopen(log_fname,"w");
  output_file = fopen(output_fname,"w");
  
  if (!report_to_screen) {
    printf("Writing further log output to %s\n", log_fname);
  }
  strbuf = new char[256];

  // Report parameters
  sprintf(strbuf, ">> Loaded parameters from %s\n", params_fname); log(strbuf);
  sprintf(strbuf, "grid_size = %i\n", grid_size); log(strbuf);
  sprintf(strbuf, "num_cycles = %i\n", num_cycles); log(strbuf);
  sprintf(strbuf, "growth_multiplier = %f\n", growth_multiplier); log(strbuf);
  sprintf(strbuf, "global_death_rate = %f\n", global_death_rate); log(strbuf);
  sprintf(strbuf, "all_die = %s\n", all_die ? "true" : "false"); log(strbuf);
  sprintf(strbuf, "report_to_screen = %s\n", report_to_screen ? "true" : "false"); log(strbuf);
  sprintf(strbuf, "log_fname = %s\n", log_fname); log(strbuf);
  sprintf(strbuf, "output_fname = %s\n", output_fname); log(strbuf);
  sprintf(strbuf, "output_interval = %i\n", output_interval); log(strbuf);
  sprintf(strbuf, "num_strains = %i\n", num_strains); log(strbuf);
  sprintf(strbuf, "strains_fname = %s\n", strains_fname); log(strbuf);
  sprintf(strbuf, "consider_antagonism = %s\n", consider_antagonism ? "true" : "false"); log(strbuf);
  sprintf(strbuf, "antagonism_table_fname = %s\n", antagonism_table_fname); log(strbuf);
  sprintf(strbuf, "use_real_growth = %s\n", use_real_growth ? "true" : "false"); log(strbuf);
  sprintf(strbuf, "growth_intercept = %f\n", growth_intercept); log(strbuf);
  sprintf(strbuf, "growth_slope = %f\n", growth_slope); log(strbuf);
  sprintf(strbuf, "food_source_rations = %i\n", food_source_rations); log(strbuf);
  sprintf(strbuf, "food_source_radius = %i\n", food_source_radius); log(strbuf);
  sprintf(strbuf, "food_source_spawn_rate = %f\n", food_source_spawn_rate); log(strbuf);
  sprintf(strbuf, "start_population_per_strain = %i\n", start_population_per_strain); log(strbuf);
  sprintf(strbuf, "start_num_food_sources = %i\n", start_num_food_sources); log(strbuf);

  // Begin initializations
  time_t raw_time = time(NULL);
  sprintf(strbuf, ">> %s", ctime(&raw_time)); log(strbuf);

  // Allocate data arrays
  do_allocations();
  log(">> Data arrays allocated\n");
  
  // Set RNG seed
  reset_seed();
  sprintf(strbuf, ">> RNG seed: %u\n", global_seed); log(strbuf);
  
  // Create strains using data from file
  create_strains();
  
  // Load antagonism matrix (if applicable)
  if (consider_antagonism) {
    log(">> Antagonisms ENABLED\n");
    load_antagonisms();
  } else {
    log(">> Antagonisms DISABLED\n");
  }

  // Compute growth_rates if not using real ones
  if (use_real_growth) {
    log(">> Using measured growth rates\n");
  } else {
    log(">> Computing growth rates using linear relation ...\n");
    compute_growth_rates();
  }
  
  // Create starting bacteria
  bac_counter = 0;
  log(">> Creating bacteria ...\n");
  create_starting_bacteria();
  sprintf(strbuf, "Population: %i\n", (int)bacteria_list.size()); log(strbuf);
  compute_strain_counts();
  
  // Create food sources
  log(">> Creating bacteria ...\n");
  create_starting_food();
  sprintf(strbuf, "Food sources: %i\n", (int)food_list.size()); log(strbuf);
  flag_bacteria_at_food_sources();

  // Ready
 
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
    
    sprintf(strbuf, "Cycle %i, %i bacteria, %i sources\n", cycle, (int)bacteria_list.size(), (int)food_list.size());
    log(strbuf);
    
    if (cycle % output_interval == 0) {
      data_output();
    }

  }
  
  raw_time = time(NULL);
  sprintf(strbuf, ">> %s", ctime(&raw_time)); log(strbuf);
  
  deallocate();
  
}
