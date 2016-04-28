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
bool  enable_antagonism;         // Simulate for antangonistic interactions?
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

// GLOBAL VARIABLES

// Current cycle number
int cycle;

// Sequential labels for all bacteria that ever lived
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

// Global debug flag
const bool global_debug = false;

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
  fprintf(log_file, "OK\n");

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
  fprintf(log_file, "OK\n");
  
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
void remove_bacterium(unsigned int idx) {
  char strbuf[128];
  if (global_debug) {
    bacteria_list[idx]->repr(strbuf);
    printf("Removing bacterium idx %i - %s\n", idx, strbuf);
  }
  delete bacteria_list[idx];
  bacteria_list[idx] = bacteria_list.back();
  bacteria_list.pop_back();
}

// =====================================================================

// Adds a new food source at (x,y)
void add_food_source(int x, int y, double radius, int rations) {
  food_list.push_back(new FoodSource(x, y, radius, rations));
}

// Removes the food source residing at 'idx' in the food_list
void remove_food_source(unsigned int idx) {
  if (global_debug) {
    printf("Removing food source idx %i - position (%i,%i)\n", idx, food_list[idx]->x, food_list[idx]->y);
  }
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
    }
  }

}

// =====================================================================

// Creates the starting food sources
void create_starting_food() {
  for (int i = 0; i < start_num_food_sources; i++) {
    add_food_source(randint(1, grid_size), randint(1, grid_size), food_source_radius, food_source_rations);
  }
}

// Spawns a food source by rolling against the spawn rate
void spawn_food_source() {
  if (randreal() <= food_source_spawn_rate) {
    add_food_source(randint(1, grid_size), randint(1, grid_size), food_source_radius, food_source_rations);
    if (global_debug) {
      printf("Spawned food source at %i,%i with %i rations\n", food_list.back()->x, food_list.back()->y, food_list.back()->rations);
    }
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
  
}

// =====================================================================

// Kill bacteria according to their death rate (only flags them as dead)
void kill_bacteria() {

  Bacterium* bac;
  unsigned int i;

  for (i = 0; i < bacteria_list.size(); i++) {
    bac = bacteria_list[i];
    if (bac->isAlive && (all_die || !bac->isAtFoodSource)) {
      if (randreal() <= bac->strain->death_rate) {
        bac->isAlive = false;
        if (global_debug) printf("Bacteria #%i died\n", bac->ID);
      }
    }
  }
  
}

// =====================================================================


// Clean up empty food sources (by removing them from the food sources list)
void cleanup_empty_food_sources() {

  unsigned int idx, i, size;

  idx = 0;
  size = food_list.size();
  for (i = 0; i < size; i++) {
    if (food_list[idx]->rations == 0) {
      remove_food_source(idx);
    } else {
      idx++;
    }    
  }

}

// Clean up dead bacteria (by removing them from the bacteria list)
void cleanup_dead_bacteria() {

  unsigned int idx, i, size;

  idx = 0;
  size = bacteria_list.size();
  for (i = 0; i < size; i++) {
    if (!bacteria_list[idx]->isAlive) {
      remove_bacterium(idx);
    } else {
      idx++;
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

  for (i = 0; i < food_list.size(); i++) {
    food_list[i]->clearBacteriaHere();
  }

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
        if (global_debug) printf("Found bacterium #%i at food source %i (%i,%i)\n", bac->ID, j, fs->x, fs->y);
      }
      
    }
  }

}

// =====================================================================

// Implements interactions that occur at food sources
// Specifically:
// - bacteria antagonize each other (if enabled)
// - if not enough food for remaining bacteria, a lottery is made
// - lucky bacteria eat and perhaps reproduce
void food_sources_interactions() {

  FoodSource* fs = NULL;
  Bacterium* bac = NULL;
  Bacterium* bac2 = NULL;
  unsigned int fs_idx, idx2, size, num_bacs, i, j, num_will_eat;
  char strbuf[256];

  const bool debug = global_debug;
  
  // Loop over all food sources
  for (fs_idx = 0; fs_idx < food_list.size(); fs_idx++) {
    fs = food_list[fs_idx];
   
    // Bacteria count at this food source
    num_bacs = (int)fs->bacteriaHere.size();
    if (debug) printf("\nSource %i at %i,%i with %i rations has %i bacteria\n", fs_idx, fs->x, fs->y, fs->rations, num_bacs);
    if (num_bacs == 0) continue;

    if (debug) {
      for (i = 0; i < num_bacs; i++) {
        fs->bacteriaHere[i]->repr(strbuf);
        printf("%i - %s\n", i, strbuf);
      }
    }

    // Antagonisms
    if (enable_antagonism) {
      
      // For each bacteria, check all other bacteria at this food
      // source. If it antagonizes any, it kills it and either
      // reproduces automatically or gets a "free" chance to do so
      for (i = 0; i < num_bacs; i++) {
        bac = fs->bacteriaHere[i];
        if (bac->isAlive) {

          for (j = 0; j < num_bacs; j++) {
            if (i == j) continue;
            bac2 = fs->bacteriaHere[j];
            if (bac2->isAlive) {
              
              // Check if bac antagonizes bac2
              if (antagonisms[bac->strain->ID][bac2->strain->ID]) {

                // Finish him!
                bac2->isAlive = false;
                
                // Spoils of combat!
                
                // Simplest mode: automatically reproduce
                add_bacterium(bac->strain, bac->x, bac->y);
                
                
                
              }
              
            }
          }

        }
      }

      // Cleanup the carnage at this food source
      // (won't eliminate the bacteria from the global list yet though)
      idx2 = 0;
      size = fs->bacteriaHere.size();
      for (i = 0; i < size; i++) {
        if (!fs->bacteriaHere[idx2]->isAlive) {
          fs->bacteriaHere[idx2] = fs->bacteriaHere.back();
          fs->bacteriaHere.pop_back();
        } else {
          idx2++;
        }
      }
      num_bacs = (int)fs->bacteriaHere.size();

    }

    // Check available rations at this food source
    if (fs->rations < num_bacs) {

      // If not enough food for everyone, randomly select who will eat
      num_will_eat = fs->rations;

      // Partial Fisher-Yates shuffle
      for (i = 0; i < num_will_eat; i++) {
        j = randint(i,num_bacs-1);
        if (debug) {
          printf("i=%i <-> j=%i ", i, j); fflush(stdout);
          printf("| %i\n", fs->bacteriaHere[j]->strain->ID);
        }
        if (j != i) {
          bac = fs->bacteriaHere[j];
          fs->bacteriaHere[j] = fs->bacteriaHere[i];
          fs->bacteriaHere[i] = bac;
        }
      }

    } else {
      
      // If enough food, everyone gets to eat
      num_will_eat = num_bacs;
      
    }
    if (debug) printf("%i bacteria will eat\n", num_will_eat);
    
    // Now feed the bacteria, reducing the food source size, and let
    // the fed bacteria have a chance to reproduce
    if (debug) printf("Bacteria that reproduced:");
    for (i = 0; i < num_will_eat; i++) {
      fs->rations--;
      bac = fs->bacteriaHere[i];
      if (randreal() <= bac->strain->growth_rate) {
        if (debug) printf(" #%i", bac->ID);
        add_bacterium(bac->strain, bac->x, bac->y);
      }
    }
    if (debug) printf("\n");
    
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

// Deallocates dynamic memory
void deallocate() {
  
  int i;
  delete strain_counts;
  delete strains;
  for (i = 0; i < num_strains; i++) {
    delete antagonisms[i];
  }
  delete antagonisms;
  delete log_fname;
  delete output_fname;
  delete strains_fname;
  delete antagonism_table_fname;
  
}

// =====================================================================

// Returns the next data-holding line, skipping comments and blank lines
std::string get_next_line(std::ifstream& file, const char* param_name) {

  char firstchar;
  std::string line;
  unsigned int i;

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
  
  line = get_next_line(params_file, "enable_antagonism");
  enable_antagonism = (atoi(line.c_str()) == 1);
  
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

// Logs a message to the logfile, echoing to screen if requested
void log(const char* message, bool echo_stdout) {
  fputs(message, log_file);
  if (echo_stdout) fputs(message, stdout); 
}

// Wrapper for the above, echoing to stdout if report_to_screen is true
void log(const char* message) {
  log(message, report_to_screen);
}

// #####################################################################
// #####################################################################

int main (int argc, char** argv) {

  char strbuf[256];
  time_t raw_time;

  // Load parameters from the specified file
  if (argc != 2) {
    printf("Error! Must provide a parameters filename as first argument\n");
    exit(EXIT_FAILURE);
  } else {
    raw_time = time(NULL);
    printf(">> %s", ctime(&raw_time));
  }
  params_fname = new char[strlen(argv[1])];
  strcpy(params_fname, argv[1]);
  printf(">> Reading parameters file %s ...\n", params_fname);
  load_parameters(params_fname);
  
  const int report_interval = (num_cycles > 100 ? num_cycles : 100) / 100; 
  
  // Open log and data files
  log_file = fopen(log_fname,"w");
  output_file = fopen(output_fname,"w");
  printf(">> Opened log file %s\n", log_fname);
  raw_time = time(NULL);
  sprintf(strbuf, ">> %s", ctime(&raw_time)); log(strbuf, false);
  sprintf(strbuf, ">> Read parameters file %s ...\n", params_fname); log(strbuf, false);

  // Report parameters
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
  sprintf(strbuf, "enable_antagonism = %s\n", enable_antagonism ? "true" : "false"); log(strbuf);
  sprintf(strbuf, "antagonism_table_fname = %s\n", antagonism_table_fname); log(strbuf);
  sprintf(strbuf, "use_real_growth = %s\n", use_real_growth ? "true" : "false"); log(strbuf);
  sprintf(strbuf, "growth_intercept = %f\n", growth_intercept); log(strbuf);
  sprintf(strbuf, "growth_slope = %f\n", growth_slope); log(strbuf);
  sprintf(strbuf, "food_source_rations = %i\n", food_source_rations); log(strbuf);
  sprintf(strbuf, "food_source_radius = %i\n", food_source_radius); log(strbuf);
  sprintf(strbuf, "food_source_spawn_rate = %f\n", food_source_spawn_rate); log(strbuf);
  sprintf(strbuf, "start_population_per_strain = %i\n", start_population_per_strain); log(strbuf);
  sprintf(strbuf, "start_num_food_sources = %i\n", start_num_food_sources); log(strbuf);
  log(">> Parameters loaded successfully\n");

  // Allocate data arrays
  do_allocations();
  log(">> Data arrays allocated\n");
  
  // Set RNG seed
  reset_seed();
  sprintf(strbuf, ">> RNG seed: %u\n", global_seed); log(strbuf);
  
  // Create strains using data from file
  create_strains();
  
  // Load antagonism matrix (if applicable)
  if (enable_antagonism) {
    log(">> Antagonisms ENABLED\n");
    load_antagonisms();
  } else {
    log(">> Antagonisms DISABLED\n");
  }

  // Compute growth_rates if not using real ones
  if (use_real_growth) {
    log(">> Using measured growth rates\n");
  } else {
    log(">> Computing growth rates from sender degrees ...\n");
    compute_growth_rates();
  }
  
  // Create starting bacteria
  bac_counter = 0;
  log(">> Creating bacteria ...\n");
  create_starting_bacteria();
  sprintf(strbuf, "Population: %i\n", (int)bacteria_list.size()); log(strbuf);
  compute_strain_counts();
  
  // Create food sources
  log(">> Creating food sources ...\n");
  create_starting_food();
  sprintf(strbuf, "Food sources: %i\n", (int)food_list.size()); log(strbuf);
  flag_bacteria_at_food_sources();

  // Ready!
  log(">> READY. Starting simulation ...\n", true);
  if (!report_to_screen) {
    printf("(no further output written to screen)\n");
  }
 
  // Main loop over cycles
  for (cycle = 1; cycle <= num_cycles; cycle++) {

    move_bacteria();

    flag_bacteria_at_food_sources();

    food_sources_interactions();

    cleanup_empty_food_sources();

    kill_bacteria();

    cleanup_dead_bacteria();
    
    spawn_food_source();

    compute_strain_counts();
    
    sprintf(strbuf, "Cycle %i, %i bacteria, %i sources\n", cycle, (int)bacteria_list.size(), (int)food_list.size());
    log(strbuf, false);
    if (report_to_screen) {
      if (cycle % report_interval == 0) {
        printf(">> Cycle %i, %i bacteria, %i sources\n", cycle, (int)bacteria_list.size(), (int)food_list.size());
      }
    }
    
    if (cycle % output_interval == 0) {
      data_output();
    }

  }
  
  raw_time = time(NULL);
  sprintf(strbuf, ">> %s", ctime(&raw_time)); log(strbuf, true);
  
  deallocate();
  
}
