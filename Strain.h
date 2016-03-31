/*==============\\
|| Strain Class ||
\\==============*/
#ifndef STRAIN_H
#define STRAIN_H
#include <string>
#include "utils.h"

class Strain {

  public:
  
  // Properties of the strain
  int ID;                 // Unique ID; use index in the antagonist table if antagonisms are used
  std::string name;       // Name of the strain
  double growth_rate;     // Growth chance per iteration (if at food source)
  double death_rate;      // Death chance per iteration
  int sender_degree;      // Sender degree, i.e., number of strains this strain antagonizes
  int receiver_degree;    // Receiver degree, i.e., number of strains this strain is antagonized by
  
  // Constructor
  Strain (int p_ID, std::string p_name, double p_growth_rate, double p_death_rate, int p_sender_degree, int p_receiver_degree) {
    ID = p_ID;
    name = p_name;
    growth_rate = p_growth_rate;
    death_rate = p_death_rate;
    sender_degree = p_sender_degree;
    receiver_degree = p_receiver_degree;
  }

  
};
  
#endif // STRAIN_H
