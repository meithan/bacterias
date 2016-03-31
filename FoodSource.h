/*==================\\
|| FoodSource Class ||
\\==================*/
#ifndef FOODSOURCE_H
#define FOODSOURCE_H
#include <string>
#include "Bacterium.h"
#include "utils.h"

class Bacterium;

class FoodSource {

  public:
  
  int x;            // x-position (integer)
  int y;            // y-position (integer)
  double radius;    // Radius of food source (can be a real)
  double radius2;   // Squared radius
  unsigned int rations;      // Food "rations" left in source
  std::vector<Bacterium*> bacteriaHere;   // Bacteria currently at this food source
  
  // Constructor
  FoodSource (int p_x, int p_y, double p_radius, int p_rations) {
    x = p_x;
    y = p_y;
    radius = p_radius;
    radius2 = radius*radius;
    rations = p_rations;
    bacteriaHere.reserve(10);
  }

  // Returns whether the point (x,y) is "inside" the food source or not
  bool contains(int x1, int y1) {
    return ((x-x1)*(x-x1) + (y-y1)*(y-y1) <= radius2);
  }

  // Clears the list of bacteria at this food source
  void clearBacteriaHere() {
    bacteriaHere.clear();
    bacteriaHere.reserve(64);   // for faster future use
  }
 
};
  
#endif // FOODSOURCE_H
