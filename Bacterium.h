/*=================\\
|| Bacterium Class ||
\\=================*/
#ifndef BACTERIUM_H
#define BACTERIUM_H
#include <stdio.h>
#include "FoodSource.h"
#include "Strain.h"
#include "utils.h"

class Strain;
class FoodSource;

class Bacterium {

  public:
  
  // Properties of the bacterium
  int ID;           // Unique sequential ID
  Strain* strain;   // Strain of the bacterium
  int x;            // x-position
  int y;            // y-position
  bool isAlive;     // Is the bacterium still alive?
  bool isAtFoodSource;    // Is the bacterium at a food source?
  
  // Constructor
  Bacterium (int p_ID, Strain* p_strain, int p_x, int p_y) {
    ID = p_ID;
    strain = p_strain;
    x = p_x;
    y = p_y;
    isAlive = true;
    isAtFoodSource = false;
  }
  
  // Moves the bacterium randomly to one of the eight neighboring positions,
  // labelled as follows:
  // 3 2 1 
  // 4 . 0
  // 5 6 7
  void move () {
    switch (randint(8)) {
      case 0: x++;      break; 
      case 1: x++; y++; break;
      case 2:      y++; break;
      case 3: x--; y++; break;
      case 4: x--;      break;
      case 5: x--; y--; break;
      case 6:      y--; break;
      case 7: x++; y--; break;
      default: printf("ERROR in move(): bad direction!"); exit(1);
    }
  }
  
};
  
#endif // BACTERIUM_H
