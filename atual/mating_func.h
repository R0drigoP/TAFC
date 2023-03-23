#include "molecula.h"
#include <ctime>

void parent_probability(vector<molecula*> pop, bool *is_parent, int *parent_order);

void generate_children(vector<molecula*> pop, int *parent_order);

void generate_children2(vector<molecula*> pop);
