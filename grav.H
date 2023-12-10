#ifndef _GRAV_H_
#define _GRAV_H_

template<bool is_preprocessing>
void hackgrav(bodyptr p, long ProcessId);
void gravsub(register nodeptr p, long ProcessId);
template<bool is_preprocessing>
void hackwalk(long ProcessId);
template<bool is_preprocessing>
void walksub(nodeptr n, real dsq, long ProcessId);
template<bool is_preprocessing>
bool subdivp(register nodeptr p, real dsq, long ProcessId);

#endif
