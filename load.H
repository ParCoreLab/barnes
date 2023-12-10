#ifndef _LOAD_H_
#define _LOAD_H_

template<bool is_preprocessing> 
void maketree(long ProcessId);

template<bool is_preprocessing>
cellptr InitCell(cellptr parent, long ProcessId);
leafptr InitLeaf(cellptr parent, long ProcessId);
void printtree(nodeptr n);

template<bool is_preprocessing> 
nodeptr loadtree(bodyptr p, cellptr root, long ProcessId);

bool intcoord(long xp[NDIM], vector rp);
long subindex(long x[NDIM], long l);

template<bool is_preprocessing>
void hackcofm(long ProcessId);

template<bool is_preprocessing>
cellptr SubdivideLeaf(leafptr le, cellptr parent, long l, long ProcessId);

template<bool is_preprocessing>
cellptr makecell(long ProcessId);
leafptr makeleaf(long ProcessId);


#endif
