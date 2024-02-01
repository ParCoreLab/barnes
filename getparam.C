/*************************************************************************/
/*                                                                       */
/*  Copyright (c) 1994 Stanford University                               */
/*                                                                       */
/*  All rights reserved.                                                 */
/*                                                                       */
/*  Permission is given to use, copy, and modify this software for any   */
/*  non-commercial purpose as long as this copyright notice is not       */
/*  removed.  All other uses, including redistribution in whole or in    */
/*  part, are forbidden without prior written permission.                */
/*                                                                       */
/*  This software is provided with absolutely no warranty and no         */
/*  support.                                                             */
/*                                                                       */
/*************************************************************************/

/*
 * GETPARAM.C:
 */
EXTERN_ENV
#include <iostream>
#define global extern

#include "stdinc.h"
//#include <stdio.h>

local string *defaults = NULL;        /* vector of "name=value" strings */
local string inputarray[12];

/*
 * INITPARAM: ignore arg vector, remember defaults.
 */

void initparam(string *defv)
{
   char buf[128];
   long leng;
   defaults = defv;
   for(int i = 0; i < 12; i++) {
	fgets(buf, 128, stdin);
	//std::cerr << "read buf: " << buf << "\n";
	size_t n = strlen(buf);
	if (n > 0 && buf[n - 1] == '\n') {
        	buf[n - 1] = '\0';
   	}
	leng = strlen(buf) + 1;
	int idx = i;
	if(idx > 0 && idx < 3)
		idx++;
	else if (idx == 3)
		idx = 1;
	if(leng > 1) {
		inputarray[idx] = strcpy(static_cast<char*>(malloc(leng)), buf);
	} else {
		inputarray[idx] = strcpy(static_cast<char*>(malloc(1)), "");
	}
	memset( (void *) buf, '\0', 128);
   }
#if 0
    std::cerr << "inputarray:\n";
    for(int i = 0; i < 12; i++)
	std::cerr << i << " " << inputarray[i] << "\n";
#endif
}

/*
 * GETPARAM: export version prompts user for value.
 */

string getparam(string name)
{
   long i, leng = 0;
   string def;
   char buf[128];

   //if (defaults == NULL)
      //error("getparam: called before initparam\n");
   i = scanbind(defaults, name);
   //std::cerr << "in getparam, idx " << i << "\n";
   if (i < 0)
      error("getparam: %s unknown\n", name);
   def = extrvalue(defaults[i]);
   //if(inputarray[i] != NULL)
   leng = strlen(inputarray[i]) + 1;
   if (leng > 1) {
   	return inputarray[i];
   } else {
	return (def);
   }
#if 0
   def = extrvalue(defaults[i]);
   fgets(buf, 128, stdin);
   size_t n = strlen(buf);

   if (n > 0 && buf[n - 1] == '\n') {
        buf[n - 1] = '\0';
   } 
   leng = strlen(buf) + 1;
   if (leng > 1) {
      return (strcpy(static_cast<char*>(malloc(leng)), buf));
   }
   else {
      return (def);
   }
#endif
}

/*
 * GETIPARAM, ..., GETDPARAM: get long, long, bool, or double parameters.
 */

long getiparam(string name)
{
   string val;

   for (val = ""; *val == '\0';) {
      val = getparam(name);
   }
   return (atoi(val));
}

long getlparam(string name)
{
   string val;

   for (val = ""; *val == '\0'; )
      val = getparam(name);
   return (atol(val));
}

bool getbparam(string name)
{
   string val;

   for (val = ""; *val == '\0'; )
      val = getparam(name);
   if (strchr("tTyY1", *val) != NULL) {
      return (TRUE);
   }
   if (strchr("fFnN0", *val) != NULL) {
      return (FALSE);
   }
   error("getbparam: %s=%s not bool\n", name, val);
}

double getdparam(string name)
{
   string val;

   for (val = ""; *val == '\0'; ) {
      val = getparam(name);
   }
   return (atof(val));
}



/*
 * SCANBIND: scan binding vector for name, return index.
 */

long scanbind(string bvec[], string name)
{
   long i;

   for (i = 0; bvec[i] != NULL; i++)
      if (matchname(bvec[i], name))
	 return (i);
   return (-1);
}

/*
 * MATCHNAME: determine if "name=value" matches "name".
 */

bool matchname(string bind, string name)
{
   char *bp, *np;

   bp = bind;
   np = name;
   while (*bp == *np) {
      bp++;
      np++;
   }
   return (*bp == '=' && *np == '\0');
}

/*
 * EXTRVALUE: extract value from name=value string.
 */

string extrvalue(string arg)
{
   char *ap;

   ap = (char *) arg;
   while (*ap != '\0')
      if (*ap++ == '=')
	 return ((string) ap);
   return (NULL);
}

