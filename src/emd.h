#ifndef _EMD_H
#define _EMD_H
/*
    emd.h

    Last update: 3/24/98

    An implementation of the Earth Movers Distance.
    Based of the solution for the Transportation problem as described in
    "Introduction to Mathematical Programming" by F. S. Hillier and
    G. J. Lieberman, McGraw-Hill, 1990.

    Copyright (C) 1998 Yossi Rubner
    Computer Science Department, Stanford University
    E-Mail: rubner@cs.stanford.edu   URL: http://vision.stanford.edu/~rubner
*/

#ifdef __cplusplus
extern "C"
{
#endif

#include <errno.h>
#include <time.h>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

/* DEFINITIONS */
#define MAX_SIG_SIZE   1000
#define MAX_ITERATIONS 1000
#define EMDINFINITY    1e20
#define EPSILON        1e-12

typedef float coord_t;

void readData ( char * inputFile, coord_t *** xyData, int * dataLen );
void *newmem(int n, int size);
void *newmem2(int number, int size, char *file, int line);

/*****************************************************************************/
/* feature_t SHOULD BE MODIFIED BY THE USER TO REFLECT THE FEATURE TYPE      */
//typedef coord_t feature_t;
typedef double feature_t;
/*****************************************************************************/


typedef struct
{
  int n;                /* Number of features in the signature */
  feature_t *Features;  /* Pointer to the features vector */
  double *Weights;      /* Pointer to the weights of the features */
} signature_t;


typedef struct
{
  int from;             /* Feature number in signature 1 */
  int to;               /* Feature number in signature 2 */
  double amount;        /* Amount of flow from "from" to "to" */
} flow_t;



double emd(signature_t *Signature1,
           signature_t *Signature2,
           double (*func)(feature_t *, feature_t *),
           flow_t *Flow, int *FlowSize);

#ifdef __cplusplus
}
#endif

#endif
