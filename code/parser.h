#ifndef __CONFIG_PARSER_H__
#define __CONFIG_PARSER_H__

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define E_FAILED        -1
#define E_SUCCESS        0
#define T_INT            0
#define T_FLOAT          1
#define T_DOUBLE         2

struct p_block {
	char name[32];
	short utype;
	int *ivalue;
	float *fvalue;
	double *dvalue;
};
typedef struct p_block p_block_t;

int config_parser( const char *filename, const int count, p_block_t *a );

#endif