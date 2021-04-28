

#if _EXTRAE_
#include "extrae_user_events.h"
// Extrae Constants
#define  PROGRAM    1000
#define  END        0
#define  SERIAL     1
#define  NWRCOST    2
#define  NWLCOST    3
#define  NWALIGN    4
#define  REVERSE    5
#endif


/*
 * The recursive step in the Needleman-Wunsch algorithm involves
 * choosing the cheapest operation out of three possibilities.  We
 * store the costs of the three operations into an array and use the
 * editop enum to index it.
 */
#define Del 0
#define Sub 1
#define Ins 2



/*
 * levenshtein is the common cost scheme for edit distance: matches cost
 * nothing, while mismatches, insertions and deletions cost one each.
 */
static int
levenshtein(char a, char b)
{
	return a != b;
}

/*
 * nwmin returns the cheapest edit operation out of three possibilities
 * when the "current" characters are a and b, the cost scheme is f,
 * and the base costs of the Del, Sub, and Ins operations are recorded
 * in the cost array in that order.  The cost array is modified by
 * adding the edit costs for a and b to the appropriate cells.
 */
static unsigned int
nwmin(int cost[3], char a, char b)
{
	unsigned int i;
	cost[Del] += levenshtein(a, 0);
	cost[Sub] += levenshtein(a, b);
	cost[Ins] += levenshtein(0, b);
	i = cost[Del] < cost[Sub] ? Del : Sub;
	i = cost[i] < cost[Ins] ? i : Ins;
	return i;
}

static void nwlcost(int *s, const char *a, unsigned int m, const char *b, unsigned int n);
static void nwrcost(int *s, const char *a, unsigned int m, const char *b, unsigned int n);
