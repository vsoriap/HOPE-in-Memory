#include <sys/types.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <err.h>
#include "hirschberg.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

double wall_time () {
        struct timeval time;
        gettimeofday(&time, NULL);
        return ((double)time.tv_sec * (double)1e6 + (double)time.tv_usec);
}



char        *hirschberg(const char *, const char*);
static void     hirschberg_recursive(char *, const char *, unsigned int, const char *, unsigned int, char right, int depth);
static char    *nwalign(char *, const char *, unsigned int, const char *, unsigned int, char right);
static void    nwlcost(int *, const char *, unsigned int, const char *, unsigned int);
static void    nwrcost(int *, const char *, unsigned int, const char *, unsigned int);
static unsigned int    nwmin(int[3], char, char);

/*
 * hirschberg calculates the global alignment of a and b with the
 * cost scheme f using Hirschberg's algorithm.  It returns the string
 * of the edit operations required to change a into b:
 * 
 *    +    insertion
 *    -    deletion
 *    =    match
 *    !    substitution
 * 
 */
char *
hirschberg(const char *a, const char *b){
    char *c;
    unsigned int m = strlen(a);
    unsigned int n = strlen(b);
    
    /*
     * The alignment string of a and b is at most as long as the
     * concatenation of the two (delete all of a + insert all of b).
     */
    c = malloc(m+n+1);
    if (c == NULL)
        return NULL;
    if (m > n) {
        /*
         * hirschberg_recursive assumes that the first string
         * passed to it is the shorter one, so if a is not, we
         * flip it and b.  The resulting alignment is otherwise
         * equivalent to the non-flipped one, with the exception
         * that insertion and deletion operators need to be
         * flipped.
         */
        #pragma omp parallel 
        {
            #pragma omp single
            {
                #pragma omp task
                hirschberg_recursive(c, b, n, a, m, 1, 1);
            }
        }
        #pragma omp parallel for simd
        for (int i = 0; i < (m+n+1); i++)
            switch (c[i]) {
            case '+': c[i] = '-'; break;
            case '-': c[i] = '+'; break;
        }
        return c;
    }
    #pragma omp parallel 
    {   
        #pragma omp single
        {   
            #pragma omp task
            hirschberg_recursive(c, a, m, b, n, 1, 1);
        }
    }
    return c;
}

/*
 * hirschberg_recursive is the recursive part of Hirschberg's algorithm.
 * The arguments are the same as hirschberg, with the exception that the
 * length m of a and the length n of b are now explicitly passed, and c
 * is a pointer to the buffer where the alignment string is to be
 * written.  hirschberg_recursive returns a pointer to the null byte
 * written after the last alignment character.
 */
void
hirschberg_recursive(char *c, const char *a, unsigned int m, const char *b, unsigned int n, char right, int depth){
    if (n > 1) {
        unsigned int mmid, nmid;
        int *lcost, *rcost;
        int min;
        lcost = malloc(sizeof(int)*(m+1));
        rcost = malloc(sizeof(int)*(m+1));
        nmid = n / 2;

        if (depth < 10000){
            #pragma omp task shared(a,b,lcost) 
            nwlcost(lcost, a, m, b, nmid);

            #pragma omp task shared(a,b,rcost) 
            nwrcost(rcost, a, m, b+nmid, n-nmid);

            #pragma omp taskwait
            mmid = 0;
            min = lcost [0] + rcost[0];

            #pragma omp parallel
            {
                int index_local = mmid;
                int min_local = min;
                #pragma omp for simd nowait
                for (int i = 1; i <= m; i++) {
                    int tmp = lcost[i] + rcost[i];
                    if (tmp < min_local) {
                        min_local = tmp;
                        index_local = i;
                    }
                }

                #pragma omp critical 
                {
                    if (min_local < min) {
                        min = min_local;
                        mmid = index_local;
                    }
                }
            }

            free(lcost);
            free(rcost);

            #pragma omp task
            hirschberg_recursive(c, a, mmid, b, nmid, 0, depth+1);
        
            #pragma omp task
            hirschberg_recursive(c+mmid+nmid, a+mmid, m-mmid, b+nmid, n-nmid, right, depth+1);

            #pragma omp taskwait

        } else {
            #pragma omp task shared(a,b,lcost) 
            nwlcost(lcost, a, m, b, nmid);

            #pragma omp task shared(a,b,rcost) 
            nwrcost(rcost, a, m, b+nmid, n-nmid);

            #pragma omp taskwait 
            mmid = 0;
            min = lcost [0] + rcost[0];
            
            #pragma omp parallel
            {   
                int index_local = mmid;
                int min_local = min;
                #pragma omp for simd nowait
                for (int i = 1; i <= m; i++) {
                    int tmp = lcost[i] + rcost[i];
                    if (tmp < min_local) {
                        min_local = tmp;
                        index_local = i;
                    }
                }
                
                #pragma omp critical 
                {
                    if (min_local < min) {
                        min = min_local;
                        mmid = index_local;
                    }
                }
            }

            free(lcost);
            free(rcost);

            hirschberg_recursive(c, a, mmid, b, nmid, 0, depth+1);

            hirschberg_recursive(c+mmid+nmid, a+mmid, m-mmid, b+nmid, n-nmid, right, depth+1);
        }
    } else {
        c = nwalign(c+m+n+right-1, a, m, b, n, right);
    }
}

/*
 * nwalign computes the Needleman-Wunsch alignment of a and b using the
 * cost scheme f and writes it into the buffer c.  It returns a pointer
 * to the null byte written after the last character in the alignment
 * string.
 * 
 * This function uses O(mn) space.  hirschberg_recursive guarantees its
 * own O(m) space usage by only calling this when n <= 1.
 */
static char *
nwalign(char *c, const char *a, unsigned int m, const char *b, unsigned int n, char right)
{
#if _EXTRAE_
    Extrae_event (PROGRAM, END);
    Extrae_event (PROGRAM, NWALIGN);
#endif
    unsigned int i, j;
    int *s;
    s = malloc(sizeof(int)*((m+1)*(n+1)));
    
    s[0] = 0;
    for (i = 1; i <= m; i++)
        s[m*i] = s[m*(i-1)] + levenshtein(a[i-1], 0);
    for (j = 1; j <= n; j++)
        s[j] = s[j-1] + levenshtein(0, b[j-1]);
    for (j = 1; j <= n; j++)
        for (i = 1; i <= m; i++) {
            int cost[3] = { s[m*(i-1)+j], s[m*(i-1)+j-1], s[m*(i)+j-1] };
            s[m*i+j] = cost[nwmin(cost, a[i-1], b[j-1])];
        }
    i = m;
    j = n;
    if(right){
        *c-- ='\0';
    }
    while (i > 0 && j > 0) {
        int cost[3] = { s[m*(i-1)+j], s[m*(i-1)+j-1], s[m*(i)+j-1] };
        switch (nwmin(cost, a[i-1], b[j-1])) {
        case Del:
            *c-- = '-';
            i--;
            break;
        case Sub:
            *c-- = a[i-1] == b[j-1] ? '=' : '!';
            *c-- = a[i-1] == b[j-1] ? '=' : '!';
            i--;
            j--;
            break;
        case Ins:
            *c-- = '+';
            j--;
            break;
        }
    }
    for (; i > 0; i--)
        *c-- = '-';
    for (; j > 0; j--)
        *c-- = '+';
#if _EXTRAE_
    Extrae_event (PROGRAM, END);
    Extrae_event (PROGRAM, SERIAL);
#endif

    return c;
}

/*
 * nwlcost stores the last column of the Needleman-Wunsch alignment
 * cost matrix of a and b into s.
 */
static void
nwlcost(int *s, const char *a, unsigned int m, const char *b, unsigned int n)
{
#if _EXTRAE_
    Extrae_event (PROGRAM, END);
    Extrae_event (PROGRAM, NWLCOST);
#endif

    int i, j, k, l, p, q;
    int ss, tmp;
    int *vecs[3];
    vecs[0] = malloc(sizeof(int)*(m+1));
    vecs[1] = malloc(sizeof(int)*(m+1));
    vecs[2] = malloc(sizeof(int)*(m+1));
    int *act, *prev, *pprev; 
    
    vecs[0][0] = 0;
    prev = vecs[0];
    act = vecs[1];
    i=1;j=0;p=0;q=0;k=1;



    while(i <= n && j <= m) {
        p = i;

        if (q == 0) { // Do diagonals until reaching the last row
            // Compute the number of iterations to be done
            // for(int l = 0; p >= 0 && q <= m; p--, q++, l++){
            int maxII = MIN(p, m-q);
            // Compute first cell that corresponds with first column
            act[0] = prev[0] + 1;
            p--; q++; l=1;

            #pragma omp simd
            for(l = 1; l < maxII; l++){
                int cost[3], tmp; 
                cost[0] = prev[l-1] + 1; 
                cost[1] = pprev[l-1] + (a[q-1] != b[p-1]); 
                cost[2] = prev[l] + 1; 
                tmp = cost[0] < cost[1] ? cost[0] : cost[1];
                tmp = tmp < cost[2] ? tmp : cost[2];
                act[l] = tmp;
                p--;
                q++;
            }

            //Compute last row
            if (maxII > 0){
                int cost[3], tmp;
                if (p == 0){ 
                    act[l] = prev[l-1] + 1;
                } else {
                    cost[0] = prev[l-1] + 1;
                    cost[1] = pprev[l-1] + (a[q-1] != b[p-1]);
                    cost[2] = prev[l] + 1;
                    tmp = cost[0] < cost[1] ? cost[0] : cost[1];
                    tmp = tmp < cost[2] ? tmp : cost[2];
                    act[l] = tmp;
                }
            }
        } else if (q == 1){ //Do transition diagonal
            //Compute the number of iterations to be done
            // for(int l = 0; p >= 0 && q <= m; p--, q++, l++){
            int maxII = MIN(p, m-q);
            
            #pragma omp simd
            for(l = 0; l < maxII; l++){
                int cost[3], tmp;
                cost[0] = prev[l] + 1;
                cost[1] = pprev[l] + (a[q-1] != b[p-1]);
                cost[2] = prev[l+1] + 1;
                tmp = cost[0] < cost[1] ? cost[0] : cost[1];
                tmp = tmp < cost[2] ? tmp : cost[2];
                act[l] = tmp; 
                p--;
                q++;
            }

            //Compute last row
            int cost[3], tmp;
            if (p == 0){
                act[l] = prev[l] + 1;
            } else {
                cost[0] = prev[l] + 1;
                cost[1] = pprev[l] + (a[q-1] != b[p-1]);
                cost[2] = prev[l+1] + 1;
                tmp = cost[0] < cost[1] ? cost[0] : cost[1];
                tmp = tmp < cost[2] ? tmp : cost[2];
                act[l] = tmp;
            }
        } else { // Do the rest diagonals
            //Compute the number of iterations to be done
            // for(int l = 0; p >= 0 && q <= m; p--, q++, l++){
            int maxII = MIN(p, m-q);

            #pragma omp simd
            for(l = 0; l < maxII; l++){
                int cost[3], tmp;
                cost[0] = prev[l] + 1; 
                cost[1] = pprev[l+1] + (a[q-1] != b[p-1]); 
                cost[2] = prev[l+1] + 1; 
                tmp = cost[0] < cost[1] ? cost[0] : cost[1];
                tmp = tmp < cost[2] ? tmp : cost[2];
                act[l] = tmp;
                p--;
                q++;
            }

            //Compute last row
            int cost[3], tmp;
            if (p == 0){
                act[l] = prev[l] + 1;
            } else {
                cost[0] = prev[l] + 1;
                cost[1] = pprev[l+1] + (a[q-1] != b[p-1]);
                cost[2] = prev[l+1] + 1;
                tmp = cost[0] < cost[1] ? cost[0] : cost[1];
                tmp = tmp < cost[2] ? tmp : cost[2];
                act[l] = tmp;
            }                                                                                                      
        }

        k++;
        pprev = prev;
        prev = act;
        act = vecs[k%3];

        if(i <= n-1) {
            i++;
            q = 0;
        } else if(j <= m-1) {
            s[j] = prev[0];
            j++;
            q = j;
        } else {
            s[j] = prev[0]; 
            break;
        }
    }
    
    free(vecs[0]);
    free(vecs[1]);
    free(vecs[2]);

#if _EXTRAE_
    Extrae_event (PROGRAM, END);
    Extrae_event (PROGRAM, SERIAL);
#endif

}

/*
 * nwrcost computes the reverse Needleman-Wunsch alignment of a and b,
 * that is, matching their suffixes rather than prefixes.  The last
 * column of this alignment cost matrix is stored into s.
 */
static void
nwrcost(int *s, const char *a, unsigned int m, const char *b, unsigned int n)
{
#if _EXTRAE_
    Extrae_event (PROGRAM, END);
    Extrae_event (PROGRAM, NWRCOST);
#endif

    int i, j, k, l, p, q;
    int ss, tmp;
    int *vecs[3];
    vecs[0] = malloc(sizeof(int)*(m+1));
    vecs[1] = malloc(sizeof(int)*(m+1));
    vecs[2] = malloc(sizeof(int)*(m+1));
    int *act, *prev, *pprev;

    vecs[0][0] = 0;
    prev = vecs[0];
    act = vecs[1];
    i=n-1;j=m;p=n-1;q=m;k=1;

    while(i >= 0 && j >= 0 ) {
        p = i;

        if (q == m) { // Do diagonals until reaching the first row
            // Compute the number of iterations to be done
            // for(int l = 0; p >= 0 && q <= m; p--, q++, l++){
            int maxII = MIN(n-p, q);
            // Compute first cell that corresponds with last column
            act[0] = prev[0] + 1;
            p++; q--; l=1;

            #pragma omp simd
            for(l = 1; l < maxII; l++){
                int cost[3], tmp;
                cost[0] = prev[l-1] + 1;
                cost[1] = pprev[l-1] + (a[q] != b[p]);
                cost[2] = prev[l] + 1;
                tmp = cost[0] < cost[1] ? cost[0] : cost[1];
                tmp = tmp < cost[2] ? tmp : cost[2];
                act[l] = tmp;
                p++;
                q--;
            }

            //Compute last row
            if (maxII > 0){
                int cost[3], tmp;
                if (p == n){
                    act[l] = prev[l-1] + 1;
                } else {
                    cost[0] = prev[l-1] + 1;
                    cost[1] = pprev[l-1] + (a[q] != b[p]);
                    cost[2] = prev[l] + 1;
                    tmp = cost[0] < cost[1] ? cost[0] : cost[1];
                    tmp = tmp < cost[2] ? tmp : cost[2];
                    act[l] = tmp;
                }
            }
        } else if (q == m-1){ //Do transition diagonal
            //Compute the number of iterations to be done
            // for(int l = 0; p >= 0 && q <= m; p--, q++, l++){
            int maxII = MIN(n-p, q);

            #pragma omp simd
            for(l = 0; l < maxII; l++){
                int cost[3], tmp;
                cost[0] = prev[l] + 1;
                cost[1] = pprev[l] + (a[q] != b[p]);
                cost[2] = prev[l+1] + 1;
                tmp = cost[0] < cost[1] ? cost[0] : cost[1];
                tmp = tmp < cost[2] ? tmp : cost[2];
                act[l] = tmp;
                p++;
                q--;
            }

            //Compute last row         
            int cost[3], tmp;
            if (p == n){
                act[l] = prev[l] + 1;
            } else {
                cost[0] = prev[l] + 1;
                cost[1] = pprev[l] + (a[q] != b[p]);
                cost[2] = prev[l+1] + 1;
                tmp = cost[0] < cost[1] ? cost[0] : cost[1];
                tmp = tmp < cost[2] ? tmp : cost[2];
                act[l] = tmp;
            }
        } else { // Do the rest diagonals
            //Compute the number of iterations to be done
            // for(int l = 0; p >= 0 && q <= m; p--, q++, l++){
            int maxII = MIN(n-p, q);

            #pragma omp simd
            for(l = 0; l < maxII; l++){
                int cost[3], tmp;
                cost[0] = prev[l] + 1;
                cost[1] = pprev[l+1] + (a[q] != b[p]); 
                cost[2] = prev[l+1] + 1;
                tmp = cost[0] < cost[1] ? cost[0] : cost[1];
                tmp = tmp < cost[2] ? tmp : cost[2];
                act[l] = tmp;
                p++;
                q--;
            }

            int cost[3], tmp;
            if (p == n){
                act[l] = prev[l] + 1;
            } else {
                cost[0] = prev[l] + 1;
                cost[1] = pprev[l+1] + (a[q] != b[p]);
                cost[2] = prev[l+1] + 1;
                tmp = cost[0] < cost[1] ? cost[0] : cost[1];
                tmp = tmp < cost[2] ? tmp : cost[2];
                act[l] = tmp;

            }
        }


        k++;
        pprev = prev;
        prev = act;
        act = vecs[k%3];

        if(i > 0) {
            i--;
            q = m;
        } else if(j > 0) {
            s[j] = prev[0];
            j--;
            q = j;
        } else {
            s[j] = prev[0];
            break;
        }
    }
        
    free(vecs[0]);
    free(vecs[1]);
    free(vecs[2]);

#if _EXTRAE_
    Extrae_event (PROGRAM, END);
    Extrae_event (PROGRAM, SERIAL);
#endif

}

#include <sys/stat.h>

/*
 * main takes two strings a and b as arguments and prints their global
 * alignment to standard output in three lines: the alignment of a, the
 * edit sequence, and the alignment of b.
 */
int
main(int argc, char *argv[]) {
    char *align, *c;

    const char *file_a=argv[1];
    struct stat st;
    const char *file_b=argv[2];
    char       *line  = NULL;

    stat(file_a, &st);
    unsigned int size_a = st.st_size;
    char            * a = malloc(size_a);
    stat(file_b, &st);
    unsigned int size_b = st.st_size;
    char            * b = malloc(size_b);
       
    FILE* fd_a = fopen(file_a, "r"); 
    a=fgets(a,size_a,fd_a);

    FILE* fd_b = fopen(file_b, "r"); 
    b=fgets(b,size_b,fd_b);

    double t_start = wall_time();

#if _EXTRAE_
    Extrae_event (PROGRAM, SERIAL);
#endif

    align = hirschberg(a, b);

#if _EXTRAE_
    Extrae_event (PROGRAM, END);
#endif

    double t_end = wall_time();

    if (align == NULL)
        err(1, "hirschberg");
    
    for (c = align; *c != '\0'; c++)
        switch (*c) {
        case '-':
            putchar(*a++);
            break;
        case '!':
        case '=':
            putchar(*a++);
            c++;
            break;
        default:
            putchar(' ');
            break;
        }
    putchar('\n');

    for (c = align; *c != '\0'; c++)
        switch (*c) {
        case '-':
        case '+':
            putchar(*c);
            break;
        case '!':
        case '=':
            putchar(*c);
            c++;
            break;
        default:
            putchar(' ');
            break;
        }
    putchar('\n');
    
    for (c = align; *c != '\0'; c++)
        switch (*c) {
        case '+':
            putchar(*b++);
            break;
        case '!':
        case '=':
            putchar(*b++);
            c++;
            break;
        default:
            putchar(' ');
            break;
        }
    putchar('\n');

    printf( "==================== RESULTS ===================== \n" );
    printf( "  Execution time (secs): %f\n", (t_end - t_start)/(double)1e6);
    printf( "================================================== \n" );
    
    return 0;
}
