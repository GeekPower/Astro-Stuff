/* 
 * File:   newmain.c
 * Author: Marian Ionita
 *
 * Created on April 16, 2011, 1:08 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <glib-2.0/glib.h>



typedef struct point{
    int id;
    float x, y;
}point;

#define sqr(x) ((x)*(x))

float max_dist, min_dist, mid_dist;
GSList *field;
point *p;


float max_distance(float distance1, float distance2, float distance3) {
    if (distance1 > distance2 && distance1 > distance3) return distance1;
    if (distance2 > distance3 && distance2 > distance1) return distance2;
    if (distance3 > distance1 && distance3 > distance2) return distance3;
    return -1;
}

float min_distance(float distance1, float distance2, float distance3) {
    if (distance1 < distance2 && distance1 < distance3) return distance1;
    if (distance2 < distance3 && distance2 < distance1) return distance2;
    if (distance3 < distance1 && distance3 < distance2) return distance3;
    return -1;
}

float mid_distance(float distance1, float distance2, float distance3) {
    if (distance1 > distance2 && distance3 > distance1) return distance1;
    if (distance1 > distance2 && distance3 > distance2) return distance2;
    if (distance1 > distance3 && distance2 > distance3) return distance3;
    return -1;
}

static inline float the_distance(GSList *na, GSList *nb)
{
    float x1 = ((point *) na->data)->x;
    float x2 = ((point *) nb->data)->x;
    float y1 = ((point *) na->data)->y;
    float y2 = ((point *) nb->data)->y;

    return sqrtf(sqr(x2 - x1) + sqr(y2 - y1));
}
static char *fname1 = "starf1.list",
            *fname2 = "test.list";

int main(int argc, char** argv) {
    char *buf;
    FILE *in, *out;
    int r;
    float x, y;

    float distance1, distance2, distance3;

    if (argc > 1)
        fname1 = argv[1];

    if (argc > 2)
        fname2 = argv[2];

    if ((in = fopen(fname1, "r")) == NULL)
        return EXIT_FAILURE;

    if ((out = fopen(fname2, "w")) == NULL)
        return EXIT_FAILURE;

    buf = malloc(80 * sizeof(char));

    int i = 0;
    while (fgets(buf, 80, in)) {
        if (!isdigit(buf[0]))
            continue;

        r = sscanf(buf, "%f %f", &x, &y);
        if (r != 2) {
            fprintf(stderr, "sscanf %d", r);
            continue;
        }

        p = malloc(sizeof(point));
        p->x = x;
        p->y = y;
        p->id = i++;

        field = g_slist_prepend(field, p);
    }
    field = g_slist_reverse(field);
    GSList *na, *nb, *nc;
    na = field->next;
 
    while(na != NULL) {
        nb = field->next->next;
        while (nb != NULL) {
            nc = field ->next->next->next;
            while (nc != NULL) {
                distance1 = the_distance(na, nb);
                distance2 = the_distance(nb, nc);
                distance3 = the_distance(na, nc);

                max_dist = max_distance(distance1, distance2, distance3);
                /*min_dist = min_distance((point *)na->data, (point *)nb->data, (point *)nc->data);
                mid_dist = mid_distance((point *)na->data, (point *)nb->data, (point *)nc->data);
                fprintf(out, "%d %d %d %f %f", ((point *)na->data)->id, ((point *)nb->data)->id,
                                    ((point *)nc->data)->id, (float)(mid_dist / max_dist), (float)(min_dist / max_dist));*/

                nc = nc->next;
            }
            nb = nb->next;
        }

        na= na->next;
    }

    fclose(in);
    fclose(out);
    return (EXIT_SUCCESS);
}

