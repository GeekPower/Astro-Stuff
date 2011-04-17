/* A would-be implementation of the algorithm in

   A Pattern-Matching Algorithm for Two-Dimensional Coordinate Lists
   Edward J. Groth, 1986
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include <glib.h>

#define N_OBJ 40
#define POS_EPSILON 0.001
#define POS_THRESHOLD 3.0		     /* times POS_TOL */
#define RATIO_LIMIT 10.0

#ifndef MAX
#define MAX(a, b) ((a) < (b) ? (b) : (a))
#endif

static double width = 0.0, height = 0.0;
static double pos_eps, pos_thresh;

#define sqr(x) ((x) * (x))

struct point {
	double x;
	double y;

	/* holds distances to previous points in array */
	double *dists;
};

struct triangle {
	/* from shortest to longest: 12, 23, 31 */
	int v[3];			  /* vertices */
	double logP;			  /* logarithm of perimeter */
	unsigned int clockwise:1;	  /* vertex order */
	double R;			  /* longest/shorter */
	double Rtolsq;			  /* ratio tolerance, squared */
	double C;			  /* cosine of angle at v1 */
	double C2;
	double Ctolsq;			  /* cosine tolerance, squared */
};

static struct point *point_new(double x, double y)
{
	struct point *pt = calloc(1, sizeof(struct point));

	pt->x = x;
	pt->y = y;

	return pt;
}

static inline double get_distance(struct point pts[], int v1, int v2)
{
	g_assert(v1 != v2);
	return (v1 > v2) ? pts[v1].dists[v2] : pts[v2].dists[v1];

}

static inline struct triangle *triangle_new (struct point pts[], int v1, int v2, int v3)
{
	struct triangle *t = calloc(1, sizeof(struct triangle));
	int v[3], vv[3];
	double r1, r2, r3;

	vv[0] = v1;
	vv[1] = v2;
	vv[2] = v3;

	/* FIXME: permute until ordered righ :-) */
	for (v[0] = 0; v[0] < 3; v[0]++) {
		for (v[1] = 0; v[1] < 3; v[1]++) {
			if (v[1] == v[0])
				continue;

			for (v[2] = 0; v[2] < 3; v[2]++) {
				if (v[2] == v[1] || v[2] == v[0])
					continue;

				/* NB: r2 latura mica, r1 latura medie, r3 latura mare */
				r2 = get_distance(pts, vv[v[0]], vv[v[1]]);
				r1 = get_distance(pts, vv[v[1]], vv[v[2]]);
				r3 = get_distance(pts, vv[v[2]], vv[v[0]]);

				if (r2 < r1 && r1 < r3)
					goto found;
			}
		}
	}

found:
	g_assert(r2 < r1 && r1 < r3 && r2 < r3);

	t->v[0] = vv[v[0]];
	t->v[1] = vv[v[1]];
	t->v[2] = vv[v[2]];

	t->R = r3 / r2;

	/* law of cosines: cos beta = (a^2 + c^2 - b^2) / (2 * a * c) */
	t->C = (sqr(r2) + sqr(r3) - sqr(r1)) / (2 * r2 * r3);

	t->C2 = ( (pts[t->v[2]].x - pts[t->v[0]].x) * (pts[t->v[1]].x - pts[t->v[0]].x) +
		  (pts[t->v[2]].y - pts[t->v[0]].y) * (pts[t->v[1]].y - pts[t->v[0]].y) ) / (r2 * r3);

	double cterm = 1.0 / sqr(r3) - t->C / (r3 * r1) + 1.0 / sqr(r1);

	t->Rtolsq = 2 * sqr(t->R) * sqr(pos_eps) * cterm;

	double sinsq = 1 - sqr(t->C);

	t->Ctolsq = 2 * sinsq * sqr(pos_eps) * cterm + 3 * sqr(t->C) * sqr(sqr(pos_eps)) * cterm;

	/* FIXME: if det. is positive, we're clockwise. Rewrite this.... */
	t->clockwise = (((pts[t->v[1]].x - pts[t->v[0]].x) * (pts[t->v[2]].y - pts[t->v[0]].y) -
			 (pts[t->v[1]].y - pts[t->v[0]].y) * (pts[t->v[2]].x - pts[t->v[1]].x)) > 0);

	printf("%d (%.2lf, %.2lf) - %d (%.2lf, %.2lf) - %d (%.2lf, %.2lf) # %.2lf %.2lf %.2lf\n",
	       t->v[0], pts[t->v[0]].x, pts[t->v[0]].y,
	       t->v[1], pts[t->v[1]].x, pts[t->v[1]].y,
	       t->v[2], pts[t->v[2]].x, pts[t->v[2]].y,
	       r2, r1, r3);

	printf("\tR = %lf, Rtolsq = %lf, C = %lf, C2 = %lf, beta = %lf, beta2 = %lf, Ctolsq = %lf, %s\n",
	       t->R, t->Rtolsq, t->C, t->C2, acos(t->C) * 180 / M_PI, acos(t->C2) * 180 / M_PI, t->Ctolsq,
	       t->clockwise ? "clockwise" : "counter-clockwise");



	return t;
}

static inline double euc_distance(double x1, double x2, double y1, double y2)
{
	return sqrt(sqr(x1 - x2) + sqr(y1 - y2));
}


static int fair_distances(struct point *pt, struct point pts[], int npts)
{
	int i;

	pt->dists = calloc(npts, sizeof(double));

	for (i = 0; i < npts; i++) {
		if ((fabs(pt->x - pts[i].x) < pos_thresh) ||
		    (fabs(pt->y - pts[i].y) < pos_thresh)) {
			return 0;
		}

		pt->dists[i] = euc_distance(pt->x, pts[i].x, pt->y, pts[i].y);
	}
	return 1;
}

/* Step 1.

   Pick up to N_OBJ while filtering points closer than POS_THRESHOLD.
   Each point will memoize distances to previous points in the resulting array.
*/
struct point *filter_point_list(GSList *l, int *npoints)
{
	struct point *pts = calloc(N_OBJ, sizeof(struct point));
	struct point *pt;
	int n;

	n = 0;
	while ((n < N_OBJ) && l) {
		pt = l->data;

		/* compute distances to previous points, tell us if we're too close */
		if (!fair_distances(pt, pts, n)) {
			l = g_slist_next(l);
			continue;
		}

		pts[n].x = pt->x;
		pts[n].y = pt->y;
		pts[n].dists = pt->dists;
		pt->dists = NULL;

		n++;
		l = g_slist_next(l);
	}

	*npoints = n;

	return pts;
}

/* Step 2.

   Generate list of triangles from a points array.
*/
GList *generate_triangles(struct point pts[], int npts)
{
	struct triangle *t;
	GList *tri = NULL;
	int i, j, k;

	for (i = 0; i < npts; i++) {
		for (j = i + 1; j < npts; j++) {
			for (k = j + 1; k < npts; k++) {
				t = triangle_new (pts, i, j, k);

				/* discard triangles with bad ratios */
				if (t->R > RATIO_LIMIT) {
					//printf("discarding ratio %.2lf\n", t->R);
					free(t);
				} else {
					tri = g_list_append(tri, t);
				}
			}
		}
	}

	return tri;
}



GSList *read_star_list(char *filename, double *maxx, double *maxy)
{
	static char buffer[1024];
	GSList *ret = NULL;
	FILE *fin;
	double x, y;
	struct point *s;

	if (!strcmp(filename, "-"))
		fin = stdin;
	else
		fin = fopen(filename, "r");

	if (fin == NULL)
		return NULL;

	while (fgets(buffer, 1024, fin)) {
		if (!isdigit(buffer[0]))
			continue;

		if (sscanf(buffer, "%lf %lf", &x, &y) == 2) {
			s = point_new (x, y);

			if (maxx && (x > *maxx))
				*maxx = x;

			if (maxy && (y > *maxy))
				*maxy = y;

			ret = g_slist_prepend (ret, s);
		}
	}
	fclose(fin);

	return g_slist_reverse (ret);
}


int main(int argc, char** argv)
{
	GSList *l1, *l2;
	GList *tri1, *tri2;
	int np1 = 0, np2 = 0;
	struct point *pt1, *pt2;

	l1 = read_star_list ((argc > 1) ? argv[1] : "starf1.list", &width, &height);
	l2 = read_star_list ((argc > 2) ? argv[2] : "starf2.list", &width, &height);

	pos_eps = POS_EPSILON * MAX(width, height);
	pos_thresh = POS_THRESHOLD * pos_eps;

	/* a) select points to be matched */
	pt1 = filter_point_list (l1, &np1);
	pt2 = filter_point_list (l2, &np2);

	printf("got %d points\n", np1);

	/* b) generate triangle lists */
	tri1 = generate_triangles (pt1, np1);
	tri2 = generate_triangles (pt2, np2);

	printf("got %d, %d triangles\n", g_list_length(tri1), g_list_length(tri2));

	return 0;
}

