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
#define COSINE_LIMIT 0.995
#define MAX_LOGM_ITER 40

#ifndef MAX
#define MAX(a, b) ((a) < (b) ? (b) : (a))
#endif

#define sqr(x) ((x) * (x))

struct point {
	int id;
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
	double Ctolsq;			  /* cosine tolerance, squared */
};

struct triangle_pair {
	struct triangle *t1;
	struct triangle *t2;
};

struct point_pair {

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

static inline struct triangle *triangle_new (struct point pts[], int v1, int v2, int v3, double epsilon)
{
	struct triangle *t = calloc(1, sizeof(struct triangle));
	int v[3], vv[3];
	double r1, r2, r3;

	vv[0] = v1;
	vv[1] = v2;
	vv[2] = v3;

	/* permute until ordered righ :-) */
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

	/* from the law of cosines */
	t->C = ( (pts[t->v[2]].x - pts[t->v[0]].x) * (pts[t->v[1]].x - pts[t->v[0]].x) +
		 (pts[t->v[2]].y - pts[t->v[0]].y) * (pts[t->v[1]].y - pts[t->v[0]].y) ) / (r2 * r3);

	double cterm = 1.0 / sqr(r3) - t->C / (r3 * r2) + 1.0 / sqr(r2);

	t->Rtolsq = 2 * sqr(t->R) * sqr(epsilon) * cterm;

	double sinsq = 1 - sqr(t->C);

	t->Ctolsq = 2 * sinsq * sqr(epsilon) * cterm + 3 * sqr(t->C) * sqr(sqr(epsilon)) * sqr(cterm);

	/* if det. is positive, we're clockwise. */
	double det = ((pts[t->v[1]].x - pts[t->v[0]].x) * (pts[t->v[2]].y - pts[t->v[0]].y) -
		      (pts[t->v[1]].y - pts[t->v[0]].y) * (pts[t->v[2]].x - pts[t->v[1]].x));

	t->clockwise = (det > 0);

	t->logP = log(r1 + r2 + r3);

	return t;
}

static inline double euc_distance(double x1, double x2, double y1, double y2)
{
	return sqrt(sqr(x1 - x2) + sqr(y1 - y2));
}


static int fair_distances(struct point *pt, struct point pts[], int npts, double epsilon)
{
	int i;

	pt->dists = calloc(npts, sizeof(double));

	for (i = 0; i < npts; i++) {
		if ((fabs(pt->x - pts[i].x) < POS_THRESHOLD * epsilon) ||
		    (fabs(pt->y - pts[i].y) < POS_THRESHOLD * epsilon)) {
			//printf("points too close, discarding: (%.2lf, %.2lf)   (%.2lf, %.2lf)\n", pt->x, pt->y, pts[i].x, pts[i].y);
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
struct point *filter_point_list(GSList *l, int *npoints, double epsilon)
{
	struct point *pts = calloc(N_OBJ, sizeof(struct point));
	struct point *pt;
	int n;
	int id;

	n = 0; id = -1;
	while ((n < N_OBJ) && l) {
		pt = l->data;
		id++;

		/* compute distances to previous points, tell us if we're too close */
		if (!fair_distances(pt, pts, n, epsilon)) {
			l = g_slist_next(l);
			continue;
		}

		pts[n].id = id;
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
GList *generate_triangles(struct point pts[], int npts, double epsilon)
{
	struct triangle *t;
	GList *tri = NULL;
	int i, j, k;

	for (i = 0; i < npts; i++) {
		for (j = i + 1; j < npts; j++) {
			for (k = j + 1; k < npts; k++) {
				t = triangle_new (pts, i, j, k, epsilon);

				/* discard triangles with bad ratios */
				if (t->R > RATIO_LIMIT || t->C > COSINE_LIMIT) {
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

/* Step 3.

   Match triangles.
*/
GList *match_triangles(GList *tri1, GList *tri2)
{
	GList *it1, *it2, *prev2 = tri2;
	struct triangle *t1, *t2, *tm;
	struct triangle_pair *pair;
	double rsq, csq, bm;
	GList *ret = NULL;

	int i1 = 0, i2 = 0, i3 = 0, i4 = 0;

	for (it1 = tri1; it1; it1 = g_list_next(it1)) {
		t1 = it1->data;

		//printf("t1 %d\n", i1);

		/* go to first usable triangle in second list */
		while (prev2) {
			t2 = prev2->data;
			rsq = sqr(t1->R - t2->R);
			if (rsq < t1->Rtolsq + t2->Rtolsq)
				break;

			prev2 = g_list_next(prev2);
			i2++;
		}

		it2 = prev2;
		tm = NULL;
		bm = 999999999999.0;

		i3 = i2;

		while (it2) {
			t2 = it2->data;

			rsq = sqr(t1->R - t2->R);
			csq = sqr(t1->C - t2->C);

			//printf("\teval %d - %d\n", i1, i3);

			if ((rsq < t1->Rtolsq + t2->Rtolsq) &&
			    (csq < t1->Ctolsq + t2->Ctolsq)) {

				//printf("match!\n");

				if (rsq + csq < bm) {
					//printf("\tmatch at %lf\n", rsq + csq);
					tm = t2;
					i4 = i3;
				}
			}

			if (rsq > t1->Rtolsq + t2->Rtolsq)
				break;

			it2 = g_list_next(it2);
			i3++;
		}

		if (tm) {
			pair = malloc(sizeof(struct triangle_pair));
			pair->t1 = t1;
			pair->t2 = tm;

			//printf("using match %d %d\n", i1, i4);
			ret = g_list_append(ret, pair);
		}

		i1++;
	}

	return ret;
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

static void print_triangle_list(struct point pts[], GList *triangles)
{
	struct triangle *t;

	for (; triangles; triangles = g_list_next(triangles)) {
		t = triangles->data;

		printf("%d (%.2lf, %.2lf) - %d (%.2lf, %.2lf) - %d (%.2lf, %.2lf)\n",
		       t->v[0], pts[t->v[0]].x, pts[t->v[0]].y,
		       t->v[1], pts[t->v[1]].x, pts[t->v[1]].y,
		       t->v[2], pts[t->v[2]].x, pts[t->v[2]].y);

		printf("\tR = %lf, Rtolsq = %lf, C = %lf, Ctolsq = %lf, logP = %lf, %s\n",
		       t->R, t->Rtolsq, t->C, t->Ctolsq, t->logP,
		       t->clockwise ? "clockwise" : "counter-clockwise");
	}
}


static void print_triangle_pairs(struct point pt1[], struct point pt2[], GList *tripairs)
{
	GList *tp;
	struct triangle_pair *t;
	FILE *fout = fopen("triangles.txt", "w+");

	for (tp = tripairs; tp; tp = g_list_next(tp)) {
		t = tp->data;

		fprintf(fout,
			"%.2lf %.2lf %.2lf %.2lf %.2lf %.2lf "
			"%.2lf %.2lf %.2lf %.2lf %.2lf %.2lf "
			"%lf %lf %lf %lf "
			"%lf %lf %lf %lf\n",

			pt1[t->t1->v[0]].x, pt1[t->t1->v[0]].y,
			pt1[t->t1->v[1]].x, pt1[t->t1->v[1]].y,
			pt1[t->t1->v[2]].x, pt1[t->t1->v[2]].y,

			pt2[t->t2->v[0]].x, pt2[t->t2->v[0]].y,
			pt2[t->t2->v[1]].x, pt2[t->t2->v[1]].y,
			pt2[t->t2->v[2]].x, pt2[t->t2->v[2]].y,

			t->t1->R, t->t1->Rtolsq, t->t1->C, t->t1->Ctolsq,
			t->t2->R, t->t2->Rtolsq, t->t2->C, t->t2->Ctolsq);


#if 0

		printf("Ra %lf, Rb %lf, dR %lf, tR %lf    Ca %lf, Cb %lf, dC %lf, tC %lf\n",
		       t->t1->R, t->t2->R,
		       sqr(t->t1->R - t->t2->R), t->t1->Rtolsq + t->t2->Rtolsq,
		       t->t1->C, t->t2->C,
		       sqr(t->t1->C - t->t2->C), t->t1->Ctolsq + t->t2->Ctolsq);
#endif
	}

	fclose(fout);

}

static int ratio_compare(struct triangle *t1, struct triangle *t2, double *maxtol)
{
	if (*maxtol < t1->Rtolsq)
		*maxtol = t1->Rtolsq;

	if (*maxtol < t2->Rtolsq)
		*maxtol = t2->Rtolsq;

	return (t1->R < t2->R) ? -1 : (t1->R > t2->R) ? 1 : 0;
}

int main(int argc, char** argv)
{
	GSList *l1, *l2;
	GList *tri1, *tri2;
	int np1 = 0, np2 = 0;
	struct point *pt1, *pt2;
	//struct triangle_pair *t;
	int i;
	double width, height, e1, e2;
	int discarded;
	GList *tp, *tmp;
	int np, nm;

	l1 = read_star_list ((argc > 1) ? argv[1] : "starf1.list", &width, &height);
	e1 = POS_EPSILON * MAX(width, height);

	l2 = read_star_list ((argc > 2) ? argv[2] : "starf2.list", &width, &height);
	e2 = POS_EPSILON * MAX(width, height);

	//printf("e1 %lf, e2 %lf\n", e1, e2);

	/* a) select points to be matched */
	pt1 = filter_point_list (l1, &np1, e1);
	pt2 = filter_point_list (l2, &np2, e2);

	//printf("got %d, %d points\n", np1, np2);

	/* b) generate triangle lists */
	tri1 = generate_triangles (pt1, np1, e1);
	tri2 = generate_triangles (pt2, np2, e2);

	printf("got %d, %d triangles\n", g_list_length(tri1), g_list_length(tri2));

	/* c) triangle match */

	/* sort in increasing order of ratio and determine max ratio tolerance */
	double rtol1 = 0.0, rtol2 = 0.0, maxtol;
	tri1 = g_list_sort_with_data (tri1, (GCompareDataFunc) ratio_compare, &rtol1);
	tri2 = g_list_sort_with_data (tri2, (GCompareDataFunc) ratio_compare, &rtol2);
	maxtol = rtol1 + rtol2;

	//print_triangle_list(pt1, tri1);
	//print_triangle_list(pt2, tri2);

	/* match the triangles */
	GList *tripairs = match_triangles(tri1, tri2);
	printf("got %d tripairs\n", g_list_length(tripairs));

	//print_triangle_pairs(pt1, pt2, tripairs);

	i = 0; discarded = 1;
	while (discarded && i < MAX_LOGM_ITER) {
		double avg, stdev, f = 0;

		np = 0; nm = 0;
		avg = 0.0;

		for (tp = tripairs; tp; tp = g_list_next(tp)) {
			struct triangle_pair *t = tp->data;


			//printf("delta log P %lf\n", t->t1->logP - t->t2->logP);
			avg += t->t1->logP - t->t2->logP;

			if (t->t1->clockwise == t->t2->clockwise)
				np++;
			else
				nm++;
		}

		avg /= (nm + np);
		stdev = 0.0;

		for (tp = tripairs; tp; tp = g_list_next(tp)) {
			struct triangle_pair *t = tp->data;
			stdev += sqr(t->t1->logP - t->t2->logP - avg);
		}

		stdev = sqrt(stdev / (nm + np));


		int mt = abs(np - nm);
		int mf = np + nm - mt;

		if (mf > mt)
			f = 1.0;
		else if (0.1 * mt > 1.0 * mf)
			f = 3.0;
		else
			f = 2.0;

		printf("iter %3d, avg %lf, stdev %lf, n+ %5d, n- %5d, factor %lf\n", i, avg, stdev, np, nm, f);

		discarded = 0;

		for (tp = tripairs; tp; ) {
			struct triangle_pair *t = tp->data;

			if (fabs(t->t1->logP - t->t2->logP - avg) > f * stdev) {

				discarded = 1;

				tmp = tp->next;
				tripairs = g_list_delete_link (tripairs, tp);
				tp = tmp;
				continue;
			}

			tp = g_list_next(tp);
		}

		i++;
	}

	if (discarded) {
		printf("failed\n");
	}

	for (tp = tripairs; tp; ) {
		struct triangle_pair *t = tp->data;

		if ((np > nm && t->t1->clockwise != t->t2->clockwise) ||
		    (nm > np && t->t1->clockwise == t->t2->clockwise)) {

			tmp = tp->next;
			tripairs = g_list_delete_link (tripairs, tp);
			tp = tmp;
			continue;
		}

		tp = g_list_next(tp);
	}

#if 0
	for (tp = tripairs; tp; tp = g_list_next(tp)) {
		struct triangle_pair *t = tp->data;

		printf("%d %d\n", pt1[t->t1->v[0]].id, pt2[t->t2->v[0]].id);
		printf("%d %d\n", pt1[t->t1->v[1]].id, pt2[t->t2->v[1]].id);
		printf("%d %d\n", pt1[t->t1->v[2]].id, pt2[t->t2->v[2]].id);
	}
#endif

	print_triangle_pairs(pt1, pt2, tripairs);



	printf("rtol1 = %lf, rtol2 = %lf, %d matches\n", rtol1, rtol2, g_list_length(tripairs));

	return 0;
}

