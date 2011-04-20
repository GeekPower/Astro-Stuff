#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <math.h>

#include <glib.h>
#include <gtk/gtk.h>

static GtkWidget *mwin, *darea;

static GSList *l1 = NULL, *l2 = NULL, *l3 = NULL;
static double max_x = 0, max_y = 0;

struct star {
	double x, y;

	struct star *pair;
};

void on_window_destroy (GtkObject *obj, gpointer data)
{
	printf("boom\n");
	gtk_main_quit ();
}

#define STARR 1.5
static void draw_star_list(cairo_t *cr, GSList *stars,
			   double xs, double ys, double radius, double r, double g, double b)
{
	GSList *l;
	struct star *s;
	double dashes[] = { 5.0, 5.0 };

	cairo_set_source_rgb (cr, r, g, b);

	for (l = stars; l; l = g_slist_next(l)) {
		s = (struct star *) l->data;

		cairo_arc (cr, s->x * xs, s->y * ys, radius, 0, 2 * M_PI);
		cairo_stroke (cr);

		/* draw pairing */
		if (s->pair) {
			cairo_save (cr);

			cairo_set_source_rgb (cr, 0.0, 0.0, 1.0);

			cairo_set_dash (cr, dashes, 2, 0);
			cairo_move_to (cr, s->x * xs, s->y * ys);
			cairo_line_to (cr, s->pair->x * xs, s->pair->y * ys);
			cairo_stroke (cr);

			cairo_restore (cr);
		}
	}
}

gboolean on_darea_expose (GtkObject *obj, GdkEventExpose *ev, gpointer data)
{
	cairo_t *cr;
	double xs = darea->allocation.width / max_x;
	double ys = darea->allocation.height / max_y;

	cr = gdk_cairo_create (darea->window);

	/* clear the area */
        cairo_set_source_rgb (cr, 1.0, 1.0, 1.0);
        cairo_rectangle (cr, 0, 0, darea->allocation.width, darea->allocation.height);
        cairo_fill (cr);

	cairo_set_line_width (cr, 1.0);

	/* first list red */
	draw_star_list (cr, l1, xs, ys, 1.0, 1.0, 0.0, 0.0);

	/* second list green */
	draw_star_list (cr, l2, xs, ys, 2.0, 0.0, 1.0, 0.0);

	cairo_destroy (cr);

	return TRUE;
}

struct star *star_new(double x, double y)
{
	struct star *s = g_malloc(sizeof(struct star));

	memset (s, 0, sizeof(struct star));

	s->x = x;
	s->y = y;

	return s;
}

GSList *read_list(char *filename, double *maxx, double *maxy)
{
	static char buffer[1024];
	GSList *ret = NULL;
	FILE *fin;
	double x, y;
	struct star *s;

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
			s = star_new (x, y);

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

/* pairs are not really stars, but indices in l1/l2 */
static void do_pairing(GSList *l1, GSList *l2, GSList *pairs)
{
	struct star *s1, *s2;

	for (;pairs;pairs = g_slist_next(pairs)) {
		struct star *p = pairs->data;

		int id1 = (int) p->x;
		int id2 = (int) p->y;

		/* yuck */
		if ((s1 = g_slist_nth_data(l1, id1)) &&
		    (s2 = g_slist_nth_data(l2, id2))) {
			s1->pair = s2;
		}
	}
}

int main(int argc, char **argv)
{
	GtkBuilder *builder;
	GtkWidget *closebutton;

	if (argc < 4)
		exit(1);

	l1 = read_list(argv[1], &max_x, &max_y);
	l2 = read_list(argv[2], &max_x, &max_y);

	/* cheat */
	l3 = read_list(argv[3], NULL, NULL);
	do_pairing(l1, l2, l3);

	gtk_init(&argc, &argv);

	builder = gtk_builder_new ();
	gtk_builder_add_from_file (builder, "matchviz.xml", NULL);

	mwin = GTK_WIDGET (gtk_builder_get_object (builder, "window1"));
	darea = GTK_WIDGET (gtk_builder_get_object (builder, "darea"));
	closebutton = GTK_WIDGET (gtk_builder_get_object (builder, "closebutton"));

	//gtk_builder_connect_signals (builder, NULL);

	g_signal_connect (G_OBJECT (mwin), "destroy",
			  G_CALLBACK (on_window_destroy), NULL);
	g_signal_connect (G_OBJECT (darea), "expose_event",
			  G_CALLBACK (on_darea_expose), NULL);
	g_signal_connect (G_OBJECT (closebutton), "clicked",
			  G_CALLBACK (on_window_destroy), NULL);

	g_object_unref (builder);

	gtk_widget_show (mwin);
	gtk_main ();

	return 0;
}
