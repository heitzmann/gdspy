/**********************************************************************
 *                                                                    *
 *  Copyright 2009-2016 Lucas Heitzmann Gabrielli                     *
 *                                                                    *
 *  This file is part of gdspy.                                       *
 *                                                                    *
 *  gdspy is free software: you can redistribute it and/or modify it  *
 *  under the terms of the GNU General Public License as published    *
 *  by the Free Software Foundation, either version 3 of the          *
 *  License, or any later version.                                    *
 *                                                                    *
 *  gdspy is distributed in the hope that it will be useful, but      *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     *
 *  GNU General Public License for more details.                      *
 *                                                                    *
 *  You should have received a copy of the GNU General Public         *
 *  License along with gdspy.  If not, see                            *
 *  <http://www.gnu.org/licenses/>.                                   *
 *                                                                    *
 **********************************************************************/

/* Compiled with:
gcc -O3 -fPIC -shared -I/usr/include/python2.7 -L/usr/lib -lpython2.7 -o boolext.so boolext.c
*/

#include <Python.h>

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

double eps = 1e-13;

long num_polygons;

struct segment_struct;

typedef struct point_struct {
	struct point_struct *nxt1;
	struct point_struct *nxt2;
	struct point_struct *nxt3;
	struct segment_struct **segment;
	double x[2];
	long max_seg;
	char flag;
} point_type;

typedef struct segment_struct {
	struct segment_struct *nxt1;
	struct segment_struct *nxt2;
	point_type *point[2];
	short *polygon;
	short output;
} segment_type;

typedef struct polygon_struct {
	struct polygon_struct *nxt;
	struct polygon_struct *inside;
	struct point_struct *point;
	double bbox[4];
} polygon_type;

#ifdef DEBUG
#define FLAG printf("DEBUG %ld\n", ++dbg);fflush(stdout);
#define POINT(m,a) printf("%s%lx: %g,%g : ",m,a,a->x[0],a->x[1]);for(itmp=0;itmp<a->max_seg&&a->segment[itmp]!=NULL;++itmp)printf("%lx,",a->segment[itmp]);printf("\b\n");fflush(stdout);
#define POINTS(m,a,b) printf(m);for(ptmp=(a);ptmp!=(b);ptmp=ptmp->nxt1){printf("%lx: %g,%g : ",ptmp,ptmp->x[0],ptmp->x[1]);for(itmp=0;itmp<ptmp->max_seg&&ptmp->segment[itmp]!=NULL;++itmp)printf("%lx,",ptmp->segment[itmp]);printf("\b\n");}fflush(stdout);
#define POINTS_ALT(m,a,b) printf(m);for(ptmp=(a);ptmp!=(b);ptmp=ptmp->nxt2){printf("%lx: %g,%g : ",ptmp,ptmp->x[0],ptmp->x[1]);for(itmp=0;itmp<ptmp->max_seg&&ptmp->segment[itmp]!=NULL;++itmp)printf("%lx,",ptmp->segment[itmp]);printf("\b\n");}fflush(stdout);
#define SEGMENT(m,a) printf("%s%lx: %g,%g ~ %g,%g (%lx ~ %lx) [",m,a,a->point[0]->x[0],a->point[0]->x[1],a->point[1]->x[0],a->point[1]->x[1],a->point[0],a->point[1]);for(itmp=0;itmp<num_polygons;++itmp)printf("%hd,",a->polygon[itmp]);printf("\b] %hd\n",a->output);fflush(stdout);
#define SEGMENTS(m,a,b) printf(m);for(stmp=(a);stmp!=(b);stmp=stmp->nxt1){printf("%lx: %g,%g	~  %g,%g (%lx ~ %lx) [",stmp,stmp->point[0]->x[0],stmp->point[0]->x[1],stmp->point[1]->x[0],stmp->point[1]->x[1],stmp->point[0],stmp->point[1]);for(itmp=0;itmp<num_polygons;++itmp)printf("%hd,",stmp->polygon[itmp]);printf("\b] %hd\n", stmp->output);}fflush(stdout);
#define SEGMENTS_ALT(m,a,b) printf(m);for(stmp=(a);stmp!=(b);stmp=stmp->nxt2){printf("%lx: %g,%g	~  %g,%g (%lx ~ %lx) [",stmp,stmp->point[0]->x[0],stmp->point[0]->x[1],stmp->point[1]->x[0],stmp->point[1]->x[1],stmp->point[0],stmp->point[1]);for(itmp=0;itmp<num_polygons;++itmp)printf("%hd,",stmp->polygon[itmp]);printf("\b] %hd\n", stmp->output);}fflush(stdout);
#define PAUSE gets(str);
point_type *ptmp;
segment_type *stmp;
char str[256];
long itmp;
long dbg = 0;
void POLYGONS(int i, polygon_type *a) {
	polygon_type *g1;
	int j;
	for(g1=a;g1!=NULL;g1=g1->nxt)
	{  
		j = i;
		while (0 < j--) printf("	");
		printf("[%lx] : ",g1);
		for(ptmp=g1->point;ptmp!=NULL;ptmp=ptmp->nxt2) printf("(%lx: %g,%g) ",ptmp,ptmp->x[0],ptmp->x[1]);
		printf("\b\n");
		POLYGONS(i+1,g1->inside);
	}
}
#endif

point_type * new_point(double x, double y)
{
	point_type *point;
	point = (point_type *) malloc(sizeof(point_type));
	point->nxt1 = NULL;
	point->nxt3 = NULL;
	point->x[0] = x;
	point->x[1] = y;
	point->flag = 0;
	point->max_seg = 2;
	point->segment = (segment_type **) calloc(2, sizeof(segment_type *));
	return point;
}


point_type * duplicate_point(point_type *p1)
{
	point_type *point;
	point = (point_type *) malloc(sizeof(point_type));
	point->nxt1 = p1->nxt1;
	point->nxt2 = p1->nxt2;
	point->nxt3 = p1;
	p1->nxt1 = point;
	point->x[0] = p1->x[0];
	point->x[1] = p1->x[1];
	point->flag = 0;
	point->max_seg = 0;
	point->segment = NULL;
	return point;
}


segment_type * new_segment(point_type *p1, point_type *p2, short *polygon)
{
	segment_type *seg;
	seg = (segment_type *) malloc(sizeof(segment_type));
	seg->nxt1 = NULL;
	seg->point[0] = p1;
	seg->point[1] = p2;
	if (polygon != NULL)
	{
		seg->polygon = (short *) malloc(num_polygons * sizeof(short));
		memcpy(seg->polygon, polygon, num_polygons * sizeof(short));
	}
	else seg->polygon = (short *) calloc(num_polygons, sizeof(short));
	return seg;
}


polygon_type * new_polygon(point_type *point)
{
	polygon_type *polygon;
	double area = 0;
	point_type *p1, *p2, *p3;
	polygon = (polygon_type *) malloc(sizeof(polygon_type));
	polygon->nxt = polygon->inside = NULL;
	polygon->bbox[0] = polygon->bbox[1] = point->x[0];
	polygon->bbox[2] = polygon->bbox[3] = point->x[1];
	for (p1 = point->nxt2, p2 = point->nxt2->nxt2; p2 != NULL; p1 = p1->nxt2, p2 = p2->nxt2)
		area += (p1->x[0] - point->x[0]) * (p2->x[1] - point->x[1]) - (p1->x[1] - point->x[1]) * (p2->x[0] - point->x[0]);
	if (area < 0)
	{
		p1 = point;
		p3 = NULL;
		while (p1 != NULL)
		{
			if (p1->x[0] < polygon->bbox[0]) polygon->bbox[0] = p1->x[0];
			else if (p1->x[0] > polygon->bbox[1]) polygon->bbox[1] = p1->x[0];
			if (p1->x[1] < polygon->bbox[2]) polygon->bbox[2] = p1->x[1];
			else if (p1->x[1] > polygon->bbox[3]) polygon->bbox[3] = p1->x[1];
			p2 = p3;
			p3 = p1;
			p1 = p1->nxt2;
			p3->nxt2 = p2;
		}
		polygon->point = p3;
	}
	else
	{
		for (p1 = point; p1 != NULL; p1 = p1->nxt2)
		{
			if (p1->x[0] < polygon->bbox[0]) polygon->bbox[0] = p1->x[0];
			else if (p1->x[0] > polygon->bbox[1]) polygon->bbox[1] = p1->x[0];
			if (p1->x[1] < polygon->bbox[2]) polygon->bbox[2] = p1->x[1];
			else if (p1->x[1] > polygon->bbox[3]) polygon->bbox[3] = p1->x[1];
		}
		polygon->point = point;
	}
	return polygon;
}


short compare_points(point_type *p1, point_type  *p2)
{
	if (p1->x[0] < p2->x[0] - eps) return -1;
	else if (p1->x[0] > p2->x[0] + eps) return 1;
	else if (p1->x[1] < p2->x[1] - eps) return -1;
	else if (p1->x[1] > p2->x[1] + eps) return 1;
	else return 0;
}


short compare_segments(segment_type *s1, segment_type *s2)
{
	double x, y1, y2, dx;
	x = 0.5 * MAX(s1->point[0]->x[0], s2->point[0]->x[0]) + 0.5 * MIN(s1->point[1]->x[0], s2->point[1]->x[0]);
	dx = s1->point[1]->x[0] - s1->point[0]->x[0];
	if (dx < eps && dx > -eps) y1 = 0.5 * s1->point[0]->x[1] + 0.5 * s1->point[1]->x[1];
	else y1 = s1->point[0]->x[1] + (s1->point[1]->x[1] - s1->point[0]->x[1]) * (x - s1->point[0]->x[0]) / dx;
	dx = s2->point[1]->x[0] - s2->point[0]->x[0];
	if (dx < eps && dx > -eps) y2 = 0.5 * s2->point[0]->x[1] + 0.5 * s2->point[1]->x[1];
	else y2 = s2->point[0]->x[1] + (s2->point[1]->x[1] - s2->point[0]->x[1]) * (x - s2->point[0]->x[0]) / dx;
	if (y1 < y2) return -1;
	else if (y1 > y2) return 1;
	return 0;
}


short compare_adjacent_segments(segment_type *s1, segment_type	*s2)
{
	if (s1->point[0] == s2->point[0]) return compare_points(s1->point[1], s2->point[1]);
	else if (s1->point[0] == s2->point[1]) return compare_points(s1->point[1], s2->point[0]);
	else if (s1->point[1] == s2->point[0]) return compare_points(s1->point[0], s2->point[1]);
	else return compare_points(s1->point[0], s2->point[0]);
}


long add_segment(point_type *p, segment_type *s)
{
	long i;
	i = p->max_seg;
	if (p->segment[i - 1] != NULL)
	{
		p->max_seg *= 2;
		p->segment = (segment_type **) realloc(p->segment, sizeof(segment_type *) * p->max_seg);
		memset(&(p->segment[i]), 0, i * sizeof(segment_type *));
	}
	while (i > 0 && p->segment[i - 1] == NULL) --i;
	p->segment[i] = s;
	return i + 1;
}


long remove_segment(point_type *p, segment_type *s)
{
	long i = 0;
	while (p->segment[i] != s) ++i;
	++i;
	while (i < p->max_seg && p->segment[i] != NULL) 
	{
		p->segment[i - 1] = p->segment[i];
		++i;
	}
	p->segment[--i] = NULL;
	return i;
}


char inside(double x, double y, polygon_type *poly)
{
	char result = 0;
	char stop = 0;
	point_type *p1, *p2;
	p1 = poly->point;
	p2 = p1->nxt2;
	while (stop == 0)
	{
		if (p2 == NULL) 
		{
			p2 = poly->point;
			stop = 1;
		}
		if (((p2->x[1] < y && p1->x[1] >= y) || (p1->x[1] < y && p2->x[1] >= y)) && 
				(p1->x[0] - p2->x[0]) * (y - p2->x[1]) / (p1->x[1] - p2->x[1]) - x + p2->x[0] < 0)
			result = 1 - result;
		p1 = p2;
		p2 = p2->nxt2;
	}
	return result;
}


void add_inside(polygon_type *in, polygon_type *out)
{
	polygon_type *g0, *g1;
	g0 = NULL;
	g1 = out->inside;
	while (g1 != NULL)
	{
		if (in->bbox[0] == g1->bbox[0] && in->bbox[1] == g1->bbox[1] && in->bbox[2] == g1->bbox[2] && in->bbox[3] == g1->bbox[3])
		{
			if (inside(0.5 * (in->point->x[0] + in->point->nxt2->x[0]), 0.5 * (in->point->x[1] + in->point->nxt2->x[1]), g1) > 0)
			{
				add_inside(in, g1);
				return;
			}
			if (inside(0.5 * (g1->point->x[0] + g1->point->nxt2->x[0]), 0.5 * (g1->point->x[1] + g1->point->nxt2->x[1]), in) > 0)
			{
				if (g0 == NULL) out->inside = g1->nxt;
				else g0->nxt = g1->nxt;
				g1->nxt = in->inside;
				in->inside = g1;
				g1 = (g0 == NULL ? out->inside : g0->nxt);
			}
		}
		else if (in->bbox[0] >= g1->bbox[0] && in->bbox[1] <= g1->bbox[1] && in->bbox[2] >= g1->bbox[2] && in->bbox[3] <= g1->bbox[3] &&
				inside(0.5 * (in->point->x[0] + in->point->nxt2->x[0]), 0.5 * (in->point->x[1] + in->point->nxt2->x[1]), g1) > 0)
		{
			add_inside(in, g1);
			return;
		}
		else if (in->bbox[0] <= g1->bbox[0] && in->bbox[1] >= g1->bbox[1] && in->bbox[2] <= g1->bbox[2] && in->bbox[3] >= g1->bbox[3] &&
				inside(0.5 * (g1->point->x[0] + g1->point->nxt2->x[0]), 0.5 * (g1->point->x[1] + g1->point->nxt2->x[1]), in) > 0)
		{
			if (g0 == NULL) out->inside = g1->nxt;
			else g0->nxt = g1->nxt;
			g1->nxt = in->inside;
			in->inside = g1;
			g1 = (g0 == NULL ? out->inside : g0->nxt);
		}
		else
		{
			g0 = g1;
			g1 = g1->nxt;
		}
	}
	in->nxt = out->inside;
	out->inside = in;
}


point_type * free_point(point_type *point)
{
	point_type *nxt1;
	nxt1 = point->nxt1;
	free(point->segment);
	free(point);
	return nxt1;
}


segment_type * free_segment(segment_type *segment)
{
	segment_type *nxt1;
	nxt1 = segment->nxt1;
	free(segment->polygon);
	free(segment);
	return nxt1;
}


void free_point_cycle(point_type *point)
{
	point_type *p;
	if (point != NULL)
	{
		p = point->nxt1;
		point->nxt1 = NULL;
		point = p;
		while (point != NULL)
		{
			p = point->nxt1;
			free(point->segment);
			free(point);
			point = p;
		}
	}
}


void free_point_list(point_type *point)
{
	point_type *p;
	while (point != NULL)
	{
		p = point->nxt1;
		free(point->segment);
		free(point);
		point = p;
	}
}


void free_segment_list(segment_type *segment)
{
	segment_type *s;
	while (segment != NULL)
	{
		s = segment->nxt1;
		free(segment->polygon);
		free(segment);
		segment = s;
	}
}


void free_polygon_list(polygon_type *poly)
{
	polygon_type *g0;
	while (poly != NULL)
	{
		g0 = poly->nxt;
		free(poly);
		poly = g0;
	}
}


point_type * mergesort_points(point_type *p1)
{
	point_type head, *tail, *p2;
	long step = 1, n1, n2;
	long merges = 2;

	head.nxt1 = p1;
	while (merges > 1)
	{
		p1 = head.nxt1;
		tail = &head;
		merges = 0;
		while (p1 != NULL)
		{
			++merges;
			p2 = p1;
			n1 = 0;
			while (p2 != NULL && n1 < step) 
			{
				++n1;
				p2 = p2->nxt1;
			}
			n2 = step;
			while (n1 > 0 || (n2 > 0 && p2 != NULL))
			{
				if (n1 == 0)
				{
					--n2;
					tail->nxt1 = p2;
					p2 = p2->nxt1;
				}
				else if (n2 == 0 || p2 == NULL)
				{
					--n1;
					tail->nxt1 = p1;
					p1 = p1->nxt1;
				}
				else if (compare_points(p2, p1) < 0)
				{
					--n2;
					tail->nxt1 = p2;
					p2 = p2->nxt1;
				}
				else 
				{
					--n1;
					tail->nxt1 = p1;
					p1 = p1->nxt1;
				}
				tail = tail->nxt1;
			}
			p1 = p2;
		}
		tail->nxt1 = NULL;
		step *= 2;
	}
	return head.nxt1;
}


point_type * merge_sorted_points(point_type *p1, point_type *p2)
{
	point_type head, *tail;
	tail = &head;
	while (p1 != NULL || p2 != NULL)
	{
		if (p1 == NULL)
		{
			tail->nxt1 = p2;
			p2 = p2->nxt1;
		}
		else if (p2 == NULL)
		{
			tail->nxt1 = p1;
			p1 = p1->nxt1;
		}
		else if (compare_points(p2, p1) < 0)
		{
			tail->nxt1 = p2;
			p2 = p2->nxt1;
		}
		else
		{
			tail->nxt1 = p1;
			p1 = p1->nxt1;
		}
		tail = tail->nxt1;
	}
	tail->nxt1 = NULL;
	return head.nxt1;
}


segment_type * sort_segments(segment_type *seg)
{
	segment_type head;
	segment_type *s0;
	head.nxt2 = NULL;
	while (seg != NULL)
	{
		head.nxt1 = seg->nxt2;
		s0 = &head;
		while (s0->nxt2 != NULL && compare_segments(seg, s0->nxt2) > 0) s0 = s0->nxt2;
		seg->nxt2 = s0->nxt2;
		s0->nxt2 = seg;
		seg = head.nxt1;
	}
	return head.nxt2;
}


polygon_type * mergesort_polygons(polygon_type *g1)
{
	polygon_type head, *tail, *g2;
	long step = 1, n1, n2;
	long merges = 2;

	head.nxt = g1;
	while (merges > 1)
	{
		g1 = head.nxt;
		tail = &head;
		merges = 0;
		while (g1 != NULL)
		{
			++merges;
			g2 = g1;
			n1 = 0;
			while (g2 != NULL && n1 < step) 
			{
				++n1;
				g2 = g2->nxt;
			}
			n2 = step;
			while (n1 > 0 || (n2 > 0 && g2 != NULL))
			{
				if (n1 == 0)
				{
					--n2;
					tail->nxt = g2;
					g2 = g2->nxt;
				}
				else if (n2 == 0 || g2 == NULL)
				{
					--n1;
					tail->nxt = g1;
					g1 = g1->nxt;
				}
				else if (g2->bbox[0] < g1->bbox[0])
				{
					--n2;
					tail->nxt = g2;
					g2 = g2->nxt;
				}
				else 
				{
					--n1;
					tail->nxt = g1;
					g1 = g1->nxt;
				}
				tail = tail->nxt;
			}
			g1 = g2;
		}
		tail->nxt = NULL;
		step *= 2;
	}
	return head.nxt;
}


long intersect_point_segment(point_type *p1)
{
	long intersections = 0;
	long i;
	double dx0, dy0, dx1, dy1;
	double seg_len2, cross2, dot;
	point_type *p2;
	segment_type *s2, *s1;
	segment_type *prev, *tail;
	segment_type beam;

	tail = &beam;
	for (i = 0; i < p1->max_seg && p1->segment[i] != NULL; ++i) tail = tail->nxt2 = p1->segment[i];
	tail->nxt2 = NULL;

	for (p2 = p1->nxt1; p2 != NULL; p1 = p1->nxt1, p2 = p2->nxt1)
	{
		prev = &beam;
		s1 = beam.nxt2;
		while (s1 != NULL)
		{
			dx0 = s1->point[1]->x[0] - s1->point[0]->x[0];
			dy0 = s1->point[1]->x[1] - s1->point[0]->x[1];
			dx1 = p2->x[0] - s1->point[0]->x[0];
			dy1 = p2->x[1] - s1->point[0]->x[1];
			seg_len2 = dx0 * dx0 + dy0 * dy0;
			cross2 = (dx1 * dy0 - dx0 * dy1);
			cross2 *= cross2;
			dot = dx1 * dx0 + dy1 * dy0;

			//if (cross2 / eps < eps * seg_len2 && dot > 0 && dot < seg_len2 && p2 != s1->point[0] && p2 != s1->point[1])
			if (cross2 < eps * eps * seg_len2 && dot > 0 && dot < seg_len2 && p2 != s1->point[0] && p2 != s1->point[1])
			{
				++intersections;
				s2 = new_segment(p2, s1->point[1], s1->polygon);
				s2->nxt1 = s1->nxt1;
				s1->nxt1 = s2;
				remove_segment(s1->point[1], s1);
				add_segment(s1->point[1], s2);
				add_segment(p2, s2);
				add_segment(p2, s1);
				s1->point[1] = p2;
			}
			else 
			{
				prev = prev->nxt2;
				s1 = s1->nxt2;
			}
		}
		/* Update scanbeam */
		for (i = 0; i < p2->max_seg && p2->segment[i] != NULL; ++i)
		{
			s2 = p2->segment[i];
			if (s2->point[0] == p2)
			{
				tail = tail->nxt2 = s2;
				tail->nxt2 = NULL;
			}
			else if (s2->point[1] == p2)
			{
				prev = &beam;
				while (prev->nxt2 != s2) prev = prev->nxt2;
				prev->nxt2 = s2->nxt2;
				if (s2 == tail) tail = prev;
			}
		}
	}
	return intersections;
}


static PyObject * clip(PyObject *self, PyObject *args)
{
	PyObject *in_polygons, *operation;
	PyObject *obj1, *obj2, *obj3, *obj4, *obj5;
	long num_points;
	long i, j;
	long result, last_result;
	long *counter;
	short left;
	short *phase;
	double den;
	double x1, dx, dx0, dx1;
	double y1, dy, dy0, dy1;
	double alpha = 0, beta = 0;
	point_type head;
	point_type *p0, *p1, *p2, *p3, *p4;
	point_type *point_first, *point_last;
	segment_type beam, removed;
	segment_type *tail;
	segment_type *s0, *s1;
	segment_type *segment_first, *segment_last;
	polygon_type polygon_tree;
	polygon_type *g0, *g1, *g2;
	polygon_type *polygon_last, *hole_prev1, *hole_prev2;

	if (!PyArg_ParseTuple(args, "OOd:clip", &in_polygons, &operation, &eps)) return NULL;
	if (!PySequence_Check(in_polygons))
	{
		PyErr_SetString(PyExc_TypeError, "First argument must be a sequence.");
		return NULL;
	}
	if (!PyCallable_Check(operation))
	{
		PyErr_SetString(PyExc_TypeError, "Second argument must be callable.");
		return NULL;
	}

	/* Add input polygons to graph structure: points/segments */
	p0 = p1 = NULL;
	segment_last = &beam;
	segment_last->nxt1 = NULL;
	point_last = &head;
	point_last->nxt1 = NULL;
	num_polygons = PySequence_Length(in_polygons);
	for (i = num_polygons - 1; i >= 0; --i)
	{
		if ((obj1 = PySequence_ITEM(in_polygons, i)) == NULL)
		{
			free_point_cycle(p0);
			free_point_list(head.nxt1);
			free_segment_list(beam.nxt1);
			return NULL;
		}
		if (!PySequence_Check(obj1))
		{
			free_point_cycle(p0);
			free_point_list(head.nxt1);
			free_segment_list(beam.nxt1);
			Py_DECREF(obj1);
			PyErr_SetString(PyExc_TypeError, "Elements of the first argument must be sequences.");
			return NULL;
		}

		/* Parse a complete polygon */
		num_points = PySequence_Length(obj1);
		for (j = num_points - 1; j >= 0; --j)
		{
			if ((obj2 = PySequence_ITEM(obj1, j)) == NULL)
			{
				free_point_cycle(p0);
				free_point_list(head.nxt1);
				free_segment_list(beam.nxt1);
				Py_DECREF(obj1);
				return NULL;
			}

			if ((obj3 = PySequence_GetItem(obj2, 0)) == NULL)
			{
				free_point_cycle(p0);
				free_point_list(head.nxt1);
				free_segment_list(beam.nxt1);
				Py_DECREF(obj2);
				Py_DECREF(obj1);
				return NULL;
			}
			x1 = PyFloat_AsDouble(obj3);
			Py_DECREF(obj3);

			if ((obj3 = PySequence_GetItem(obj2, 1)) == NULL)
			{
				free_point_cycle(p0);
				free_point_list(head.nxt1);
				free_segment_list(beam.nxt1);
				Py_DECREF(obj2);
				Py_DECREF(obj1);
				return NULL;
			}
			y1 = PyFloat_AsDouble(obj3);
			Py_DECREF(obj3);
			Py_DECREF(obj2);

			if (p0 == NULL)
			{
				p0 = new_point(x1, y1);
				p1 = p0;
			}
			else 
			{
				p1->nxt1 = new_point(x1, y1);
				p1 = p1->nxt1;
			}
		}
		Py_DECREF(obj1);
		if (num_points > 0) p1->nxt1 = p0;

		/* Remove 0-area triangles formed by any 3 consecutive vertices */
		result = num_points - 2;
		while (result > 0)
		{
			result = 0;
			p1 = p0;
			do
			{
				p2 = p1->nxt1;
				alpha = (p1->x[0] - p2->x[0]) * (p2->nxt1->x[1] - p2->x[1]) - (p2->nxt1->x[0] - p2->x[0]) * (p1->x[1] - p2->x[1]);
				if (alpha < 0.5 * eps && alpha > -0.5 * eps)
				{
					p1->nxt1 = p2->nxt1;
					if (p2 == p0) p0 = p0->nxt1;
					free(p2);
					result = --num_points - 2;
				}
				p1 = p1->nxt1;
			}
			while (p1 != p0 && num_points > 2);
		}

		/* Append polygon to graph structure */
		if (num_points > 2)
		{
			for (p1 = p0, j = num_points; j > 0; p1 = p1->nxt1, --j)
			{
				segment_last->nxt1 = new_segment(p1->nxt1, p1, NULL);
				segment_last = segment_last->nxt1;
				segment_last->polygon[i] = compare_points(p1, p1->nxt1);
				if (segment_last->polygon[i] < 0)
				{
					segment_last->point[0] = p1;
					segment_last->point[1] = p1->nxt1;
				}
				add_segment(p1, segment_last);
				add_segment(p1->nxt1, segment_last);
			}

			point_last->nxt1 = p0->nxt1;
			p0->nxt1 = NULL;
			point_last = p0;
		}
		else free_point_cycle(p0);
		p0 = NULL;
	}

	if (head.nxt1 == NULL)
	{
		Py_INCREF(Py_None);
		return Py_None;
	}

	/* Sort points by x-coordinate (y, if same x) */
	point_first = mergesort_points(head.nxt1);
	segment_first = beam.nxt1;

	/* Merge overlaping points */
	p1 = point_first;
	p2 = point_first->nxt1;
	while (p2 != NULL)
	{
		if (compare_points(p1, p2) == 0)
		{
			for (i = 0; i < p2->max_seg && p2->segment[i] != NULL; ++i)
			{
				s1 = p2->segment[i];
				add_segment(p1, s1);
				if (s1->point[0] == p2) s1->point[0] = p1;
				else s1->point[1] = p1;
			}
			p1->nxt1 = p2 = free_point(p2);
		}
		else
		{
			p1 = p2;
			p2 = p2->nxt1;
		}
	}

	/* Calculate intersections between segments and points */
	result = 1;
	while (result > 0)
	{
		result = intersect_point_segment(point_first);
		for (p1 = point_first; p1 != NULL; p1 = p1->nxt1)
		{
			x1 = p1->x[0];
			p1->x[0] = p1->x[1];
			p1->x[1] = x1;
		}
		point_first = mergesort_points(point_first);
		for (s1 = segment_first; s1 != NULL; s1 = s1->nxt1) if (compare_points(s1->point[0], s1->point[1]) > 0)
		{
			p1 = s1->point[1];
			s1->point[1] = s1->point[0];
			s1->point[0] = p1;
		}

		result += intersect_point_segment(point_first);
		for (p1 = point_first; p1 != NULL; p1 = p1->nxt1)
		{
			x1 = p1->x[0];
			p1->x[0] = p1->x[1];
			p1->x[1] = x1;
		}
		point_first = mergesort_points(point_first);
		for (s1 = segment_first; s1 != NULL; s1 = s1->nxt1) if (compare_points(s1->point[0], s1->point[1]) > 0)
		{
			p1 = s1->point[1];
			s1->point[1] = s1->point[0];
			s1->point[0] = p1;
		}
	}

	/* Add points to segment crossings */
	p2 = NULL;
	tail = &beam;
	tail->nxt2 = NULL;
	for (p1 = point_first; p1 != NULL; p1 = p1->nxt1)
	{
		/* Update scanbeam */
		for (i = 0; i < p1->max_seg && p1->segment[i] != NULL; ++i)
		{
			s1 = p1->segment[i];
			if (s1->point[0] == p1)
			{
				tail->nxt2 = s1;
				tail = tail->nxt2;
				tail->nxt2 = NULL;
			}
			else if (s1->point[1] == p1)
			{
				s0 = &beam;
				while (s0->nxt2 != s1) s0 = s0->nxt2;
				s0->nxt2 = s1->nxt2;
				if (s1 == tail) tail = s0;
				
				/* Check for intersections */
				dx1 = s1->point[1]->x[0] - s1->point[0]->x[0];
				dy1 = s1->point[1]->x[1] - s1->point[0]->x[1];
				for (s0 = beam.nxt2; s0 != NULL; s0 = s0->nxt2)
				{
					dx0 = s0->point[1]->x[0] - s0->point[0]->x[0];
					dy0 = s0->point[1]->x[1] - s0->point[0]->x[1];
					den = dx0 * dy1 - dx1 * dy0;
					if (den > eps * eps || den < -eps * eps)
					{
						dx = s1->point[1]->x[0] - s0->point[0]->x[0];
						dy = s1->point[1]->x[1] - s0->point[0]->x[1];
						alpha = (dy1 * dx - dx1 * dy) / den;
						beta = (dx0 * dy - dy0 * dx) / den;
						if (alpha > 0 && beta > 0 && alpha < 1 && beta < 1)
						{
							alpha = (s0->point[0]->x[1] - s1->point[0]->x[1]) * ((dx0 * dx1)	/ den) + s1->point[0]->x[0] * ((dy1 * dx0) / den) - s0->point[0]->x[0] * ((dy0 * dx1) / den);
							beta = -(s0->point[0]->x[0] - s1->point[0]->x[0]) * ((dy0 * dy1)	/ den) - s1->point[0]->x[1] * ((dx1 * dy0) / den) + s0->point[0]->x[1] * ((dx0 * dy1) / den);
							if (p2 == NULL)
							{
								p2 = new_point(alpha, beta);
								point_last = p2;
							}
							else
							{
								point_last->nxt1 = new_point(alpha, beta);
								point_last = point_last->nxt1;
							}
						}
					}
				}
			}
		}
	}

	/* Add crossings to graph */
	p2 = mergesort_points(p2);
	point_first = merge_sorted_points(point_first, p2);

	/* Recalculate intersections between segments and points */
	result = 1;
	while (result > 0)
	{
		result = intersect_point_segment(point_first);
		for (p1 = point_first; p1 != NULL; p1 = p1->nxt1)
		{
			x1 = p1->x[0];
			p1->x[0] = p1->x[1];
			p1->x[1] = x1;
		}
		point_first = mergesort_points(point_first);
		for (s1 = segment_first; s1 != NULL; s1 = s1->nxt1) if (compare_points(s1->point[0], s1->point[1]) > 0)
		{
			p1 = s1->point[1];
			s1->point[1] = s1->point[0];
			s1->point[0] = p1;
		}

		result += intersect_point_segment(point_first);
		for (p1 = point_first; p1 != NULL; p1 = p1->nxt1)
		{
			x1 = p1->x[0];
			p1->x[0] = p1->x[1];
			p1->x[1] = x1;
		}
		point_first = mergesort_points(point_first);
		for (s1 = segment_first; s1 != NULL; s1 = s1->nxt1) if (compare_points(s1->point[0], s1->point[1]) > 0)
		{
			p1 = s1->point[1];
			s1->point[1] = s1->point[0];
			s1->point[0] = p1;
		}
	}

	/* Merge overlaping points and segments */
	tail = &beam;
	tail->nxt2 = NULL;
	p1 = point_first;
	p2 = point_first->nxt1;
	while (p1 != NULL)
	{
		if (p2 != NULL && compare_points(p1, p2) == 0)
		{
			for (i = 0; i < p2->max_seg && p2->segment[i] != NULL; ++i)
			{
				s1 = p2->segment[i];
				if (s1->point[1] != p1 && s1->point[0] != p1)
				{
					add_segment(p1, s1);
					if (s1->point[0] == p2) s1->point[0] = p1;
					else s1->point[1] = p1;
				}
				else
				{
					remove_segment(p1, s1);
					if (s1 == segment_first)
						segment_first = free_segment(segment_first);
					else
					{
						s0 = segment_first;
						while (s0->nxt1 != s1) s0 = s0->nxt1;
						s0->nxt1 = free_segment(s1);
					}
				}
			}
			p1->nxt1 = p2 = free_point(p2);
		}
		else
		{
			/* Update scanbeam */
			removed.nxt2 = NULL;
			for (i = 0; i < p1->max_seg && p1->segment[i] != NULL; ++i)
			{
				s1 = p1->segment[i];
				if (s1->point[0] == p1)
				{
					tail = tail->nxt2 = s1;
					tail->nxt2 = NULL;
				}
				else if (s1->point[1] == p1)
				{
					/* Search for duplicate edges only among the ones removed from the scanbeam so far in this iteration */
					s0 = &beam;
					while (s0->nxt2 != s1)
					{
						s0 = s0->nxt2;
					}
					s0->nxt2 = s1->nxt2;
					if (s1 == tail) tail = s0;

					s0 = &removed;
					while (s0->nxt2 != NULL && compare_adjacent_segments(s0->nxt2, s1) < 0) s0 = s0->nxt2;
					if (s0->nxt2 == NULL || compare_adjacent_segments(s0->nxt2, s1) > 0)
					{
						s1->nxt2 = s0->nxt2;
						s0->nxt2 = s1;
					}
					else
					{
						s0 = s0->nxt2;
						if (s0->point[0] == s1->point[0])
							for (j = num_polygons - 1; j >= 0; --j) s0->polygon[j] += s1->polygon[j];
						else
							for (j = num_polygons - 1; j >= 0; --j) s0->polygon[j] -= s1->polygon[j];
						remove_segment(s0->point[0], s1);
						remove_segment(s0->point[1], s1);
						--i;

						if (s1 == segment_first)
							segment_first = free_segment(segment_first);
						else
						{
							s0 = segment_first;
							while (s0->nxt1 != s1) s0 = s0->nxt1;
							s0->nxt1 = free_segment(s1);
						}
					}
				}
			}
			p1 = p2;
			p2 = (p2 != NULL ? p2->nxt1 : NULL);
		}
	}

#ifdef DEBUG
POINTS("\nAll pts:\n",point_first,NULL);
SEGMENTS("\nAll segs:\n",segment_first,NULL);
#endif

	/* Perform operation */
	counter = (long *) malloc(num_polygons * sizeof(long));
	phase = (short *) malloc(num_polygons * sizeof(short));
#if PY_MAJOR_VERSION >= 3
	obj3 = PyLong_FromLong(0L);
#else
	obj3 = PyInt_FromLong(0L);
#endif
	obj4 = PyDict_New();
	tail = &beam;
	for (i = 0; i < point_first->max_seg && point_first->segment[i] != NULL; ++i)
	{
		tail = tail->nxt2 = point_first->segment[i];
		//tail->output = 0;
	}
	tail->nxt2 = NULL;
	beam.nxt2 = sort_segments(beam.nxt2);
	for (p1 = point_first->nxt1; p1 != NULL; p1 = p1->nxt1)
	{
		memset(counter, 0, num_polygons * sizeof(long));
		memset(phase, 0, num_polygons * sizeof(short));
		last_result = 0;
#if DEBUG > 0
POINT("\nBeam at: ",p1);
#endif
		if (beam.nxt2 != NULL)
		{
#if DEBUG > 1
SEGMENTS_ALT("", beam.nxt2, NULL);
#endif
			for (s0 = beam.nxt2; s0->nxt2 != NULL; s0 = s0->nxt2)
			{
				obj1 = PyTuple_New(num_polygons);
				for (i = num_polygons - 1; i >= 0; --i)
				{
					if (s0->polygon[i] != 0)
					{
						if (counter[i] == 0)
						{
							phase[i] = s0->polygon[i];
							++counter[i];
						}
						else if (s0->polygon[i] * (long) phase[i] > 0) ++counter[i];
						else --counter[i];
					}
#if PY_MAJOR_VERSION >= 3
					PyTuple_SetItem(obj1, i, PyLong_FromLong(counter[i]));
#else
					PyTuple_SetItem(obj1, i, PyInt_FromLong(counter[i]));
#endif
				}
				if (PyDict_Contains(obj4, obj1) > 0)
				{
					obj2 = PyDict_GetItem(obj4, obj1);
					result = PyObject_IsTrue(obj2);
				}
				else
				{
					obj2 = PyObject_CallObject(operation, obj1);
					obj5 = PyObject_RichCompare(obj2, obj3, Py_GT);
					PyDict_SetItem(obj4, obj1, obj5);
					result = PyObject_IsTrue(obj5);
					Py_DECREF(obj5);
					Py_DECREF(obj2);
				}
				if (result != last_result)
				{
					if (result > 0) s0->output = 1;
					else s0->output = -1;
					last_result = result;
				}
				else s0->output = 0;
#if DEBUG > 0
printf("Op = %+d : ", result);
SEGMENT("", s0);
#endif
				Py_DECREF(obj1);
			}
			s0->output = (last_result > 0 ? -1 : 0);
#if DEBUG > 0
printf("Op last = %+d : ", last_result);
SEGMENT("", s0);
#endif
		}
		/* Update scanbeam */
		i = 0;
		while (i < p1->max_seg && p1->segment[i] != NULL)
		{
			s1 = p1->segment[i];
			if (s1->point[0] == p1)
			{
				s0 = &beam;
				while (s0->nxt2 != NULL && compare_segments(s1, s0->nxt2) > 0) s0 = s0->nxt2;
				//s1->output = 0;
				s1->nxt2 = s0->nxt2;
				s0->nxt2 = s1;
				++i;
			}
			else
			{
#if DEBUG > 0
SEGMENT("Test: ",s1);
#endif
				s0 = &beam;
				while (s0->nxt2 != s1) s0 = s0->nxt2;
				s0->nxt2 = s1->nxt2;
				if (s1->output != 0) ++i;
				else
				{
#if DEBUG > 0
SEGMENT("REM: ", s1);
#endif
					remove_segment(s1->point[0], s1);
					remove_segment(s1->point[1], s1);
					if (s1 == segment_first) segment_first = free_segment(segment_first);
					else
					{
						s0 = segment_first;
						while (s0->nxt1 != s1) s0 = s0->nxt1;
						s0->nxt1 = free_segment(s1);
					}
				}
			}
		}
	}
	Py_DECREF(obj3);
	Py_DECREF(obj4);
	free(counter);
	free(phase);

	if (segment_first == NULL)
	{
		free_point_list(point_first);
		Py_INCREF(Py_None);
		return Py_None;
	}

#ifdef DEBUG
SEGMENTS("\nResulting segs:\n",segment_first,NULL);
#endif

	/* Find resulting polygons */
	polygon_tree.inside = NULL;
	while (segment_first != NULL)
	{
		s0 = segment_first;
		segment_first = segment_first->nxt1;
		remove_segment(s0->point[0], s0);
		remove_segment(s0->point[1], s0);
		p0 = s0->point[0];
		point_last = p0;
		point_last->nxt2 = NULL;
		p1 = s0->point[1];
		left = s0->output;
		head.nxt3 = NULL;
#if DEBUG > 0
POINT("\nStart: ",p0);
#endif
		while (p1 != p0)
		{
#if DEBUG > 0
POINT("  + ",p1);
#endif
			if (p1->flag != 0)
			{
				for (p3 = NULL, point_last = p1; point_last != NULL; p3 = (p3 == NULL ? point_last : p3->nxt2), point_last = point_last->nxt2)
					if (point_last->flag != 0)
					{
						p2 = &head;
						while (p2->nxt3 != point_last) p2 = p2->nxt3;
						p2->nxt3 = point_last->nxt3;
						point_last->nxt3 = NULL;
						point_last->flag = 0;
						if (p3 != NULL) p3->nxt2 = duplicate_point(point_last);
					}
				g0 = new_polygon(duplicate_point(p1));
#if DEBUG > 0
printf("Remove [%lx]\n", g0);
POINTS_ALT("",g0->point,NULL);
#endif
				add_inside(g0, &polygon_tree);
				for (point_last = p0; point_last->nxt2 != p1; point_last = point_last->nxt2);
			}
			i = 0;
			while (i < p1->max_seg && p1->segment[i] != NULL &&
					(p1->segment[i]->point[0] != p1 || p1->segment[i]->output != left) && 
					(p1->segment[i]->point[1] != p1 || p1->segment[i]->output == left)) ++i;
			s1 = p1->segment[i];
			if (s1 == segment_first) segment_first = segment_first->nxt1;
			else
			{
				s0 = segment_first;
				while (s0->nxt1 != s1) s0 = s0->nxt1;
				s0->nxt1 = s1->nxt1;
			}
			remove_segment(s1->point[0], s1);
			remove_segment(s1->point[1], s1);
			if (p1->segment[0] != NULL)
			{
				p1->nxt3 = head.nxt3;
				head.nxt3 = p1;
				p1->flag = 1;
			}
			point_last = point_last->nxt2 = p1;
			point_last->nxt2 = NULL;
			if (s1->point[0] == p1) p1 = s1->point[1];
			else p1 = s1->point[0];
			free_segment(s1);
		}
		for (p3 = p0, point_last = p0->nxt2; point_last != NULL; p3 = p3->nxt2, point_last = point_last->nxt2)
			if (point_last->flag != 0) p3->nxt2 = duplicate_point(point_last);
		if (p0->segment[0] != NULL) p0 = duplicate_point(p0);
		g0 = new_polygon(p0);
#if DEBUG > 0
printf("Final [%lx]\n", g0);
POINTS_ALT("",g0->point,NULL);
#endif
		add_inside(g0, &polygon_tree);
		p2 = head.nxt3;
		while (p2 != NULL)
		{
			p3 = p2->nxt3;
			p2->nxt3 = NULL;
			p2->flag = 0;
			p2 = p3;
		}
	}

#ifdef DEBUG
printf("\nAll polygons:\n");
POLYGONS(0, polygon_tree.inside);
#endif

	/* Connect holes to their outer polygons */
	polygon_last = polygon_tree.inside;
	while (polygon_last->nxt != NULL) polygon_last = polygon_last->nxt;
	for (g0 = polygon_tree.inside; g0 != NULL; g0 = g0->nxt)
	{
#if DEBUG > 0
printf("\nPiercing [%lx]\n", g0);
fflush(stdout);
#endif
		hole_prev1 = NULL;
		g1 = g0->inside;
		while (g1 != NULL)
		{
			/* Polygons inside holes are not holes */
			if (g1->inside != NULL)
			{
				polygon_last->nxt = g1->inside;
				g1->inside = NULL;
				while (polygon_last->nxt != NULL) polygon_last = polygon_last->nxt;
			}
			/* Check for existing connections with outer polygon */
			left = 1;
			for (p0 = g0->point; p0 != NULL && left > 0; p0 = p0->nxt2)
			{
				if (p0->nxt3 != NULL || (p0->nxt1 != NULL && p0->nxt1->nxt3 == p0))
				{
					for (p1 = g1->point; p1 != NULL && left > 0; p1 = p1->nxt2)
					{
						if (p1->nxt3 == p0 || p0->nxt3 == p1 || (p1->nxt3 != NULL && p0->nxt3 == p1->nxt3))
						{
							p2 = (p1->nxt2 != NULL ? p1->nxt2 : g1->point);
							p1->nxt2 = p0->nxt2;
							p4 = p1;
							while (p2 != p1)
							{
								p0->nxt2 = p2;
								p3 = p2->nxt2;
								p2->nxt2 = p4;
								p4 = p2;
								p2 = (p3 != NULL ? p3 : g1->point);
							}
							if (hole_prev1 == NULL) g0->inside = g1->nxt;
							else hole_prev1->nxt = g1->nxt;
							free(g1);
							g1 = hole_prev1;
							left = 0;
						}
					}
				}
			}
			if (left > 0)
			{
				/* Check for existing connections with other holes */ 
				hole_prev2 = g1;
				g2 = g1->nxt;
				while (g2 != NULL)
				{
					for (p1 = g1->point; p1 != NULL && left > 0; p1 = p1->nxt2)
					{
						if (p1->nxt3 != NULL || (p1->nxt1 != NULL && p1->nxt1->nxt3 == p1))
						{
							for (p2 = g2->point; p2 != NULL && left > 0; p2 = p2->nxt2)
							{
								if (p2->nxt3 == p1 || p1->nxt3 == p2 || (p1->nxt3 != NULL && p1->nxt3 == p2->nxt3))
								{
									if (g2->inside != NULL) polygon_last->nxt = g2->inside;
									while (polygon_last->nxt != NULL) polygon_last = polygon_last->nxt;
									p3 = p2;
									while (p3->nxt2 != NULL) p3 = p3->nxt2;
									p3->nxt2 = g2->point;
									p3 = p1->nxt2;
									p1->nxt2 = p2->nxt2;
									p2->nxt2 = p3;
									hole_prev2->nxt = g2->nxt;
									free(g2);
									g2 = hole_prev2;
									left = 0;
								}
							}
						}
					}
					hole_prev2 = g2;
					g2 = g2->nxt;
				}
			}
			hole_prev1 = g1;
			g1 = (g1 == NULL ? g0->inside : g1->nxt);
		}
		/* Connect to outer polygon */
		g1 = mergesort_polygons(g0->inside);
		g0->inside = NULL;
		while (g1 != NULL)
		{
			p1 = g1->point;
			for (p2 = p1->nxt2; p2 != NULL; p2 = p2->nxt2) if (compare_points(p2, p1) < 0) p1 = p2;
			/* Intersect outer polygon */
			i = 0;
			for (p2 = g0->point, p3 = g0->point->nxt2; p2 != NULL; p2 = p2->nxt2, p3 = p3->nxt2)
			{
				if (p3 == NULL) p3 = g0->point;
				beta = p1->x[0];
				if (p2->x[1] >= p1->x[1] - eps && p2->x[1] <= p1->x[1] + eps)
					beta = p2->x[0];
				else if ((p2->x[1] > p1->x[1] - eps && p3->x[1] < p1->x[1] - 2 * eps) || (p2->x[1] < p1->x[1] + eps && p3->x[1] > p1->x[1] + 2 * eps))
					beta = p2->x[0] + (p3->x[0] - p2->x[0]) * ((p1->x[1] - p2->x[1]) / (p3->x[1] - p2->x[1]));
				if (beta < p1->x[0] && (i == 0 || beta > alpha))
				{
					i = 1;
					p0 = p2;
					alpha = beta;
				}
			}
			beta = p0->x[1] - p1->x[1];
			if (beta > eps || beta < -eps)
			{
				/* Create new vertex */
				p2 = new_point(alpha, p1->x[1]);
				p2->nxt1 = point_first;
				point_first = p2;
				p2->nxt2 = p0->nxt2;
				p0->nxt2 = p2;
				p0 = p2;
			}
			/* Connect at p0 */
			p2 = duplicate_point(p1);
			p1 = p1->nxt2 = p4 = duplicate_point(p0);
			while (p2 != p1)
			{
				p0->nxt2 = p2;
				p3 = p2->nxt2;
				p2->nxt2 = p4;
				p4 = p2;
				p2 = (p3 != NULL ? p3 : g1->point);
			}
			g2 = g1;
			g1 = g1->nxt;
			free(g2);
		}
	}

#ifdef DEBUG
printf("\nPierced polygons:\n");
POLYGONS(0, polygon_tree.inside);
#endif

	for (i = 0, g0 = polygon_tree.inside; g0 != NULL; ++i, g0 = g0->nxt);
	if ((obj1 = PyTuple_New(i)) == NULL)
	{
		free_polygon_list(polygon_tree.inside);
		free_point_list(point_first);
		return NULL;
	}
	for (i = 0, g0 = polygon_tree.inside; g0 != NULL; ++i, g0 = g0->nxt)
	{
		for (j = 0, p0 = g0->point; p0 != NULL; ++j, p0 = p0->nxt2);
		if ((obj2 = PyTuple_New(j)) == NULL)
		{
			Py_DECREF(obj1);
			free_polygon_list(polygon_tree.inside);
			free_point_list(point_first);
			return NULL;
		}
		for (j = 0, p0 = g0->point; p0 != NULL; ++j, p0 = p0->nxt2)
		{
			if ((obj3 = PyTuple_New(2)) == NULL)
			{
				Py_DECREF(obj1);
				Py_DECREF(obj2);
				free_polygon_list(polygon_tree.inside);
				free_point_list(point_first);
				return NULL;
			}
			if ((obj4 = PyFloat_FromDouble(p0->x[0])) == NULL)
			{
				Py_DECREF(obj1);
				Py_DECREF(obj2);
				Py_DECREF(obj3);
				free_polygon_list(polygon_tree.inside);
				free_point_list(point_first);
				return NULL;
			}
			PyTuple_SET_ITEM(obj3, 0, obj4);
			if ((obj4 = PyFloat_FromDouble(p0->x[1])) == NULL)
			{
				Py_DECREF(obj1);
				Py_DECREF(obj2);
				Py_DECREF(obj3);
				free_polygon_list(polygon_tree.inside);
				free_point_list(point_first);
				return NULL;
			}
			PyTuple_SET_ITEM(obj3, 1, obj4);
			PyTuple_SET_ITEM(obj2, j, obj3);
		}
		PyTuple_SET_ITEM(obj1, i, obj2);
	}

	free_polygon_list(polygon_tree.inside);
	free_point_list(point_first);

	return obj1;
}

static const char doc[] = "Boolext is a Python C extension that implements boolean operations\n\
between polygons (polygon clipping).\n\n\
It is distributed as part of the *gdspy* module to help the creation of\n\
complex structures in the GDSII stream format, but the function\n\
interfaces should allow its usage in more general situations.";

static PyMethodDef booleanMethods[] = {
	{"clip", clip, METH_VARARGS,\
"Perform a boolean operation (clipping) on a set of polygons.\n\n\
Parameters\n\
----------\n\
polygons : list of array-like[N][2]\n\
        List of polygons. Each polygon is an array-like[N][2] object\n\
        with the coordinates of the vertices of the polygon.\n\
operation : function\n\
        This function should accept N of variables as arguments and\n\
        output an integer or boolean representing the desired operation\n\
        to be performed on the polygons. Each input variable\n\
        corresponds to one polygon of ``polygons``.\n\
eps : positive number\n\
        Small number to be used as a tolerance for intersection and overlap\n\
        calculations.\
Returns\n\
-------\n\
out : list of array-like[N][2]\n\
        List of polygons resulting from the boolean operation."},
	{NULL, NULL, 0, NULL}
};


#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef booleanModule = {
	PyModuleDef_HEAD_INIT,
	"boolext",
	doc,
	-1,
	booleanMethods
};

PyMODINIT_FUNC PyInit_boolext(void)
{
	return PyModule_Create(&booleanModule);
}
#else
PyMODINIT_FUNC initboolext(void)
{
	(void) Py_InitModule3("boolext", booleanMethods, doc);
}
#endif
