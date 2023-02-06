/*
 * This program is intended to analyze the strain energy density of elastic
 * networks that have been relaxed after the imposition of displacements at
 * their boundaries. The program accepts as input the initial and final states
 * of a network, as well as the bonds' stretching and bending stiffnesses.
 * The convex hull of the network's point set is then partitioned into a
 * Voronoi tessellation, the generating points of which are the vertices
 * of the network. The energy density is defined cell-wise, by assigning
 * the bending energy of each bond triple to the cell whose vertex is the
 * hinge of that triple, and each bond's stretching energy to those cells
 * through which it passes, in proportion with the fraction of the bond
 * in each cell. After each cell's strain energy is computed, its strain
 * energy density is computed by dividing by cell area.
 *
 * Following determination of cells' energy densities, the Cellular Fourier
 * Transform for the strain energy density is computed, according to the
 * method described by Frueleux and Boudaoud in "Cellular Fourier analysis for 
 * geometrically disordered materials", Phys. Rev. Research 3, p. 023036, 2021.
 */

#include <boost/geometry.hpp>
#include <boost/polygon/voronoi.hpp>
#include <boost/polygon/voronoi_builder.hpp>
#include <boost/polygon/voronoi_diagram.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
//#include <boost/geometry/index/rtree.hpp>
#include <tuple>
#include <string>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <boost/regex.hpp>
#include <stdio.h>
#include <cmath>
//#include "src/triangle.h"
#include <map>

using namespace std;
using boost::polygon::voronoi_builder;
using boost::polygon::voronoi_diagram;

//Shorthand for important namespaces
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

//Shorthand for data types used in spatial indexing
typedef bg::model::point<double, 2, bg::cs::cartesian> point;
typedef bg::model::box<point> box;
typedef std::pair<point, unsigned> value;

//Data structures for representing bonds in the network and 
//interacting triples.
typedef tuple<double, double> vec2d;
typedef tuple<int, int, double, char> edge_datum;
typedef tuple<int, int, int, char> triple_datum;

//Structure for splitting strings according to a user-defined delimiter
typedef boost::tokenizer<boost::char_separator<char>> tokenizer;

//Constants for converting between integer and floating point formats for
//coordinates
#define INT_FLOOR 0xF8000000
#define INT_CEIL 0x07FFFFFF
#define DIFLOOR -134217728
#define DIRANGE 268435455
#define MIN_SLOPE 10e-20
#define MAX_SLOPE 10e20
#define FLOAT_TOL 1e-8

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

//Check whether a point is in bounds
//#define in_xy_bounds(a, b, c, d, e, f) (a > c - FLOAT_TOL && a < e + FLOAT_TOL && b > d - FLOAT_TOL && b < f + FLOAT_TOL) ? true : false
#define in_y_bounds(a, b, c) (a > b - FLOAT_TOL && a < c + FLOAT_TOL) ? true : false

//Strucutres for points with integer coordinates and line segments bounded by 
//such points
struct IPoint {
    int x, y;

    IPoint(int myx, int myy) : x(myx), y(myy) {}
};

struct Segment {
    IPoint p0, p1;

    Segment(int x1, int y1, int x2, int y2) : p0(x1, y1), p1(x2, y2) {}
};

//Double-precision floating-point coordinate type
struct Point {

    double x, y;

    Point(){
	    x = 0;
	    y = 0;
    }

    Point(double myx, double myy) : x(myx), y(myy) {}

    //Define equality as coordinates differing by less than a given
    //floating point arithmetic tolerance. While in principle this
    //places all points in the same equivalence class, it should work
    //in practice as long as the precision of floating point arithmetic
    //is smaller than the minimum spacing between points.
    bool operator == (const Point& point) const{
        return(abs(point.x - x) < FLOAT_TOL && abs(point.y - y) < FLOAT_TOL);
    }

    //Stipulate that points be sorted first according to their x 
    //coordinates, then according to their y coordiantes. Include a
    //tolerance so that almost bit-for-bit equal pairs of coordinates
    //are treated as equal.
    friend bool operator < (const Point &p1, const Point &p2){
        if(p2.x > p1.x + FLOAT_TOL) return true;
        else if(p1.x > p2.x + FLOAT_TOL) return false;
        else if(p2.y > p1.y + FLOAT_TOL) return true;
        else return false;
    }
};

//A means of sorting points according to their distance from some common point,
//in order of least to most distant.
struct PSort {

    Point base;

    PSort(Point mybase) : base(mybase) {}

    bool operator() (Point p1, Point p2){
                double sqdist1, sqdist2;

                sqdist1 = (p1.x-base.x)*(p1.x-base.x)+(p1.y-base.y)*(p1.y-base.y);
                sqdist2 = (p2.x-base.x)*(p2.x-base.x)+(p2.y-base.y)*(p2.y-base.y);
                return sqdist1 < sqdist2;
        }
};

//Boilerplate code to ensure IPoint and Segment are compatible with Boost
//prototypes
namespace boost {
namespace polygon {

template <>
struct geometry_concept<IPoint> {
        typedef point_concept type;
};

template <>
struct point_traits<IPoint> {
        typedef int coordinate_type;

        static inline coordinate_type get(const IPoint& point, orientation_2d orient) {
                return (orient == HORIZONTAL) ? point.x : point.y;
        }
};

template <>
struct geometry_concept<Segment> {
        typedef segment_concept type;
};

template <>
struct segment_traits<Segment> {
        typedef int coordinate_type;
        typedef IPoint point_type;

        static inline point_type get(const Segment& segment, direction_1d dir) {
                return dir.to_int() ? segment.p1 : segment.p0;
        }
};
}
}

//Read the reference state of a network
bool read_network(ifstream &input, vector<Point> &points, vector<edge_datum> &edges, double &minx, double &maxx, double &miny, double &maxy, vector<int> &left_set){

    string nextline;
    int num_read, num_points, num_edges, idx1, idx2;
    double x, y, dx, dy, offset;
    char mult;

    minx = FLT_MAX;
    miny = FLT_MAX;
    maxy = FLT_MIN;

    //Attempt to read the header for the network description file. If too little
    //information or invalid data are found, close the input stream and report
    //a failure to read a network.
    getline(input, nextline);
    num_read = sscanf(nextline.c_str(), "%d %d %lf", &num_points, &num_edges, &offset);
    if(num_read < 3 || offset < 0 || num_points <= 0 || num_edges < 0){
        cerr << "Reading failed. The provided file has an invalid header.\n";
        input.close();
        return false;
    }

    //If a valid header was read, attempt to read the specified number of points
    while(! input.eof() && points.size() < num_points){
        getline(input, nextline);
        num_read = sscanf(nextline.c_str(), "%lf %lf", &x, &y);
        if(num_read == 2){
            points.push_back(Point(x, y));
            if(x < minx - FLOAT_TOL){
	        minx = x;
    	        left_set.clear();
		left_set.push_back(points.size() - 1);
            }
	    else if(x < minx + FLOAT_TOL) left_set.push_back(points.size()-1);

            if(y < miny) miny = y;
            if(y > maxy) maxy = y;
        }
    }
    maxx = minx + offset;

    //If the program has reached the end of the input stream, too little
    //information has been provided, and the file is invalid.
    if(input.eof()){
        cerr << "Reading failed. Too few points were read.\n";
        input.close();
        return false;
    }


    while(!input.eof() && edges.size() < num_edges){
        getline(input, nextline);
        num_read = sscanf(nextline.c_str(),"%d %d %hhd\n", &idx1, &idx2, &mult);

        //Ensure two point indices and an offset are given. If so, and the point
        //indices are in bounds, find the length of the edge, and add it to
        //the list.
        if(num_read == 3){
            if(min(idx1, idx2) > -1 && max(idx1, idx2) < num_points){
                dx = points[idx1].x - points[idx2].x;
                dy = points[idx1].y - points[idx2].y;
                edges.push_back(make_tuple(idx1, idx2, sqrt(dx*dx+dy*dy), mult));
            }
            else{
                cerr << "Reading failed due to an out-of-bounds point index.\n";
                input.close();
                return false;
            }
        }
    }

    input.close();

    if(edges.size() < num_edges) cerr << "Too few edges were found.\n";

    return edges.size() == num_edges;
}

//Find all triples in the network
vector<triple_datum> get_triples(vector<Point> points, vector<edge_datum> edges, double offset){

    //Map from each point to the points with which it shares a bond
    map<int, vector<tuple<int, char>>> adjacency_lists;
    vector<triple_datum> triples;
    int idx1, idx2, idx3, i, j;
    double x, y, dx2, dy2, dx3, dy3, dot, mag_sq;
    tuple<int, char> pair1, pair2;
    vector<tuple<int, char>> neighbors;

    for(edge_datum next_edge : edges){
	idx1 = get<0>(next_edge);
	idx2 = get<1>(next_edge);

	if(adjacency_lists.find(idx1) == adjacency_lists.end()){
	    adjacency_lists.insert(make_pair(idx1, vector<tuple<int, char>>()));
	}
	if(adjacency_lists.find(idx2) == adjacency_lists.end()){
	    adjacency_lists.insert(make_pair(idx2, vector<tuple<int, char>>()));
	}

	adjacency_lists[idx1].push_back(make_tuple(idx2, get<2>(next_edge)));
	adjacency_lists[idx2].push_back(make_tuple(idx1, -get<2>(next_edge)));
    }

    //Once the adjacency lists are built, identify all antiparallel bonds
    //sharing a common vertex. Whenever such a pair of bonds is identified,
    //create a triple, represented by the integer index of the common point,
    //or "hinge", followed by the integer indices of the two outside points,
    //followed by the offset that should be added for the bond from the hinge
    //to the first outside point.
    for(auto iter=adjacency_lists.begin();iter != adjacency_lists.end();iter++){
        
	idx1 = iter->first;
	neighbors = iter->second;
	x = points[idx1].x;
	y = points[idx1].y;

	for(i = 0; i < neighbors.size(); i ++){
	    pair1 = neighbors[i];
            idx2 = get<0>(pair1);
	    dx2 = points[idx2].x + offset*get<1>(pair1) - x;
	    dy2 = points[idx2].y - y;

	    for(j = i+1; j < neighbors.size(); j++){
                pair2 = neighbors[j];
		idx3 = get<0>(pair2);

		dx3 = points[idx3].x + offset*get<1>(pair2) - x;
		dy3 = points[idx3].y - y;

		dot = dx2*dx3 + dy2*dy3;
		mag_sq = (dx2*dx2 + dy2*dy2)*(dx3*dx3 + dy3*dy3);

		//If the dot product of the vectors pointing from the hinge
		//to either end point is negative, and the square of this
		//dot product is equal to the product of the square of the
		//magnitudes of the two displacement vectors, the two 
		//displacement vectors are antiparallel, and a valid triple has
		//been identified.
		if(dot < 0 && abs(dot*dot  - mag_sq) < FLOAT_TOL){
		    triples.push_back(make_tuple(idx1,idx2,idx3,get<1>(pair1)));
		}
	    }
	}
    }

    return triples;
}

//Read a displacement field
vector<vec2d> read_displacements(ifstream &input){

    int num_read;
    double x, y;
    string nextline;

    vector<vec2d> disps;

    while(! input.eof()){
        getline(input, nextline);
        num_read = sscanf(nextline.c_str(), "%lf %lf\n", &x, &y);
    if(num_read == 2){
        disps.push_back(make_tuple(x, y));
    }
    };

    input.close();

    return disps;
}

//Snap points with floating point coordinates to an integer grid to enable
//use of the Boost library's implementation of Voronoi tessellation
vector<IPoint> snap_to_grid(vector<Point> in, double xmin, double xrange, double ymin, double yrange){

        vector<IPoint> out_points;
        int int_x, int_y;

        for(Point in_point : in){
                int_x = (int) ((in_point.x - xmin) / xrange * (INT_CEIL-INT_FLOOR)) + INT_FLOOR;
                int_y = (int) ((in_point.y - ymin) / yrange * (INT_CEIL-INT_FLOOR)) + INT_FLOOR;
                out_points.push_back(IPoint(int_x, int_y));
        }

        return out_points;
}

//Boost's Voronoi diagram implementation reorders input points and their
//associated Voronoi cells. This function determines the cycles of the
//permutation mapping input points' initial sequence positions to their
//sequence positions according to Boost's sorting, and reverses the
//permutation.
/*void reverse_permutation(vector<vector<int>> &polygons, vector<int> positions){

    vector<vector<int>> cycles;
    set<int> unvisited;
    vector<int> temp;
    int start, curr, iter;

    for(int iter = 0; iter < polygons.size(); iter++) unvisited.insert(iter);

    //Enumerate all cycles of the permutation
    while(unvisited.size() > 0){

	vector<int> next_cycle;

        start = *(unvisited.begin());
	curr = start;
	do{
            unvisited.erase(curr);
	    next_cycle.push_back(curr);
	    curr = positions[curr];
	}while(curr != start);

        cycles.push_back(next_cycle);
    }

    //Undo the permutation by traversing each cycle, transposing successive
    //pairs of elements.
    for(vector<int> next_cycle : cycles){
        for(iter = 0; iter < next_cycle.size() - 1; iter++){
	    temp = polygons[next_cycle[iter + 1]];
	    polygons[next_cycle[iter + 1]] = polygons[next_cycle[iter]];
	    polygons[next_cycle[iter]] = temp;
	}
    }

}*/

//Compute the Voronoi diagram for a point set. The input is a set of points
//with integer coordinates, and the output is a set of polygons specified
//as a set of points with double precision floating point coordinates.
void get_voronoi_cells(vector<Point> dpts, set<int> keep_indices, vector<Point> &outPoints, vector<vector<int>> &polygons, double minx, double xrange, double miny, double yrange){

    //Point p0, p1;
    map<Point, int> pmap;
    int p_index = 0, idx1, idx2, curr_idx, neigh_idx, pos = 0;
    double x0, y0, x1, y1, dx, dy, len, xc = 0, yc = 0, rx1, ry1, rx2,ry2,cross;
    vector<Segment> segments;
    double long_len = sqrt(xrange*xrange + yrange*yrange);
    vector<IPoint> int_pts;
    //vector<int> positions(keep_indices.size(), 0);
    map<int, int> index_map;

    for(int iter = 0; iter < keep_indices.size(); iter++){
        polygons.push_back(vector<int>());
    }

    curr_idx = 0;
    for(auto it = keep_indices.begin(); it != keep_indices.end(); it++){
        index_map.insert(make_pair(*it, curr_idx++));
    }

    //Find the mean x and y coordinates for the point set
    for(Point next_point : dpts){
        xc += next_point.x;
	yc += next_point.y;
    }
    xc /= dpts.size();
    yc /= dpts.size();

    int_pts = snap_to_grid(dpts, minx, xrange, miny, yrange);

    //Create a Voronoi diagram of the point packing
    voronoi_diagram<double> vd;
    construct_voronoi(int_pts.begin(), int_pts.end(), segments.begin(), segments.end(), &vd);

    //Iterate over the cells in the voronoi diagram, and enumerate their
    //vertices. When an infinite edge is encountered, replace that
    //edge with a long, but finite, edge with the same orientation.
    for(voronoi_diagram<double>::const_cell_iterator it = vd.cells().begin(); it != vd.cells().end(); ++it){

        curr_idx = (*it).source_index();
        if(keep_indices.find(curr_idx) == keep_indices.end()) continue;
	//positions[curr_idx] = pos++;
	//cout << "Current index: " << curr_idx << ", Size: " << keep_indices.size() << "\n";

        const voronoi_diagram<double>::cell_type &cell = *it;
        const voronoi_diagram<double>::edge_type *edge = cell.incident_edge();

        //vector<int> vertex_list;
        //cout << "Start of the do loop.\n";

        do{
            if(edge->is_primary()){
                if(edge->is_finite()){
                    x0 = (edge->vertex0()->x()-DIFLOOR)/DIRANGE*xrange + minx;
                    y0 = (edge->vertex0()->y()-DIFLOOR)/DIRANGE*yrange + miny;
                    Point p0 = Point(x0, y0);

                    if(pmap.find(p0) == pmap.end()){
                        pmap.insert(make_pair(p0, p_index++));
                        outPoints.push_back(p0);
                    }

                    //vertex_list.push_back(pmap[p0]);
                    polygons[index_map[curr_idx]].push_back(pmap[p0]);
                }

                else{
                    neigh_idx = edge->twin()->cell()->source_index();
                    dx = dpts[neigh_idx].x - dpts[curr_idx].x;
                    dy = dpts[neigh_idx].y - dpts[curr_idx].y;
                    len = sqrt(dx*dx + dy*dy);

		    rx1 = dpts[curr_idx].x - xc;
                    ry1 = dpts[curr_idx].y - yc;
		    rx2 = dpts[neigh_idx].x - xc;
		    ry2 = dpts[neigh_idx].y - yc;
		    cross = rx1*ry2 - rx2*ry1;
		    cross /= abs(cross);

                    if(edge->vertex0()){
                        x0 = (edge->vertex0()->x()-DIFLOOR)/DIRANGE*xrange + minx;
                        y0 = (edge->vertex0()->y()-DIFLOOR)/DIRANGE*yrange + miny;
                        Point p0 = Point(x0, y0);

                        if(pmap.find(p0) == pmap.end()){
                            pmap.insert(make_pair(p0, p_index++));
                            outPoints.push_back(p0);
                        }

                        //vertex_list.push_back(pmap[p0]);
			polygons[index_map[curr_idx]].push_back(pmap[p0]);

                        x1 = x0 + long_len * dy / len * cross;
                        y1 = y0 - long_len * dx / len * cross;
                        Point p1 = Point(x1, y1);

                        if(pmap.find(p1) == pmap.end()){
                            pmap.insert(make_pair(p1, p_index++));
                            outPoints.push_back(p1);
                        }

                        //vertex_list.push_back(pmap[p1]);
                        polygons[index_map[curr_idx]].push_back(pmap[p1]);
                    }      

                    else{
                        x1 = (edge->vertex1()->x()-DIFLOOR)/DIRANGE*xrange + minx;
                        y1 = (edge->vertex1()->y()-DIFLOOR)/DIRANGE*yrange + miny;

                        x0 = x1 + long_len * dy / len * cross;
                        y0 = y1 - long_len * dx / len * cross;
                        Point p0 = Point(x0, y0);

                        if(pmap.find(p0) == pmap.end()){
                            pmap.insert(make_pair(p0, p_index++));
                            outPoints.push_back(p0);
                        }

                        //vertex_list.push_back(pmap[p0]);
                        polygons[index_map[curr_idx]].push_back(pmap[p0]);
                    }
                }
	    }

	    edge = edge->next();
        }while(edge != cell.incident_edge());

	//cout << "End of the do loop.\n";
    }

    //cout << "Made it here.\n";

}

//Given two points, find any points at which the line segment between them
//passes through a given y coordinate.
vector<Point> find_y_intersections(Point p1, Point p2, double minx, double miny, double maxx, double maxy){

    double low, high, left, right, x_int, y_int, slope;
    vector<Point> intersections;

    low = min(p1.y, p2.y);
    high = max(p1.y, p2.y);
    left = min(p1.x, p2.x);
    right = max(p1.x, p2.x);

    if(abs(p2.x - p1.x) > MIN_SLOPE){
        slope = abs(p2.y - p1.y) > MIN_SLOPE ? (p2.y - p1.y) / (p2.x - p1.x) : MIN_SLOPE;
    }
    else{
        slope = copysign(MAX_SLOPE, p2.y - p1.y);
    }

    //Check for intersection with the bottom edge
    if(miny >= low && miny <= high){
        x_int = (miny - p1.y) / slope + p1.x;
        //if(x_int >= minx && x_int <= maxx){
            Point new_point(x_int, miny);
            if(!(new_point == p1 || new_point == p2)){
                intersections.push_back(new_point);
            }
        //}
    }
    
    //Check for intersection with the top edge
    if(maxy >= low && maxy <= high){
        x_int = (maxy - p1.y) / slope + p1.x;
        //if(x_int >= minx && x_int <= maxx){
            Point new_point(x_int, maxy);
            if(!(new_point == p1 || new_point == p2)){
                intersections.push_back(new_point);
            }
        //}
    }

    /*
    //Check for intersection with the left edge
    if(minx >= left && minx <= right){
        y_int = slope*(p1.x - minx) + p1.y;
        if(y_int >= miny && y_int <= maxy){
            Point new_point(minx, y_int);
            if(!(new_point == p1 || new_point == p2)){
                intersections.push_back(new_point);
            }
        }
    }
   
    //Check for intersection with the right edge
    if(maxx >= left && maxx <= right){
        y_int = slope*(p1.x - maxx) + p1.y;
        if(y_int >= miny && y_int <= maxy){
            Point new_point(maxx, y_int);
            if(!(new_point == p1 || new_point == p2)){
                intersections.push_back(new_point);
            }
        }
    }
    */

    //Sort intersection points in order of least to greatest distance from p1
    if(intersections.size() > 1){
        sort(intersections.begin(), intersections.end(), PSort(p1));
    }

    return intersections;
}

//Trim Voronoi cells to fit within a set of y bounds. It is assumed
//that cells should only be trimmed to comply with y bounds, but not x bounds,
//as x bounds have been imposed by constructing the Voronoi diagram with
//multiple translated copies of the original point set.
void get_trimmed_polygons(vector<Point> in_vertices, vector<vector<int>> in_polygons, vector<Point> &out_vertices, vector<vector<int>> &out_polygons, double minx, double miny, double maxx, double maxy){

    map<Point, int> pmap;
    int pindex = 0;
    Point curr_p, next_p;
    bool curr_in, next_in;

    for(vector<int> next_poly : in_polygons){

        //if(next_poly.size() < 3) continue;

        curr_p = in_vertices[next_poly[0]];
        vector<int> new_polygon;

        curr_in = in_y_bounds(curr_p.y, miny, maxy);

        for(int iter = 0; iter < next_poly.size(); iter ++){
            next_p = in_vertices[next_poly[(iter + 1) % next_poly.size()]];

            next_in = in_y_bounds(next_p.y, miny, maxy);

            //If the current vertex is in bounds, it should be added to the 
            //trimmed polygon.
            if(curr_in){
                if(pmap.find(curr_p) == pmap.end()){
                    pmap.insert(make_pair(curr_p, pindex));
                    out_vertices.push_back(curr_p);
                    pindex ++;
                }

                new_polygon.push_back(pmap[curr_p]);
            }

            //If either the current or the next vertex is not within y bounds
            //those points along the segment from the current to the next
            //vertex where this segment intersects a boundary should be
            //included in the trimmed polygon.
	    if(!curr_in || !next_in){
                for(Point next_cross : find_y_intersections(curr_p, next_p, minx, miny, maxx, maxy)){
                    if(pmap.find(next_cross) == pmap.end()){
                        pmap.insert(make_pair(next_cross, pindex));
                        out_vertices.push_back(next_cross);
                        pindex ++;
                    }

                    new_polygon.push_back(pmap[next_cross]);
                }
	    }

            curr_in = next_in;
	    curr_p = next_p;
        }

        if(new_polygon.size() > 0){
            out_polygons.push_back(new_polygon);
        }
    }
}

//Break a polygon, specified as a set of ordered vertices, into a triangulation
//for use in numerical integration.

//Given a reference network configuration, a deformed configuration, and a set
//of mechanical attributes, find the strain energy associated with each node in
//the network. 

//Parse a job file. The file should provide the following information:
//    -A pattern for network file names
//    -A pattern for displacement file names
//    -A stretching stiffness, a bending stiffness, and a shear modulus
//    -A kernel type (either exponential or Gaussian), a kernel width, and
//    a length cutoff for evaluating integrals for cells' kernels
//    -The fraction of discrete Laplacian eigenmodes to obtain
//All blank lines, and lines beginning with a "#" will be ignored. A return type
//of true indicates all necessary information was obtained, while a return type
//of false indicates some information was missing.
bool parse_input_args(string in_name, string &net_pat, string &disp_pat, double &ks, double &kb, double &mu, char &ktype, double &width, double &cutoff, double &frac){

    ifstream infile;
    vector<string> lines;
    int num_read;
    boost::regex blank_regex("^\\s*$");
    string nextline;

    //Read all lines of the input file, and scrub blank lines and comments
    infile.open(in_name, ifstream::in);
    if(! infile.is_open()){
        cerr << "The job file could not be opened for reading.\n";
        return false;
    }

    while(! infile.eof()){
        getline(infile, nextline);
        if((! regex_match(nextline, blank_regex)) && !( nextline[0] == '#')){
            lines.push_back(nextline);
        }
    }

    infile.close();

    //If too few lines were read, report failure.
    if(lines.size() < 5){
        cerr << "Too little information was supplied.\n";
        return false;
    }

    //Obtain file name patterns
    net_pat = lines[0];
    disp_pat = lines[1];
    
    //Obtain mechanical attributes
    num_read = sscanf(lines[2].c_str(), "%lf %lf %lf", &ks, &kb, &mu);
    if(num_read < 3){
        cerr << "To few mechanical attributes were supplied.\n";
        return false;
    }

    //Obtain kernel information
    num_read = sscanf(lines[3].c_str(), "%c %lf %lf", &ktype, &width, &cutoff);
    if(num_read < 3){
        cerr << "Too little kernel information was obtained.\n";
        return false;
    }

    //Obtain the eigenmode fraction
    num_read = sscanf(lines[4].c_str(), "%lf", &frac);
    if(num_read < 1){
        cerr << "No eigen mode fraction could be obtained.\n";
        return false;
    }

    return true;
}

vector<string> split(string input, char delim){
    vector<string> result;
    size_t start, iter, len;
    bool reading_token = false;

    for(iter = 0; iter < input.size(); iter++){
        if(delim != input[iter]){
            if(!reading_token){ 
                reading_token = true;
                start = iter;
            }
        }

        else{
            if(reading_token){
                reading_token = false;
                result.push_back(input.substr(start, iter - start));
            }
        }
    }

    if(reading_token) result.push_back(input.substr(start, iter - start));


    return result;
}

//Append a directory name to a file within that directory
string get_full_name(string name, string dir){
	ostringstream full_name;
	full_name << dir << name;
	return full_name.str();
}

//Find a list of files matching a pattern
vector<string> get_file_list(string pattern){

    string dir_name, base_name, name;
    boost::filesystem::path p;
    vector<string> tokens, matches;
    ostringstream oss;
    boost::regex expr;
    boost::smatch result;
    int iter;

    if(pattern.compare("") == 0) return matches;

    //Split the response into tokens, delimited by a "/"
    tokens = split(pattern, '/');
    base_name = *tokens.rbegin();

    if(tokens.size() > 1){
        if(pattern[0] == '/') oss << "/";
        for(iter = 0; iter < tokens.size() - 1; iter ++){
            oss << tokens[iter] << "/";
        }
        dir_name = oss.str();
    }
    else{
        dir_name = "./";
    }
    p = boost::filesystem::path(dir_name);

    //Create the regular expression
    expr = boost::regex{base_name};

    if(exists(p)){
        for(auto i=boost::filesystem::directory_iterator(p); i!=boost::filesystem::directory_iterator(); i++){
            if(! boost::filesystem::is_directory(i->path())){
                name = i->path().filename().string();
                if(boost::regex_search(name, result, expr)){
                    matches.push_back(get_full_name(result[0], dir_name));
                }
            }
        }
    }

    //Sort the results
    sort(matches.begin(), matches.end());

    return matches;

}

//Given a displacement file name, obtain the associated .vor and .cft file 
//names.
vector<string> get_vor_cft_names(string disp_name){

    ostringstream base_name, vor_oss, cft_oss;
    boost::char_separator<char> sep{"."};
    tokenizer tok{disp_name, sep};
    vector<string> tokens;
    vector<string> file_names;

    //Obtain the base name by tokenizing the displacement file name, using
    //periods as delimeters, and exclude everything from the final
    //period onward.
    base_name.str("");
    for(tokenizer::iterator it = tok.begin(); it != tok.end(); it++){
        tokens.push_back(*it);
    }

    for(int i = 0; i < tokens.size() - 1; i++){
        base_name << tokens[i] << ".";
    }

    //Create the Voronoi file name by appending the .vor file
    //extension.
    vor_oss << base_name.str() << "vor";
    file_names.push_back(vor_oss.str());

    //Create the CFT file by appending the .cft extension
    cft_oss << base_name.str() << "cft";
    file_names.push_back(cft_oss.str());

    return file_names;
}

//Given an initial set of points, and an offset, make copies of the original
//point set offset by the specified amount to both the left and right.
vector<Point> copy_right_left(vector<Point> original, double offset){

    vector<Point> output, left, right;
    
    for(Point p : original){
        left.push_back(Point(p.x - offset, p.y));
        right.push_back(Point(p.x + offset, p.y));
    }

    output.insert(output.end(), original.begin(), original.end());
    output.insert(output.end(), right.begin(), right.end());
    output.insert(output.end(), left.begin(), left.end());

    return output;
}

//Find the energy due to stretching and compression of bonds
void get_stretch_shear_energy(vector<Point> points, vector<edge_datum> edges, double offset, vector<vec2d> disps, double ks, double mu, vector<double> &energies){

    double dx, dy, length, dux, duy, de, estr, eshear;
    int iter, idx1, idx2;
    edge_datum curr_edge;

    for(iter = 0; iter < edges.size(); iter++){
        curr_edge = edges[iter];
        idx1 = get<0>(curr_edge);
        idx2 = get<1>(curr_edge);
        length = get<2>(curr_edge);

        dx = points[idx2].x + get<3>(curr_edge)*offset - points[idx1].x;
        dy = points[idx2].y - points[idx1].y;
        dux = get<0>(disps[idx2]) - get<0>(disps[idx1]);
        duy = get<1>(disps[idx2*2]) - get<1>(disps[idx1*2]);
        de = (dux*dx + duy*dy) / length;
        estr = .25 * ks * de*de;
        eshear = .25 * mu * (dux*dux + duy*duy - de*de);
	energies[idx1] += estr;
	energies[idx2] += estr;
	energies[idx1] += eshear;
	energies[idx2] += eshear;
    }
}


//Find the energy due to bending of bonds
void get_bending_energy(vector<Point> points, vector<triple_datum> tdata, double offset, vector<vec2d> disps, double kb, vector<double> &energies){

    double dx, dy, length_sq, dux, duy, du_c_r;
    int iter, idx1, idx2, idx3;
    triple_datum tdat;

    for(iter = 0; iter < tdata.size(); iter ++){
        tdat = tdata[iter];
        idx1 = get<0>(tdat);
        idx2 = get<1>(tdat);
        idx3 = get<2>(tdat);

        dx = points[idx2].x + offset*get<3>(tdat) - points[idx1].x;
        dy = points[idx2].y - points[idx1].y;
        length_sq = dx*dx + dy*dy;
        dux = get<0>(disps[idx2])+get<0>(disps[idx3])-2*get<0>(disps[idx1*2]);
        duy = get<1>(disps[idx2])+get<1>(disps[idx3])-2*get<1>(disps[idx1*2+1]);
        du_c_r = dux * dy - duy * dx;
        energies[idx1] += .5 * kb * du_c_r * du_c_r / length_sq;
    }
}

//Given a set of points, and sets of indices into those points indicating
//polygon vertices, find the area of each polygon.
void find_areas(vector<Point> pts, vector<vector<int>> polys, vector<double> &areas){

    double sum;
    size_t len;

    areas.clear();

    for(vector<int> next_poly : polys){
        sum = 0;
	len = next_poly.size();

	for(int next_index : next_poly){
            sum += pts[next_index].x * pts[(next_index + 1)%len].y;
            sum -= pts[(next_index + 1) % len].x * pts[next_index].y;
        }
    
	areas.push_back(.5 * abs(sum));
    }
}

int main(int argc, char **argv){

    FILE *vor_file, *cft_file;
    string net_pat, disp_pat;
    vector<string> net_names, disp_names, report_names;
    double ks, kb, mu, width, cutoff, fraction, minx, miny, maxx, maxy,frac;
    double maxrange;
    char ktype;
    ifstream input;
    vector<double> areas;

    if(argc < 2){
        cerr << "No job file was specified.\n";
	return -1;
    }

    //Attempt to read a batch file specifying one or more network description
    //files, one or more companion displacement files, mechanical attributes,
    //kernel information, and the fraction of eigenmodes to obtain.
    if(! parse_input_args(argv[1], net_pat, disp_pat, ks, kb, mu, ktype, width, cutoff, frac)){
        cerr << "An invalid input file was specified.\n";
        return -1;
    }

    //Upon successfully obtaining parameters, iterate over the list of network
    //and displacement files. Form a Voronoi tiling for each network, find the
    //energy for each tile, and obtain cellular Fourier transforms. Write
    //Voronoi and CFT data to output files. 
    net_names = get_file_list(net_pat);
    disp_names = get_file_list(disp_pat);

    for(int iter = 0; iter < net_names.size(); iter ++){
	//Obtain points, edges and displacements. If the bending stiffness
	//is non-zero, also obtain bending triples.
        vector<Point> points, ext_points, v_points, v_points_trimmed;
        vector<vector<int>> v_polygons, v_polygons_trimmed;
        vector<edge_datum> edges;
        vector<triple_datum> triples;
        vector<vec2d> disps;
        vector<int> left_set;
        set<int> keep_set;
	map<int, int> index_map;
	int count;
        
	input.open(net_names[iter]);
        read_network(input, points, edges, minx, maxx, miny, maxy, left_set);
        if(kb != 0){
            triples = get_triples(points, edges, maxx - minx);
        }

	input.open(disp_names[iter]);
        disps = read_displacements(input);

	//Make copies of points to the left and right of the original set
	ext_points = copy_right_left(points, maxx - minx);

	//Make a record of the indices of Voronoi cells that should be kept
	keep_set.clear();
	for(int k_iter = 0; k_iter < points.size(); k_iter++){
		keep_set.insert(k_iter);
	}
	count = points.size();
	for(int next_point_idx : left_set){
		//cout << "Next index: " << next_point_idx << "\n";
		keep_set.insert(next_point_idx + points.size());
		index_map.insert(make_pair(next_point_idx + points.size(), count ++));
	}
	
	//Obtain the raw, and subsequently the cropped, Voronoi tessellation
	maxrange = max(3*(maxx-minx), maxy - miny);
        get_voronoi_cells(ext_points, keep_set, v_points, v_polygons, 2*minx - maxx, maxrange, miny, maxrange);

        get_trimmed_polygons(v_points, v_polygons, v_points_trimmed, v_polygons_trimmed,  minx, miny, maxx, maxy);

	//Obtain areas and triangulations of all polygons
        find_areas(v_points_trimmed, v_polygons_trimmed, areas);

	//Find energy densities
	vector<double> energies(keep_set.size(), 0);
        get_stretch_shear_energy(points, edges, maxx-minx, disps, ks, mu, energies);
	if(kb > 0){
            get_bending_energy(points, triples, maxx-minx, disps, kb, energies);
        }

	//Duplicate energies along the left edge on its translated image along
	//the right edge of the network, then normalize energies by area.
	cout << "Number of cells: " << energies.size() << "\n";
	for(int next_idx : left_set){
            energies[index_map[next_idx + points.size()]] = energies[next_idx];
	}
	for(int eiter = 0; eiter < energies.size(); eiter++){
	    energies[eiter] /= areas[eiter];
	}
	
	//Carry out numerical integrals to find the elements of the discrete
	//Laplace operator, represented as a sparse matrix
	
	//Compute the desired number of eigenvalues of the discrete Laplace
	//operator, and infer the corresponding wave number.
	
	//Obtain report names, and write the Voronoi tessellation and CFT to 
	//files
	report_names = get_vor_cft_names(disp_names[iter]);
	vor_file = fopen(report_names[0].c_str(), "w");
	if(! vor_file){
	    cerr << "The Vorononoi file could not be opened for writing.\n";
	}
	else{
            //long int tally = 0;
	    //for(vector<int> poly : v_polygons) if(poly.size() >= 3) tally++;
            fprintf(vor_file, "%ld %ld\n", v_points_trimmed.size(), v_polygons_trimmed.size());

            for(Point p : v_points_trimmed){
                fprintf(vor_file, "%12.8lf\t%12.8lf\n", p.x, p.y);
            }

	    for(vector<int> next_poly : v_polygons_trimmed){
                for(int p_iter = 0; p_iter < next_poly.size() - 1; p_iter++){
                    fprintf(vor_file, "%d ", next_poly[p_iter]);
		}
		fprintf(vor_file, "%d\n", *next_poly.rbegin());
            }

            for(double energy : energies){
		fprintf(vor_file, "%12.8le\n", energy);
	    }

	    fclose(vor_file);
        }
	cft_file = fopen(report_names[1].c_str(), "w");
	for(Point p : points) fprintf(cft_file, "%lf\t%lf\n", p.x, p.y);
	fclose(cft_file);
    }

    return 0;
}
