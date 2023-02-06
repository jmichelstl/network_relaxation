/*
 * This program computes strains associated with a displacement field of a
 * strained, filamentous, elastic network.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <tuple>
#include <map>
#include <cfloat>
#include <cstring>
#include <stdio.h>
#include <vector>
#include <cmath>

using namespace std;

//Simple structure for representing a Cartesian point or vector in 2D
#define vec2d tuple<double, double>
#define edge tuple<int, int, short int>

//Read a set of vertices from a file
void read_points_edges(string file_name, double &offset, vector<vec2d> &points, vector<edge> &edges){

	double x, y, dummy;
	ifstream file_stream;
	string nextline;
	int num_points, num_edges, idx1, idx2;
	short int mult;

	file_stream.open(file_name);
	if(! file_stream.is_open()){
		cerr << "Error: Could not open " << file_name << ".\n";
		return;
	}

	getline(file_stream, nextline);
	if(sscanf(nextline.c_str(), "%d %d %lf\n", &num_points, &num_edges, &offset) < 3){
		cerr << "Invalid header detected.\n";
		file_stream.close();
		return;
	}

	while(! file_stream.eof() && points.size() < num_points){
		getline(file_stream, nextline);
		if(sscanf(nextline.c_str(),"%lf %lf",&x,&y) == 2){
			points.push_back(make_tuple(x, y));
		}
	}

	while(! file_stream.eof() && edges.size() < num_edges){
		getline(file_stream, nextline);
		if(sscanf(nextline.c_str(),"%d %d %hd",&idx1, &idx2, &mult)==3){
			edges.push_back(make_tuple(idx1, idx2, mult));
		}
	}

	file_stream.close();

	if(points.size() < num_points){
		cerr << "The indiciated number of points could not be read.\n";
	}

	if(edges.size() < num_edges){
		cerr << "The indicated number of edges could not be read.\n";
	}
}

//Read a displacement field from a file
void get_displacements(string file_name, vector<vec2d> &disps){

	double dx, dy;
	ifstream file_stream;
	string nextline;
	int num_read;

	file_stream.open(file_name);

	if(! file_stream.is_open()){
		cerr << file_name << " could not be opened.\n";
		return;
	}

	while(! file_stream.eof()){
		getline(file_stream, nextline);
		if((num_read = sscanf(nextline.c_str(),"%lf %lf",&dx, &dy))==2){
			disps.push_back(make_tuple(dx, dy));
		}
	}

	file_stream.close();

}

vector<double> get_strains(vector<vec2d> points, vector<edge> edges, double offset, vector<vec2d> disps){

	//Initial and final differences in x and y coordinates, and initial
	//and final lengths of edges
	double dxi, dyi, dxf, dyf, li, lf;
	int idx1, idx2, mult;
	vector<double> strains;

	for(edge next_edge : edges){
		idx1 = get<0>(next_edge);
		idx2 = get<1>(next_edge);
		mult = get<2>(next_edge);

		dxi = get<0>(points[idx2]) + mult*offset - get<0>(points[idx1]);
		dyi = get<1>(points[idx2]) - get<1>(points[idx1]);
		dxf = dxi + get<0>(disps[idx2]) - get<0>(disps[idx1]);
		dyf = dyi + get<1>(disps[idx2]) - get<1>(disps[idx1]);

		li = sqrt(dxi*dxi + dyi*dyi);
		lf = sqrt(dxf*dxf + dyf*dyf);
		strains.push_back((lf - li) / li);
	}

	return strains;
}

/*
 * Next, obtain the points and displacements for the reference case, and find
 * the minimum and maximum x and y values, for the purpose of computing the
 * area. Compute the non-affine parameter for each network as the mean square
 * difference between the displacement field at each vertex of a network and
 * the displacement field for the same point in the reference state. If a point
 * is missing from a given network, this vertex does not mechanically 
 * participate in the network, and its displacement field will be set to zero.
 */
int main(int argc, char **argv){

	vector<vec2d> points, disps;
	vector<edge> edges;
	vector<double> strains;
	double offset;
	ifstream datfile, dispfile;
	FILE *report;

	if(argc < 4){
		cerr << "Usage: dat file, disp file, strain file";
		return -1;
	}

	//Obtain points, edges and displacements
	read_points_edges(string(argv[1]), offset, points, edges);
	get_displacements(string(argv[2]), disps);

	if(points.size() != disps.size()){
		cerr << "Point and displacement count mismatch.\n";
		return -1;
	}
	
	strains = get_strains(points, edges, offset, disps);

	report = fopen(argv[3], "w");
	if(! report){
		cerr << "The file " << argv[3] << " could not be opened.\n";
		return -1;
	}

	for(double next_strain : strains){
		fprintf(report, "%2.10le\n", next_strain);
	}

	fclose(report);

	return 0;
}
