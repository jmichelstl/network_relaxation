#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <map>
#include <string>
#include <stdio.h>
#include <cmath>
#include <sstream>
#include <cfloat>
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>
#include <algorithm>
#include <set>
#include "file_lists.hpp"
#include <time.h>
#include <cuda_runtime.h>
#include "cublas_v2.h"
#include "float.h"
#include <iomanip>

#define FX_MASK 0x1
#define FY_MASK 0x2
#define FX_FY 0x3

using namespace std;

//Data structures for representing vectors, bonds in the network, and 
//interacting triples.
//typedef tuple<double, double> vec2d;
typedef tuple<int, int, double, char> edge_datum;
typedef tuple<int, int, int, char, char> triple_datum;

//Type for specifying the minimum value, maximum value, and increment for
//straining a network a network in a particular manner.
typedef tuple<double, double, double> strain_range;

#define FLT_TOL 1e-10

//Simple data structure for a Cartesian point
struct Point {

    Point(){}

    Point(double xval, double yval) : x(xval), y(yval) {} 

    double x, y;
};

//Struct to compare two Cartesian 2D points
//This struct sorts points first by x index, and then by y index, as the
//intended use case is to enumerate all points along the left edge of a
//network.
struct PointComp {
    bool operator()(const Point& lhs, const Point& rhs) const { 
        if(lhs.x < rhs.x - FLT_TOL) return true;
	else if(lhs.x > rhs.x + FLT_TOL) return false;
	else if(lhs.y < rhs.y - FLT_TOL) return true;
	else return false;
    }
};

bool in_bounds(int idx, int min, int max){
    return idx >= min && idx <= max;
}

//Read a network-mesh composite file, which specifies mesh points, network
//edges, mesh facets, and facet areas. Return a record of points, edges, facets,
//areas and pairs of points along the left and right edges that should be
//treated as duplicates of one another.
bool read_network_mesh(string file_name, double &offset, vector<double> &points, vector<edge_datum> &edges, vector<int> &facets, vector<double> &areas, vector<int> &mesh_pairs, double &min_x, double &min_y, double &max_x, double &max_y){

    ifstream input;
    string nextline;
    int num_read, num_points, num_edges, num_facets, idx1, idx2, idx3;
    double x, y, dx, dy, area;
    char mult;
    map<Point, int, PointComp> pmap;

    input.open(file_name);
    if(! input.is_open()) return false;

    min_x = FLT_MAX;
    min_y = FLT_MAX;
    max_x = FLT_MIN;
    max_y = FLT_MIN;

    //Attempt to read the header for the network description file. If too little
    //information or invalid data are found, close the input stream and report
    //a failure to read a network.
    getline(input, nextline);
    num_read = sscanf(nextline.c_str(), "%d %d %d %lf", &num_points, &num_edges,&num_facets, &offset);
    if(num_read < 4 || offset < 0 || num_points <= 0 || num_edges < 0 || num_facets < 0){
        cerr << "Reading failed. The provided file has an invalid header.\n";
        input.close();
        return false;
    }

    //If a valid header was read, attempt to read the specified number of points
    while(! input.eof() && points.size() < 2*num_points){
        getline(input, nextline);
        num_read = sscanf(nextline.c_str(), "%lf %lf", &x, &y);
        if(num_read == 2){
            pmap.insert(make_pair(Point(x, y), points.size() / 2));
            points.push_back(x);
            points.push_back(y);
            if(x < min_x) min_x = x;
            if(y < min_y) min_y = y;
            if(x > max_x) max_x = x;
            if(y > max_y) max_y = y;
        }
    }

    //If the program has reached the end of the input stream, too little
    //information has been provided, and the file is invalid.
    if(input.eof()){
        cerr << "Reading failed. Too few points were read.\n";
        input.close();
        return false;
    }

    //Attempt to read edges
    while(!input.eof() && edges.size() < num_edges){
        getline(input, nextline);
        num_read = sscanf(nextline.c_str(),"%d %d %hhd\n", &idx1, &idx2, &mult);

        //Ensure two point indices and an offset are given. If so, and the point
        //indices are in bounds, find the length of the edge, and add it to
        //the list.
        if(num_read == 3){
            if(min(idx1, idx2) > -1 && max(idx1, idx2) < num_points){
                dx = points[2*idx2] + offset*mult -points[2*idx1];
                dy = points[2*idx2 + 1] - points[2*idx1 + 1];
                edges.push_back(make_tuple(idx1, idx2,sqrt(dx*dx+dy*dy), mult));
            }
            else{
                cerr << "Reading failed due to an out-of-bounds point index.\n";
                input.close();
                return false;
            }
        }
    }

    //If the number of edges is less than expected, an error has occurred,
    //and parsing should be aborted.
    if(edges.size() < num_edges){
        cerr << "Too few edges were read.\n";
	input.close();
	return false;
    }

    //If a non-zero number of facets is expected, but the end of the file has
    //been reached, an error has occurred.
    if(input.eof() && num_facets > 0){
        cerr << "The file ended before facets could be read.\n";
	input.close();
	return false;
    }

    //Attempt to read the specified number of facets
    while(!input.eof() && facets.size() < 3 * num_facets){
        getline(input, nextline);
        num_read = sscanf(nextline.c_str(),"%d %d %d\n", &idx1, &idx2, &idx3);

        //Ensure three facet indices are given. If so, and the point indices are
       	//in bounds, add a record of the facet to the running list.
        if(num_read == 3){
            if(in_bounds(idx1,0,num_points-1) && in_bounds(idx2,0,num_points-1) && in_bounds(idx3, 0, num_points-1)){
                dx = points[2*idx2] + offset*mult -points[2*idx1];
                dy = points[2*idx2 + 1] - points[2*idx1 + 1];
		facets.push_back(idx1);
		facets.push_back(idx2);
		facets.push_back(idx3);
            }
            else{
                cerr << "Reading failed due to an out-of-bounds point index.\n";
                input.close();
                return false;
            }
        }
    }

    //If the number of facets is less than expected, an error has occurred,
    //and parsing should be aborted.
    if(facets.size() < 3 * num_facets){
        cerr << "Too few facets were read.\n";
	input.close();
	return false;
    }

    //If the number of facets is non-zero, but no facet areas are specified,
    //report an error.
    if(input.eof() && num_facets > 0){
        cerr << "The file ended before facets could be read.\n";
	input.close();
	return false;
    }

    //Attempt to read areas
    while(! input.eof() && areas.size() < num_facets){
        getline(input, nextline);
        if(sscanf(nextline.c_str(), "%lf", &area) == 1){
            if(area > 0) areas.push_back(area);
	    else{
                cerr << "An illegal area was specified.\n";
		input.close();
		return false;
            }
	}
    }

    input.close();

    if(areas.size() < num_facets){
        cerr << "Too few areas were specified.\n";
	return false;
    }

    //Find pairs of duplicate points along the left and right edges that should
    //be subject to the same forces.
    for(auto iter = pmap.begin(); iter != pmap.end() && iter->first.x < min_x + FLT_TOL; iter ++){
        if(pmap.find(Point(iter->first.x+offset, iter->first.y)) != pmap.end()){
            mesh_pairs.push_back(iter->second);
	    mesh_pairs.push_back(pmap[Point(iter->first.x+offset,iter->first.y)]);
        }
    }

    return true;
}

vector<triple_datum> get_triples(vector<double> points, vector<edge_datum> edges, double offset){

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

        adjacency_lists[idx1].push_back(make_tuple(idx2, get<3>(next_edge)));
        adjacency_lists[idx2].push_back(make_tuple(idx1, -get<3>(next_edge)));
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
        x = points[2*idx1];
        y = points[2*idx1 + 1];

        for(i = 0; i < neighbors.size(); i ++){
            pair1 = neighbors[i];
            idx2 = get<0>(pair1);
            dx2 = points[2*idx2] + offset*get<1>(pair1) - x;
            dy2 = points[2*idx2+1] - y;

            for(j = i+1; j < neighbors.size(); j++){
                pair2 = neighbors[j];
                idx3 = get<0>(pair2);

                dx3 = points[2*idx3] + offset*get<1>(pair2) - x;
                dy3 = points[2*idx3+1] - y;

                dot = dx2*dx3 + dy2*dy3;
                mag_sq = (dx2*dx2 + dy2*dy2)*(dx3*dx3 + dy3*dy3);

                //If the dot product of the vectors pointing from the hinge
                //to either end point is negative, and the square of this
                //dot product is equal to the product of the square of the
                //magnitudes of the two displacement vectors, the two 
                //displacement vectors are antiparallel, and a valid triple has
                //been identified.
                if(dot < 0 && abs(dot*dot  - mag_sq) < FLT_TOL){
                    triples.push_back(make_tuple(idx1,idx2,idx3,get<1>(pair1),get<1>(pair2)));
                }
            }
        }
    }

    return triples;
}

//Utility function to find the number of digits in the base-ten representation
//of an integer, for the purpose of appending a zero-padded integer to the
//base name for a report file
int get_dec_digits(int value){
    int num_digits = 0;

    do{
        num_digits++;
        value /= 10;
    }while(value > 0);

    return num_digits;
}

//Utility function to create a file name, given a base, the number of decimal
//digits to use in an appended count, and the position of the report file in
//a series.
string report_name(string base, int num_digits, int count, string ext){
    int padding, iter;
    ostringstream oss;

    padding = num_digits - get_dec_digits(count);
    oss << base << "_";
    for(iter = 0; iter < padding; iter++) oss << 0;
    oss << count << "." << ext;

    return oss.str();
}

//Report displaced points to a file
void report_displaced(string file_name, vector<double> points){

    FILE *report_file;
    int iter;

    report_file = fopen(file_name.c_str(), "w");
    for(iter = 0; iter < points.size() / 2; iter++){
        fprintf(report_file, "%2.10le\t%2.10le\n", points[iter*2], points[iter*2 + 1]);
    }

    fclose(report_file);
}

//Split a line into pieces, using whitespace characters as delimiters
vector<string> split_line(string line){
    vector<string> result;
    char *full_line, *token;

    full_line = (char *) malloc((1 + line.size()) * sizeof(char));
    strcpy(full_line, line.c_str());

    token = strtok(full_line, " \t");
    while(token != NULL){
        result.push_back(string(token));
        token = strtok(NULL, " \t");
    }

    free(full_line);

    return result;
}

vector<double> parse_doubles(vector<string> numstring){
    double nextnum;
    vector<double> result;

    for(string s : numstring){
        if(sscanf(s.c_str(), "%lf", &nextnum)) result.push_back(nextnum);
    }

    return result;
}

//Attempt to read a line from an input stream. First, make sure there is
//something to read. Next, read until either a non-comment and non-blank line
//is encountered, or the end of the file is reached. Return true upon success,
//and false upon failure.
bool read_line(ifstream &input, string &nextline){

    if(input.eof()){
        cerr << "Not enough information was provided in the batch file.\n";
        input.close();
        return false;
    }

    while(! input.eof()){
        getline(input, nextline);
        if(nextline.compare("") == 0 || nextline[0] == '#') continue;
        else break;
    }

    if(nextline.compare("") == 0 || nextline[0] == '#'){
        input.close();
        return false;
    }

    return true;
}

//Parse a file containing instructions for a batch of simulations. The input
//file should countain the following:
//    -base name for network description (.dat) files
//    -mechanical attributes of network fibers and the background gel
//    -the range of compressions to apply
//    -the range of shears to apply
//Information should be stated in the order given above. Blank lines, and lines
//beginning with a pound sign, will be ignored.
//
//Return value:
//    -true, if all necessary information could be obtained
//    -false, otherwise
bool parse_batch_file(string file_name, string &base_name, string &tag, double *ks, double *kb, double *bmod, double *smod, double *dilation, double fp_list[5], int *nmin, double *mass, double *fcut, int *nmax, strain_range &compr, strain_range &shearr){

    string nextline;
    ifstream input;
    vector<string> tokens;
    int num_read, iter;
    double min, max, inc;
    vector<double> fp_vector;
    bool all_found;

    //Make sure the batch file can be opened for reading
    input.open(file_name);
    if(! input.is_open()){
        cerr << "The file " << file_name << " could not be read.\n";    
        return false;
    }

    //Read the base name and an optional tag
    if(! read_line(input, nextline)) return false;

    //Initialize the base name, and, optionally, a tag to append to the base
    //name.
    tokens = split_line(nextline);
    base_name = tokens[0];
    if(tokens.size() > 1) tag = tokens[1];

    //Obtain the mechanical attributes of the network and gel
    if(! read_line(input, nextline)) return false;
        
    num_read = sscanf(nextline.c_str(), "%lf %lf %lf %lf %lf", ks, kb, bmod, smod, dilation);
    if(num_read < 3 || *ks < 0 || *kb < 0 || *bmod < 0 || *smod < 0 || *dilation < 0){
        cerr << "Error: invalid format for mechanical attributes.";
        input.close();
        return false;
    }

    //Obtain the parameters for structural relaxation according to the FIRE
    //method
    if(! read_line(input, nextline)) return false;
    
    fp_vector = parse_doubles(split_line(nextline));
    if(fp_vector.size() < 5){
        cerr << "Too few FIRE parameters were supplied. 8 are needed.\n";
        input.close();
        return false;
    }
    else{
        for(iter = 0; iter < 5; iter++) fp_list[iter] = fp_vector[iter];
    }

    if(! read_line(input, nextline)) return false;
    num_read = sscanf(nextline.c_str(),"%d %lf %lf %d", nmin, mass, fcut, nmax);
    if(num_read < 4 || mass < 0 || fcut < 0 || nmax < 0){
        cerr<<"On or more of nmin, mass, cutoff force and step limit were invalid.\n";
        input.close();
        return false;
    }

    //Obtain the parameters for compressing the network
    if(! read_line(input, nextline)) return false;

    num_read = sscanf(nextline.c_str(), "%lf %lf %lf", &min, &max, &inc);
    if(num_read < 3 || min < 0 || max < 0 || inc < 0 || max < min){
           cerr << "Error: invalid specification of compression.";
          input.close();
           return false;
    }
    else compr = make_tuple(min, max, inc);

    //Obtain the parameters for shearing the network
    if(! read_line(input, nextline)) return false;

    num_read = sscanf(nextline.c_str(), "%lf %lf %lf", &min, &max, &inc);
    if(num_read < 3 || min < 0 || max < 0 || inc < 0 || max < min){
        cerr << "Error: invalid specification of shear.";
           input.close();
        return false;
    }

    else{
        shearr = make_tuple(min, max, inc);
           all_found = true;
    }

    input.close();
    return all_found;
}

//Pre-compute inverse matrices used in finding affine transformations from
//reference to the deformed state
__global__ void get_inv_mats(int n_facets, double *pos, int *facets, double *inv_mats){

    int fIdx, idx1, idx2, idx3;
    double x1, y1, x2, y2, x3, y3, det;

    fIdx = blockIdx.x * blockDim.x * blockDim.y;
    fIdx += threadIdx.y * blockDim.x + threadIdx.x;
    if(fIdx >= n_facets) return;

    idx1 = facets[3*fIdx];
    idx2 = facets[3*fIdx + 1];
    idx3 = facets[3*fIdx + 2];

    x1 = pos[2*idx1];
    y1 = pos[2*idx1 + 1];
    x2 = pos[2*idx2];
    y2 = pos[2*idx2 + 1];
    x3 = pos[2*idx3];
    y3 = pos[2*idx3 + 1];

    det = x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2);

    inv_mats[6*fIdx] = (y2 - y3) / det; 
    inv_mats[6*fIdx + 1] = (x3 - x2) / det; 
    inv_mats[6*fIdx + 2] = (y3 - y1) / det; 
    inv_mats[6*fIdx + 3] = (x1 - x3) / det; 
    inv_mats[6*fIdx + 4] = (y1 - y2) / det; 
    inv_mats[6*fIdx + 5] = (x2 - x1) / det;
}

//Find forces due to stretching and compression of bonds
__global__ void axial_forces(double ks, int npairs, double width, double *pos, unsigned char *fcodes, int *pair_list, double *lengths, char *wrap_codes, double *forces){
    int pair_index, index1, index2;
    double x1, x2, y1, y2, inv_dist, fx, fy;

    pair_index = blockIdx.x * blockDim.x * blockDim.y;
    pair_index += threadIdx.y * blockDim.x + threadIdx.x;
    if(pair_index >= npairs) return;

    index1 = pair_list[pair_index * 2];
    index2 = pair_list[pair_index * 2 + 1];

    x1 = pos[index1*2];
    y1 = pos[index1*2 + 1];
    x2 = pos[2*index2];
    y2 = pos[2*index2 + 1];
    inv_dist = rhypot(x1 - x2 - wrap_codes[pair_index]*width, y1 - y2);
    fx = ks*(1-lengths[pair_index]*inv_dist)*(x2 + wrap_codes[pair_index]*width-x1);
    fy = ks*(1-lengths[pair_index]*inv_dist)*(y2 - y1);

    //If a node is not at the bottom or top of the network, relax
    //the x and y coordinates. If a node is on the top or bottom, determine
    //whether the node may relax in the x direction, and if so, relax in the x
    //direction.
    atomicAdd(forces + index1*2, fx*(fcodes[index1] & FX_MASK));
    atomicAdd(forces + index1*2 + 1, fy*((fcodes[index1] & FY_MASK) / 2));

    atomicAdd(forces + index2*2, -fx*(fcodes[index2] & FX_MASK));
    atomicAdd(forces + index2*2 + 1, -fy*((fcodes[index2] & FY_MASK) / 2));
}

//Find forces due to bending of bonds
__global__ void bending_forces(double kb, int ntriples, double width, double *pos, unsigned char *fcodes, int *triple_list, char *wrap_codes, double *forces){

    //Positions, displacement vector components for computing bending forces,
    //and geometric quantities associated with those components
    double x1, y1, x2, y2, x3, y3, rijx, rijy, rikx, riky, invLen1, invLen2;
    double ang, f2x, f2y, f3x, f3y;
    int triple_idx, idx1, idx2, idx3;

    triple_idx = blockIdx.x * blockDim.x * blockDim.y;
    triple_idx += threadIdx.y * blockDim.x + threadIdx.x;
    if(triple_idx >= ntriples) return;

    idx1 = triple_list[triple_idx * 3];
    idx2 = triple_list[triple_idx * 3 + 1];
    idx3 = triple_list[triple_idx * 3 + 2];

    x1 = pos[idx1*2];
    y1 = pos[idx1*2 + 1];
    x2 = pos[idx2*2];
    y2 = pos[idx2*2 + 1];
    x3 = pos[idx3*2];
    y3 = pos[idx3*2 + 1];

    //printf("Triple Index: %d, Index 1: %d, Index 2: %d, Index 3: %d, x1: %lf, y1: %lf, x2: %lf, y2: %lf, x3: %lf, y3: %lf\n", triple_idx, idx1, idx2, idx3, x1, y1, x2, y2, x3, y3);

    rijx = x2 + wrap_codes[2*triple_idx]*width - x1;
    rijy = y2 - y1;
    rikx = x3 + wrap_codes[2*triple_idx + 1]*width - x1;
    riky = y3 - y1;

    invLen1 = rhypot(rijx, rijy);
    invLen2 = rhypot(rikx, riky);

    ang = atan2(rijx*riky - rijy*rikx, rijx*rikx + rijy*riky);
    if(ang < 0) ang += 2*M_PI;

    f2x = kb * (ang - M_PI) * invLen1 * invLen1;
    f2y = f2x * rijx;
    f2x *= -rijy;

    f3x = kb * (ang - M_PI) * invLen2 * invLen2;
    f3y = -f3x * rikx;
    f3x *= riky;

    atomicAdd(forces + idx1*2, (-f2x - f3x)*(fcodes[idx1] & FX_MASK));
    atomicAdd(forces + idx1*2 + 1, (-f2y - f3y)*((fcodes[idx1] & FY_MASK) / FY_MASK));

    atomicAdd(forces + idx2*2, f2x*(fcodes[idx2] & FX_MASK));
    atomicAdd(forces + idx2*2 + 1, f2y*((fcodes[idx2] & FY_MASK)/FY_MASK));

    atomicAdd(forces + idx3*2, f3x*(fcodes[idx3] & FX_MASK));
    atomicAdd(forces + idx3*2 + 1, f3y*((fcodes[idx3] & FY_MASK)/FY_MASK));
}

//Find energies due to stretching and compression of bonds
__global__ void axial_energies(double ks, int npairs, double width, double *pos, int *pair_list, double *lengths, char *wrap_codes, double *d_spe){
    int pair_idx, index1, index2;
    double x1, x2, y1, y2, dist;

    pair_idx = blockIdx.x * blockDim.x * blockDim.y;
    pair_idx += threadIdx.y * blockDim.x + threadIdx.x;
    if(pair_idx >= npairs) return;

    index1 = pair_list[pair_idx * 2];
    index2 = pair_list[pair_idx * 2 + 1];

    x1 = pos[index1*2];
    y1 = pos[index1*2 + 1];
    x2 = pos[2*index2];
    y2 = pos[2*index2 + 1];
    dist = hypot(x1 - x2 - wrap_codes[pair_idx]*width, y1 - y2);

    d_spe[pair_idx] = .5*ks*(dist-lengths[pair_idx])*(dist-lengths[pair_idx]);
}

//Find the strain energy contribution due to bending of bonds
__global__ void bending_energies(double kb, int ntriples, double width, double *pos, int *triple_list, char *wrap_codes, double *d_bpe){

    //Positions, displacement vector components for computing bending forces,
    //and geometric quantities associated with those components
    double x1, y1, x2, y2, x3, y3, rijx, rijy, rikx, riky, ang;
    int triple_idx, idx1, idx2, idx3;

    triple_idx = blockIdx.x * blockDim.x * blockDim.y;
    triple_idx += threadIdx.y * blockDim.x + threadIdx.x;
    if(triple_idx >= ntriples) return;

    idx1 = triple_list[triple_idx * 3];
    idx2 = triple_list[triple_idx * 3 + 1];
    idx3 = triple_list[triple_idx * 3 + 2];

    x1 = pos[idx1*2];
    y1 = pos[idx1*2 + 1];
    x2 = pos[idx2*2];
    y2 = pos[idx2*2 + 1];
    x3 = pos[idx3*2];
    y3 = pos[idx3*2 + 1];

    rijx = x2 + wrap_codes[2*triple_idx]*width - x1;
    rijy = y2 - y1;
    rikx = x3 + wrap_codes[2*triple_idx + 1]*width - x1;
    riky = y3 - y1;

    ang = atan2(rijx*riky - rijy*rikx, rijx*rikx + rijy*riky);
    if(ang < 0) ang += 2*M_PI;
    d_bpe[triple_idx] = .5 * kb * (ang - M_PI)*(ang - M_PI);
}

//Find energies due to deformation of FEM elements
__global__ void gel_energy(double smod, double bmod, double dilation, int nelems, double *pos, int *elems, double *inv_mats, double *areas, double *energies){

    int fIdx, idx1, idx2, idx3;
    double f11, f12, f21, f22, det;
    
    fIdx = blockIdx.x * blockDim.x * blockDim.y;
    fIdx += threadIdx.y * blockDim.x + threadIdx.x;
    if(fIdx >= nelems) return;

    //Find the elements of the deformation gradient tensor
    idx1 = elems[fIdx*3];
    idx2 = elems[fIdx*3 + 1];
    idx3 = elems[fIdx*3 + 2];


    f11 = pos[2*idx1]*inv_mats[fIdx*6] + pos[2*idx2]*inv_mats[fIdx*6+2] + pos[2*idx3]*inv_mats[fIdx*6 + 4];
    f12 = pos[2*idx1]*inv_mats[fIdx*6+1] + pos[2*idx2]*inv_mats[fIdx*6+3] + pos[2*idx3]*inv_mats[fIdx*6 + 5];
    f21 = pos[2*idx1+1]*inv_mats[fIdx*6] + pos[2*idx2+1]*inv_mats[fIdx*6+2] + pos[2*idx3+1]*inv_mats[fIdx*6 + 4];
    f22 = pos[2*idx1+1]*inv_mats[fIdx*6+1] + pos[2*idx2+1]*inv_mats[fIdx*6+3] + pos[2*idx3+1]*inv_mats[fIdx*6 + 5];

    //Find the determinant, and compute the energy
    det = f11*f22 - f12*f21;

    energies[fIdx] = smod*areas[fIdx]*((f11*f11+f12*f12+f21*f21+f22*f22)/det - 2)/2;
    energies[fIdx] += bmod*areas[fIdx]*(det-1-dilation)*(det-1-dilation)/2;
}

//Find the forces on the vertices of a finite element
__global__ void gel_forces(double smod, double bmod, double dilation, int nelems, double *pos, int *elems, double *inv_mats, double *areas, unsigned char *fcodes, double *forces){

    int fIdx, idx1, idx2, idx3;
    double f11, f12, f21, f22, det, coeff1, coeff2;
    double du_df11, du_df12, du_df21, du_df22;
    
    fIdx = blockIdx.x * blockDim.x * blockDim.y;
    fIdx += threadIdx.y * blockDim.x + threadIdx.x;
    if(fIdx >= nelems) return;

    //Find the elements of the deformation gradient tensor
    idx1 = elems[fIdx*3];
    idx2 = elems[fIdx*3 + 1];
    idx3 = elems[fIdx*3 + 2];


    f11 = pos[2*idx1]*inv_mats[fIdx*6] + pos[2*idx2]*inv_mats[fIdx*6+2] + pos[2*idx3]*inv_mats[fIdx*6 + 4];
    f12 = pos[2*idx1]*inv_mats[fIdx*6+1] + pos[2*idx2]*inv_mats[fIdx*6+3] + pos[2*idx3]*inv_mats[fIdx*6 + 5];
    f21 = pos[2*idx1+1]*inv_mats[fIdx*6] + pos[2*idx2+1]*inv_mats[fIdx*6+2] + pos[2*idx3+1]*inv_mats[fIdx*6 + 4];
    f22 = pos[2*idx1+1]*inv_mats[fIdx*6+1] + pos[2*idx2+1]*inv_mats[fIdx*6+3] + pos[2*idx3+1]*inv_mats[fIdx*6 + 5];

    det = f11*f22 - f12*f21;
    coeff1 = -smod*areas[fIdx] / det;
    coeff2 = bmod * (dilation + 1 - det) + smod*(f11*f11+f12*f12+f21*f21+f22*f22)/(2*det*det);
    coeff2 *= areas[fIdx];

    du_df11 = f11*coeff1 + f22*coeff2;
    du_df12 = f12*coeff1 - f21*coeff2;
    du_df21 = f21*coeff1 - f12*coeff2;
    du_df22 = f22*coeff1 + f11*coeff2;

    //Forces on vertex 1
    atomicAdd(forces+idx1*2, du_df11*inv_mats[fIdx*6]*(fcodes[idx1] & FX_MASK));
    atomicAdd(forces+idx1*2, du_df12*inv_mats[fIdx*6+1]*(fcodes[idx1]&FX_MASK));
    atomicAdd(forces+idx1*2+1, du_df21*inv_mats[fIdx*6]*((fcodes[idx1]&FY_MASK)/FY_MASK));
    atomicAdd(forces+idx1*2+1, du_df22*inv_mats[fIdx*6+1]*((fcodes[idx1]&FY_MASK)/FY_MASK));

    //Forces on vertex 2
    atomicAdd(forces+idx2*2, du_df11*inv_mats[fIdx*6+2]*(fcodes[idx2]&FX_MASK));
    atomicAdd(forces+idx2*2, du_df12*inv_mats[fIdx*6+3]*(fcodes[idx2]&FX_MASK));
    atomicAdd(forces+idx2*2+1, du_df21*inv_mats[fIdx*6+2]*((fcodes[idx2]&FY_MASK)/FY_MASK));
    atomicAdd(forces+idx2*2+1, du_df22*inv_mats[fIdx*6+3]*((fcodes[idx2]&FY_MASK)/FY_MASK));

    //Forces on vertex 3
    atomicAdd(forces+idx3*2, du_df11*inv_mats[fIdx*6+4]*(fcodes[idx3]&FX_MASK));
    atomicAdd(forces+idx3*2, du_df12*inv_mats[fIdx*6+5]*(fcodes[idx3]&FX_MASK));
    atomicAdd(forces+idx3*2+1, du_df21*inv_mats[fIdx*6+4]*((fcodes[idx3]&FY_MASK)/FY_MASK));
    atomicAdd(forces+idx3*2+1, du_df22*inv_mats[fIdx*6+5]*((fcodes[idx3]&FY_MASK)/FY_MASK));
}

//Given a set of vertex pairs, set the force on each vertex to be the sum
//of the forces on that vertex and its duplicate
__global__ void duplicate_forces(unsigned char *fcodes, double *forces, int *m_pairs, int npairs){

    int pIdx, idx1, idx2;
    double temp;

    pIdx = blockIdx.x * blockDim.x * blockDim.y;
    pIdx += threadIdx.y * blockDim.x + threadIdx.x;
    if(pIdx >= npairs) return;

    idx1 = m_pairs[2*pIdx];
    idx2 = m_pairs[2*pIdx + 1];

    //Copy x components
    temp = forces[2 * idx1];
    //atomicAdd(forces + idx1*2, forces[2*idx2]*(fcodes[idx1] & FX_MASK));
    atomicAdd(forces + idx1*2, forces[2*idx2]);
    //atomicAdd(forces + idx2*2, temp*(fcodes[idx2] & FX_MASK));
    atomicAdd(forces + idx2*2, temp);

    //Copy y components
    temp = forces[2 * idx1 + 1];
    //atomicAdd(forces+idx1*2+1, forces[2*idx2+1]*((fcodes[idx1]&FY_MASK)/FY_MASK));
    atomicAdd(forces+idx1*2+1, forces[2*idx2+1]);
    //atomicAdd(forces+idx2*2+1, temp*((fcodes[idx2]&FY_MASK)/FY_MASK));
    atomicAdd(forces+idx2*2+1, temp);
}

//Given a set of vertices, and minimum and maximum y values, apply an affine
//displacement to each vertex.
__global__ void affine_displace(int npoints, double *points, double miny, double yrange, double xdisp, double ydisp){

    int pIdx;
    
    pIdx = blockIdx.x * blockDim.x * blockDim.y;
    pIdx += threadIdx.y * blockDim.x + threadIdx.x;
    if(pIdx >= npoints) return;

    points[2 * pIdx] += xdisp * (points[2 * pIdx + 1] - miny) / yrange;
    points[2 * pIdx + 1] += ydisp * (points[2 * pIdx + 1] - miny) / yrange;
}

//Uniformly expand a network
__global__ void dilate(int npoints, double *points, double minx, double miny, double dilation){

    int pIdx;

    pIdx = blockIdx.x * blockDim.x * blockDim.y;
    pIdx += threadIdx.y * blockDim.x + threadIdx.x;
    if(pIdx >= npoints) return;

    points[2 * pIdx] = minx + dilation * (points[2 * pIdx] - minx);
    points[2 * pIdx + 1] = miny + dilation * (points[2 * pIdx + 1] - miny);
}

//Ensure a mission-critical CUDA operation is successful. If it fails,
//print an error message, and indicate the program cannot proceed.
bool check(cudaError_t error){
    if(error != cudaSuccess){
        printf("code: %d, reason: %s\n", error, cudaGetErrorString(error));
        return false;
    }

    return true;
}

//Copy data from CPU-side representations of bonds to GPU-side representations
//Return true if all copies were succesfully made, and false otherwise.
void copy_bond_data(vector<edge_datum> edges, int *d_pairs, double *d_lengths, char *d_offsets){

    //Structures to unpack data about edges in a manner amenable to 
    //GPU-accelerated computation
    vector<int> h_pairs;
    vector<double> h_lengths;
    vector<char> h_offsets;

    //Copy data from the host to the device
    for(edge_datum edat : edges){
        h_pairs.push_back(get<0>(edat));
        h_pairs.push_back(get<1>(edat));
        h_lengths.push_back(get<2>(edat));
        h_offsets.push_back(get<3>(edat));
    }

    cudaMemcpy(d_pairs, &h_pairs[0], 2*sizeof(int)*edges.size(), cudaMemcpyHostToDevice);
    cudaMemcpy(d_lengths, &h_lengths[0], sizeof(double)*edges.size(),cudaMemcpyHostToDevice);
    cudaMemcpy(d_offsets, &h_offsets[0], sizeof(char)*edges.size(),cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();
}

//Copy data from GPU-side representations of bending triples to GPU-side
//representations
void copy_triple_data(vector<triple_datum> triples, int *d_triples, char *d_offsets){

    //Structures to unpack data about triples into a GPU-friendly form
    vector<int> h_indices;
    vector<char> h_offsets;

    for(triple_datum tdat : triples){
        h_indices.push_back(get<0>(tdat));
        h_indices.push_back(get<1>(tdat));
        h_indices.push_back(get<2>(tdat));
        h_offsets.push_back(get<3>(tdat));
        h_offsets.push_back(get<4>(tdat));
    }

    //Copy data from the CPU to the GPU
    cudaMemcpy(d_triples, &h_indices[0], 3*sizeof(int)*triples.size(), cudaMemcpyHostToDevice);
    cudaMemcpy(d_offsets, &h_offsets[0], 2*sizeof(char)*triples.size(), cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();
}

bool run_fire_md(int num_dof, int npairs, int ntriples, int nelems, double width, double ks, double kb, double smod, double bmod, double dilation, double *d_pos, double *d_vel, double *d_forces, int *d_pair_list, double *d_lengths, char *d_e_offsets, int *d_triples, char *d_t_offsets, int *d_elems, double *d_areas, double *d_inv_mats, int *d_mpairs, int n_mpairs, unsigned char* d_fcodes, double fire_params[5], double mass, double fcut, int max_steps, int nmin, dim3 e_grid, dim3 block, dim3 t_grid, dim3 m_grid, dim3 mp_grid, double report_data[4], int &step_count){

    double alpha, finc, fdec, alpha_start, falpha, fire_dt_max, sqrt_dof;
    double power, dt, dt_mult, dt_sq_mult, vmag, fmag, a_vm_o_fm, neg_alpha;
    double *d_spe, *d_bpe, *d_gpe;
    int since_leq_0 = 0;
    cublasHandle_t handle;

    step_count = 0;

    //Initialize handle for CUDA-accelerated blas calculations
    cublasCreate(&handle);

    //Unpack parameters of FIRE minimization scheme and initialize values
    alpha_start = fire_params[0];
    falpha = fire_params[1];
    fire_dt_max = fire_params[2];
    finc = fire_params[3];
    fdec = fire_params[4];

    alpha = alpha_start;
    neg_alpha = -alpha;
    dt = fire_dt_max;
    dt_sq_mult = .5 * dt*dt / mass;
    dt_mult = .5 * dt / mass;
    sqrt_dof = sqrt(num_dof);

    cudaMemset(d_vel, 0, sizeof(double) * num_dof);

    //Find the forces at the outset
    cudaMemset(d_forces, 0, sizeof(double) * num_dof);
    cudaDeviceSynchronize();

    if(ks > 0){
        axial_forces<<<e_grid, block>>>(ks, npairs, width, d_pos, d_fcodes, d_pair_list, d_lengths, d_e_offsets, d_forces);
        cudaDeviceSynchronize();
    }

    if(kb > 0){
        bending_forces<<<t_grid, block>>>(kb, ntriples, width, d_pos, d_fcodes, d_triples, d_t_offsets, d_forces);
        cudaDeviceSynchronize();
    }

    if(bmod > 0){
        gel_forces<<<m_grid, block>>>(smod, bmod, dilation, nelems, d_pos, d_elems, d_inv_mats, d_areas, d_fcodes, d_forces);
        cudaDeviceSynchronize();
        duplicate_forces<<<mp_grid, block>>>(d_fcodes, d_forces, d_mpairs, n_mpairs);
        cudaDeviceSynchronize();
    }

    //cublasDnrm2(handle, num_dof, d_vel, 1, &vmag);
    vmag = 0;

    cublasDnrm2(handle, num_dof, d_forces, 1, &fmag);
    cudaDeviceSynchronize();

    //cout << "Starting RMS Force: " << fmag / sqrt_dof << "\n";
 
    //Perform molecular dynamics steps using velocity verlet method until the
    //kinetic energy cutoff is reached, or the maximum number of steps have
    //taken place
    while(step_count < max_steps){

        step_count ++;        

        //Update positions and velocities
        cublasDaxpy(handle, num_dof, &dt, d_vel, 1, d_pos, 1);
        cudaDeviceSynchronize();
        cublasDaxpy(handle, num_dof, &dt_sq_mult, d_forces, 1, d_pos, 1);
        cudaDeviceSynchronize();
        cublasDaxpy(handle, num_dof, &dt_mult, d_forces, 1, d_vel, 1);
        cudaDeviceSynchronize();

        //Calculate forces
        cudaMemset(d_forces, 0, sizeof(double) * num_dof);
        cudaDeviceSynchronize();

	if(ks > 0){
            axial_forces<<<e_grid, block>>>(ks, npairs, width, d_pos, d_fcodes, d_pair_list, d_lengths, d_e_offsets, d_forces);
            cudaDeviceSynchronize();
	}

        if(kb > 0){
            bending_forces<<<t_grid, block>>>(kb, ntriples, width, d_pos, d_fcodes, d_triples, d_t_offsets, d_forces);
            cudaDeviceSynchronize();
        }

	if(bmod > 0){
            gel_forces<<<m_grid, block>>>(smod, bmod, dilation, nelems, d_pos, d_elems, d_inv_mats, d_areas, d_fcodes, d_forces);
            cudaDeviceSynchronize();
            duplicate_forces<<<mp_grid, block>>>(d_fcodes, d_forces, d_mpairs, n_mpairs);
            cudaDeviceSynchronize();
        }

        //Update velocities and calculate power
        cublasDaxpy(handle, num_dof, &dt_mult, d_forces, 1, d_vel, 1);
        cudaDeviceSynchronize();
        cublasDdot(handle, num_dof, d_vel, 1, d_forces, 1, &power);
        cudaDeviceSynchronize();

        //Adjust velocities according to FIRE algorithm
        cublasDnrm2(handle, num_dof, d_vel, 1, &vmag);
        cublasDnrm2(handle, num_dof, d_forces, 1, &fmag);
        cudaDeviceSynchronize();
        a_vm_o_fm = alpha * vmag / fmag;

        cublasDaxpy(handle, num_dof, &neg_alpha, d_vel, 1, d_vel, 1);
        cudaDeviceSynchronize();
        cublasDaxpy(handle, num_dof, &a_vm_o_fm, d_forces, 1, d_vel, 1);
        cudaDeviceSynchronize();
        
        //Adjust FIRE parameters according to current power
        if(power > 0){
            since_leq_0 ++;
            if(since_leq_0 > nmin){
                dt = min(dt*finc, fire_dt_max);
                dt_mult = .5 * dt / mass;
                dt_sq_mult = .5*dt*dt / mass;
                alpha *= falpha;
                neg_alpha = -alpha;
            }
        }
        else{
            since_leq_0 = 0;
            dt *= fdec;
            dt_mult = .5 * dt / mass;
            dt_sq_mult = dt * dt * .5 / mass;
            alpha = alpha_start;
            neg_alpha = -alpha;
            cudaMemset(d_vel, 0, sizeof(double) * num_dof);
            cudaDeviceSynchronize();
        }

        //Check for RMS force convergence convergence
        if(fmag / sqrt_dof < fcut){
            break;
        }
    }

    report_data[0] = fmag / sqrt_dof;

    if(ks > 0){
        if(! check(cudaMalloc((void **) &d_spe, sizeof(double) * npairs))){
            cublasDestroy(handle);
            return false;
        }
        cudaDeviceSynchronize();
    }

    if(kb > 0){
        if(! check(cudaMalloc((void **) &d_bpe, sizeof(double) * ntriples))){
            cublasDestroy(handle);
	    if(ks > 0) cudaFree(d_spe);
            return false;
        }
        cudaDeviceSynchronize();
    }

    if(bmod > 0){
        if(! check(cudaMalloc((void **) &d_gpe, sizeof(double) * nelems))){
            cublasDestroy(handle);
	    if(ks > 0) cudaFree(d_spe);
	    if(kb > 0) cudaFree(d_bpe);
            return false;
        }
        cudaDeviceSynchronize();
    }

    //If the stretching stiffness is non-zero, compute the total stretching
    //energy.
    if(ks > 0){
        axial_energies<<<e_grid, block>>>(ks, npairs, width, d_pos, d_pair_list, d_lengths, d_e_offsets, d_spe);
	cublasDasum(handle, npairs, d_spe, 1, report_data + 1);
        cudaDeviceSynchronize();
    }
    else report_data[1] = 0;

    //If the bending energy is non-zero, compute the total bending energy.
    if(kb > 0){
        bending_energies<<<t_grid, block>>>(kb, ntriples, width, d_pos, d_triples, d_t_offsets, d_bpe);
	cublasDasum(handle, ntriples, d_bpe, 1, report_data + 2);
        cudaDeviceSynchronize();
    }
    else report_data[2] = 0;

    //Find the gel elastic energy if the moduli are non-zero
    if(bmod > 0){
        gel_energy<<<m_grid, block>>>(smod, bmod, dilation, nelems, d_pos, d_elems, d_inv_mats, d_areas, d_gpe);
	cublasDasum(handle, nelems, d_gpe, 1, report_data + 3);
        cudaDeviceSynchronize();
    }
    else report_data[3] = 0;

    if(ks > 0) cudaFree(d_spe);
    if(kb > 0) cudaFree(d_bpe);
    if(bmod > 0) cudaFree(d_gpe);

    cublasDestroy(handle);

    return true;
}

//Utility function to free a list of arrays allocated on the GPU
void free_device_arrays(vector<void *> allocated){

    for(void *ptr : allocated){
        cudaFree(ptr);
    }
}

//Attempt to allocate an array of a given size on the GPU. If successful,
//assign a specified pointer to the first byte of allocated data. If allocation
//is unsuccessful, free the accompanying list of pointers, print an error
//message, and indicate failure.
bool attempt_device_alloc(void **ptr,size_t size,vector<void *> &allocated){
    
    if(!check(cudaMalloc(ptr, size))){
        free_device_arrays(allocated);
        return false;
    }

    cudaDeviceSynchronize();
    allocated.push_back(*ptr);
    return true;
}

//Given a set of points, bonds and triples, a desired compression and shear,
//a set of point constraints, and FIRE relaxation parameters, relax a network
//to a minimum energy configuration, subject to Lee-Edwards boundary conditions,
//with an affine starting guess.
void do_relaxation_run(string net_file, string log_name, string edata_name, string disp_base, double ks, double kb, double bmod, double smod, double dilation, double fire_params[5], int nmin, double mass, double fcut, int n_steps, dim3 block, strain_range compr, strain_range shearr, bool strain_mode, bool itemize){

    //Points, edges, and bending triples
    vector<double> points, displaced_points, areas;
    vector<edge_datum> edges;
    vector<triple_datum> triples;
    vector<int> triangles, mesh_pairs;
    vector<unsigned char> fcodes;
    int iter, num_dof, steps, num_digits, disp_count = 1;
    double min_disp, max_disp, curr_disp, inc, smult, max_comp, net_bmod;
    bool success;

    //Data for carrying out calculations on the GPU
    double *d_pos, *d_forces, *d_vel, *d_lengths, *d_areas, *d_inv_mats, *d_spe;
    int *d_pairs, *d_triples, *d_triangles, *d_mpairs;
    unsigned char *d_fcodes;
    char *d_e_offsets, *d_t_offsets;

    //A running list of all currently allocated arrays on the GPU
    vector<void *> allocated;

    //Minimum and maximum y coordinates for all points in a network, network
    //area, and preferred_dilation
    double offset, minx, miny, maxx, maxy, range, area, pref_d;

    //Structures for distributing calculations of forces and energies due
    //to stretching and bending, and for applying displacements.
    dim3 e_grid, t_grid, d_grid, m_grid, mp_grid;

    //Means of collecting and reporting results
    double report_data[4];
    FILE *log_file, *edata_file;

    //Structure to accelerate vectorized calculations
    cublasHandle_t handle;

    clock_t ct;

    ct = clock();

    //Read the network from a file, obtain the vertices at the top and bottom
    //of a network, and acquire bonds, triples mesh facets, facet areas, and
    //pairs of duplicate mesh vertices that should be subject to identical
    //forces.
    read_network_mesh(net_file, offset, points, edges, triangles, areas, mesh_pairs, minx, miny, maxx, maxy);
    area = (maxx - minx) * (maxy - miny);
    offset *= sqrt(1 + dilation);

    if(kb > 0){
        triples = get_triples(points, edges, offset);
    }

    num_dof = points.size();
    range = (maxy - miny)*sqrt(1 + dilation);
    smult = strain_mode ? range : 1;
    num_digits = get_dec_digits(1 + (int) ((get<1>(shearr)-get<0>(shearr))/get<2>(shearr)));

    //Allocate GPU data

    //Positions, velocities and forces
    if(! attempt_device_alloc((void**)&d_pos,sizeof(double)*num_dof,allocated)){
        return;
    }
    
    if(! attempt_device_alloc((void**)&d_vel,sizeof(double)*num_dof,allocated)){
        return;
    }
    
    if(! attempt_device_alloc((void**)&d_forces,sizeof(double)*num_dof,allocated)){
        return;
    }

    //Information about constraints
    if(!attempt_device_alloc((void **)&d_fcodes,sizeof(char)*num_dof/2,allocated)){
        return;
    }

    //Information about bonds
    if(ks > 0){
        if(! attempt_device_alloc((void **) &d_pairs,sizeof(int)*edges.size()*2,allocated)){
            return;
        }

        if(! attempt_device_alloc((void **) &d_lengths,sizeof(double)*edges.size(),allocated)){
            return;
        }

        if(! attempt_device_alloc((void **) &d_e_offsets,sizeof(char)*edges.size(),allocated)){
            return;
        }
    }

    //Information about bending triples
    if(kb > 0){
        if(! attempt_device_alloc((void **) &d_triples,sizeof(int)*triples.size()*3,allocated)){
            return;
        }

        if(! attempt_device_alloc((void **) &d_t_offsets,2*sizeof(char)*triples.size(),allocated)){
            return;
        }
    }

    //Facet triples and areas, pre-computed inverse matrices, and pairs of 
    //companion mesh vertices
    if(! attempt_device_alloc((void**)&d_triangles,sizeof(int)*triangles.size(),allocated)){
        return;
    }
    if(! attempt_device_alloc((void**)&d_areas,sizeof(double)*areas.size(),allocated)){
        return;
    }
    if(! attempt_device_alloc((void**)&d_inv_mats,2*sizeof(double)*triangles.size(),allocated)){
        return;
    }
    if(! attempt_device_alloc((void**)&d_mpairs,sizeof(int)*mesh_pairs.size(),allocated)){
        return;
    }

    //Transfer preliminary data from the CPU to the GPU

    //Positions
    cudaMemcpy(d_pos,&points[0],sizeof(double)*num_dof, cudaMemcpyHostToDevice);

    //Bond and triple data
    if(ks > 0) copy_bond_data(edges, d_pairs, d_lengths,  d_e_offsets);
    if(kb > 0) copy_triple_data(triples, d_triples, d_t_offsets);

    //Codes indicating whether or not to constrain various DOF
    for(iter = 0; iter < points.size()/2; iter++){
        if(points[2*iter+1] < miny+FLT_TOL || points[2*iter+1] > maxy-FLT_TOL){
            fcodes.push_back(0);
        }
        else fcodes.push_back(FX_FY);
    }
    cudaMemcpy(d_fcodes, &fcodes[0], sizeof(char)*num_dof/2, cudaMemcpyHostToDevice);

    //FEM mesh information
    cudaMemcpy(d_triangles,&triangles[0],sizeof(int)*triangles.size(), cudaMemcpyHostToDevice);
    cudaMemcpy(d_areas,&areas[0],sizeof(double)*areas.size(), cudaMemcpyHostToDevice);
    cudaMemcpy(d_mpairs,&mesh_pairs[0],sizeof(int)*mesh_pairs.size(), cudaMemcpyHostToDevice);


    //Set dimensions of grids
    d_grid = dim3(num_dof / (block.x * block.y) + 1, 1);
    e_grid = dim3(edges.size() / (block.x * block.y) + 1, 1);
    if(kb > 0) t_grid = dim3(triples.size() / (block.x * block.y) + 1, 1);
    m_grid = dim3(areas.size() / (block.x * block.y) + 1, 1);
    mp_grid = dim3(mesh_pairs.size() / 2 / (block.x * block.y) + 1, 1);

    //Prepare output files
    log_file = fopen(log_name.c_str(), "w");
    edata_file = fopen(edata_name.c_str(), "w");

    //Pre-compute inverse matrices used to compute affine transformations
    //from the reference state to the deformed state of the gel.
    get_inv_mats<<<m_grid, block>>>(areas.size(), d_pos,d_triangles,d_inv_mats);

    //Apply pre-swelling, determine the energy required to deform the network,
    //and compute the bulk modulus as the energy per unit area.
    dilate<<<d_grid, block>>>(num_dof, d_pos, minx, miny, sqrt(1+dilation));

    if(ks > 0){
        if(! check(cudaMalloc((void **) &d_spe, sizeof(double)*edges.size()))){
            cerr << "The initial stretching energy could not be computed.\n";
            free_device_arrays(allocated);
            return;
        }
        cublasCreate(&handle);

        axial_energies<<<e_grid, block>>>(ks, edges.size(), offset, d_pos, d_pairs, d_lengths, d_e_offsets, d_spe);

	cublasDasum(handle, edges.size(), d_spe, 1, &net_bmod);

	net_bmod *= 2 / area / (dilation*dilation);
        cublasDestroy(handle);
	cudaFree(d_spe);
    }
    else net_bmod = 0;

    pref_d = dilation * (bmod + net_bmod) / bmod;
    cout << "Preferred dilation: " << pref_d << "\n";

    //Apply a set of compression steps, followed by a set of shear steps. Start
    //from an affine guess, then relax about this guess to a minimum energy 
    //state.

    //Apply pre-compression.
    //If the program is running in "strain mode", rescale the specified
    //displacements by the height of the network.
    min_disp = get<0>(compr) * smult;
    max_disp = get<1>(compr) * smult;
    inc = get<2>(compr) * smult;

    if(max_disp > 0){
        //Apply the initial displacement
        affine_displace<<<d_grid, block>>>(num_dof,d_pos,miny,range,0,-min_disp);
	range -= min_disp;

        for(curr_disp = min_disp; curr_disp <= max_disp + FLT_TOL; curr_disp += inc){

            //Relax vertices
	    success = run_fire_md(num_dof, edges.size(), triples.size(), areas.size(), offset, ks, kb, smod, bmod, pref_d, d_pos, d_vel, d_forces, d_pairs, d_lengths, d_e_offsets, d_triples, d_t_offsets, d_triangles, d_areas, d_inv_mats, d_mpairs, mesh_pairs.size()/2, d_fcodes, fire_params, mass, fcut, n_steps, nmin, e_grid, block, t_grid, m_grid, mp_grid, report_data, steps);

            //cout << "Relaxation finished.\n";

	    if(! success){
                cerr << "Relaxation failed.\n";
                free_device_arrays(allocated);
	    }

            //Report a measure of convergence success to the log file, and the
	    //residual strain energy.
            if(log_file != NULL){
                if(report_data[0] <= fcut){
                    fprintf(log_file, "Convergence reached: ");
	        }
                else fprintf(log_file, "Convergence not reached: ");

	        fprintf(log_file, "Steps: %d, RMS Force: %12.8le\n", steps, report_data[0]);
	    }

	    printf("Steps: %d, RMS Force: %12.8le\n", steps, report_data[0]);
	    if(edata_file != NULL){
	        fprintf(edata_file, "%12.0lf\t%12.10lf\t%12.8le\t%12.8le\t%12.8le\t%12.8le\n", 0., curr_disp, report_data[1], report_data[2], report_data[3], report_data[1]+report_data[2]+report_data[3]);
	    }

	    //If the maximum displacement has not been reached, apply an
	    //additional compression
	    if(curr_disp < max_disp){
                affine_displace<<<d_grid, block>>>(num_dof,d_pos,miny,range,0,-inc);
		range -= inc;
            }
        }
    }
    max_comp = max_disp > min_disp ? curr_disp - inc : max_disp;

    //Apply the desired combination of shear displacements.
    //Start from an affine guess, then relax about this guess to a minimum
    //energy state.
    min_disp = get<0>(shearr) * smult;
    max_disp = get<1>(shearr) * smult;
    inc = get<2>(shearr) * smult;

    if(max_disp > 0){
        //Apply the initial displacement
        affine_displace<<<d_grid, block>>>(num_dof,d_pos,miny,range,min_disp,0);
        for(curr_disp = min_disp; curr_disp <= max_disp + FLT_TOL; curr_disp += inc){

            //Relax vertices
	    success = run_fire_md(num_dof, edges.size(), triples.size(), areas.size(), offset, ks, kb, smod, bmod, pref_d, d_pos, d_vel, d_forces, d_pairs, d_lengths, d_e_offsets, d_triples, d_t_offsets, d_triangles, d_areas, d_inv_mats, d_mpairs, mesh_pairs.size()/2, d_fcodes, fire_params, mass, fcut, n_steps, nmin, e_grid, block, t_grid, m_grid, mp_grid, report_data, steps);

	    if(! success){
                cerr << "Relaxation failed.\n";
                free_device_arrays(allocated);
	    }

	    //Report the success of the fitting attempt, find final energies,
	    //copy the final positions from the GPU to the CPU, and report them
	    //to a file.
            if(log_file != NULL){
                if(report_data[0] <= fcut){
                    fprintf(log_file, "Convergence reached: ");
	        }
                else fprintf(log_file, "Convergence not reached: ");

	        fprintf(log_file, "Steps: %d, RMS Force: %12.8le\n", steps, report_data[0]);
	    }
	    printf("Steps: %d, RMS Force: %12.8le\n", steps, report_data[0]);
	    if(edata_file != NULL){
                fprintf(edata_file, "%12.8lf\t%12.8lf\t%12.8le\t%12.8le\t%12.8le\t%12.8le\n", curr_disp, max_comp, report_data[1], report_data[2], report_data[3], report_data[1]+report_data[2]+report_data[3]);
	    }

            cudaMemcpy(&points[0],d_pos,sizeof(double)*num_dof,cudaMemcpyDeviceToHost);
            report_displaced(report_name(disp_base, num_digits, disp_count++, "pos"), points);

	    //If the maximum displacement has not been reached, apply an
	    //additional shear
	    if(curr_disp < max_disp){
                affine_displace<<<d_grid, block>>>(num_dof,d_pos,miny,range,inc,0);
            }
        }
    }

    ct = clock() - ct;

    if(log_file != NULL){
        fprintf(log_file,"Elapsed time: %f seconds\n",((float)ct)/CLOCKS_PER_SEC);
    }

    free_device_arrays(allocated);
    if(log_file != NULL) fclose(log_file);
    if(edata_file != NULL) fclose(edata_file);
}

//Given a set of command line arguments, attempt to read instructions for
//relaxing one or more fiber networks. Relax each network in turn, reporting
//results to a set of files giving a measure of convergence success, residual
//strain energies, and the displacements of all points.
void relaxation_batch(string batch_file, dim3 block, bool s_mode, bool itemize){

    string base, tag;
    double ks, kb, bmod, smod, dilation;
    //Parameters for structural relxation using FIRE
    double fp[5];
    strain_range compr, shearr;
    vector<string> net_files, log_files, disp_files, edat_files;
    bool code;
    int iter, nmin, ncut;
    double mass, fcut;

    if(!parse_batch_file(batch_file, base,tag,&ks,&kb,&bmod,&smod,&dilation,fp,&nmin,&mass,&fcut,&ncut,compr,shearr)){
        cerr << "Parsing of the batch file failed.\n";
        return;
    }

    //Generate lists of dat file and output file names
    code = get_file_lists(base, "nmc", tag, net_files, log_files, disp_files, edat_files);

    if(! code){
        cerr << "File lists could not be read.\n";
        return;
    }

    for(iter = 0; iter < net_files.size(); iter++){
        do_relaxation_run(net_files[iter], log_files[iter], edat_files[iter], disp_files[iter], ks, kb, bmod, smod, dilation, fp, nmin, mass, fcut, ncut, block, compr, shearr, s_mode, itemize);
    }
}

//Process command line arguments, then hand over control to relaxation_batch
int main(int argc, char **argv){

    char c;
    int num_read;
    bool strain_mode = true, itemize = false;
    char *batch_file = NULL;
    int block_x = 128, block_y = 1;

    while((c = getopt(argc, argv, "ab:ix:y:")) != -1){
        switch(c) {
            //Use absolute, rather than relative, displacements
            case 'a':
                strain_mode = false;
            case 'b':
                batch_file = (char *) malloc((1 + strlen(optarg))*sizeof(char));
                strcpy(batch_file, optarg);
                break;
            case 'i':
                itemize = true;
                break;
	    case 'x':
                num_read = sscanf(optarg,"%d",&block_x);
                if(num_read < 1){
                    cerr << "Option \"x\" requires an integer.\n";
                }
                else if(block_x < 1){
                    cerr << "Block dimensions must be positive.\n";
                    block_x = 128;
                }
                else{
                    printf("Block x dimension: %d\n", block_x);
                }
		break;
            case 'y':
                num_read = sscanf(optarg,"%d",&block_y);
                if(num_read < 1){
                    cerr << "Option \"y\" requires an integer.\n";
                }
                else if(block_y < 1){
                    cerr << "Block dimensions must be positive.\n";
                    block_y = 1;
                }
                else{
                    printf("Block y dimension: %d\n", block_y);
                }
                break;
            case '?' :
                if(optopt == 'b'){
                    cerr << "Option \"b\" requires a file name.\n";
                    return -1;
                }
		if(optopt == 'x'){
                    cerr << "Option \"x\" requires a block dimension.\n";
                }
                if(optopt == 'y'){
                    cerr << "Option \"y\" requires a block dimension.\n";
                }
                else if(isprint(optopt)){
                    fprintf(stderr, "Unrecognized option: %c\n", optopt);
                }
                else{
                    cerr << stderr, "Unknown option character.\n";
                }
                break;
            default:
                break;
        }
    }

    if(batch_file == NULL){
        cerr << "No batch file was specified.\n";
	exit(-1);
    }

    dim3 block(block_x, block_y);

    relaxation_batch(batch_file, block, strain_mode, itemize);

    return 0;
}
