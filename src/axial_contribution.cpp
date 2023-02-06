/*
 * Given an initial network configuration, a set of displacements, and a
 * spring stiffness, this code finds the energy of deformation assoicated
 * with the stretching and compression of fibers.
 */

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <stdio.h>
#include <string.h>
#include <cmath>

using namespace std;

//Simple structure to hold information about a bond. Required quantities are
//the indices of the bond's end points, and the Cartesian coordinates of
//the bond's unit direction vector.
typedef tuple<int, int, double, double> bond_vec;
typedef tuple<double, double> vec2d;

//Read a network description file, and create a set of unit vectors
bool read_network(string file_name, int &num_pts, vector<bond_vec> &bvecs){

    string nextline;
    int num_read, num_edges, idx1, idx2;
    double x, y, dx, dy, offset, length;
    char mult;
    ifstream input;
    vector<vec2d> points;

    input.open(file_name, ifstream::in);
    if(! input.is_open()){
        cerr << "The network file name could not be read.\n";
        return false;
    }

    //Attempt to read the header for the network description file. If too 
    //little information or invalid data are found, close the input stream 
    //and report a failure to read a network.
    getline(input, nextline);
    num_read = sscanf(nextline.c_str(), "%d %d %lf", &num_pts, &num_edges, &offset);
    if(num_read < 3 || offset < 0 || num_pts <= 0 || num_edges < 0){
        cerr << "Reading failed. The provided file has an invalid header.\n";
        input.close();
        return false;
    }

    //If a valid header was read, attempt to read the specified number of 
    //points
    while(! input.eof() && points.size() < num_pts){
        getline(input, nextline);
        num_read = sscanf(nextline.c_str(), "%lf %lf", &x, &y);
        if(num_read == 2){
            points.push_back(make_tuple(x, y));
        }
    }

    //If the program has reached the end of the input stream, too little
    //information has been provided, and the file is invalid.
    if(input.eof()){
        cerr << "Reading failed. Too few points were read.\n";
        input.close();
        return false;
    }


    while(!input.eof() && bvecs.size() < num_edges){
        getline(input, nextline);
        num_read = sscanf(nextline.c_str(),"%d %d %hhd\n", &idx1, &idx2, &mult);

        //Ensure two point indices and an offset are given. If so, and 
        //the point indices are in bounds, find the length of the edge,
        //and add a record of it to the list of unit vectors.
        if(num_read == 3){
            if(min(idx1, idx2) > -1 && max(idx1, idx2) < num_pts){
                dx = get<0>(points[idx2])+offset*mult - get<0>(points[idx1]);
                dy = get<1>(points[idx2]) - get<1>(points[idx1]);
                length = sqrt(dx*dx + dy*dy);
		dx /= length;
		dy /= length;
                bvecs.push_back(make_tuple(idx1, idx2, dx, dy));
            }

            else{
                cerr << "Reading failed due to an out-of-bounds point index.\n";
                input.close();
                return false;
            }
        }
    }

    input.close();

    return bvecs.size() == num_edges;
}

//Read a displacement file
bool read_displacements(string file_name, vector<vec2d> &disps){

    double ux, uy;
    int num_read;
    ifstream input;
    string nextline;

    input.open(file_name, ifstream::in);
    if(! input.is_open()){
        cerr << "The specified displacement file could not be read.\n";
	return false;
    }

    while(! input.eof()){
	    getline(input, nextline);
            num_read = sscanf(nextline.c_str(), "%lf %lf", &ux, &uy);
	    if(num_read == 2){
                disps.push_back(make_tuple(ux, uy));
	    }
    }

    input.close();

    return true;
}

//Check that a valid network description and displacement file have been passed
//at the command line. Read in the network to obtain unit vectors, obtain
//the displacement field, and tabulate the contribution to the energy from
//stretching and compression.
int main(int argc, char **argv){

    string net_file, disp_file, usage_message;
    double k, energy, rx, ry, ux1, uy1, ux2, uy2, proj;
    int num_read, num_pts, idx1, idx2;
    vector<bond_vec> bvecs;
    vector<vec2d> disps;

    usage_message = "Usage: network file, displacement file, stiffness\n";

    //Check the input arguments. Ensure there are two strings and a
    //positive number.
    if(argc < 4){
        cerr << usage_message;
        exit(-1);
    }
    net_file = string(argv[1]);
    disp_file = string(argv[2]);
    num_read = sscanf(argv[3], "%lf", &k);
    if(num_read < 1){
	    cerr << usage_message;
	    exit(-1);
    }
    
    //Attempt to read a network and obtain bond unit vectors. Abort
    //if the file cannot be read.
    if(! read_network(net_file, num_pts, bvecs)){
        cerr << "A valid network could not be read.\n";
	exit(-2);
    }
    
    //Attempt to read displacement vectors. Abort if the file cannot be
    //read, or if too few displacement vectors are obtained.
    if(! read_displacements(disp_file, disps)){
        cerr << "A valid set of displacements could not be read.\n";
	exit(-2);
    }
    else if(disps.size() < num_pts){
        cerr << "Too few displacement vectors were read.\n";
	exit(-2);
    }

    energy = 0;
    //Tabulate energy, according to the form k/2 (u_ij dot r_ij)^2
    for(bond_vec bond : bvecs){
        idx1 = get<0>(bond);
	idx2 = get<1>(bond);
	rx = get<2>(bond);
	ry = get<3>(bond);

	ux1 = get<0>(disps[idx1]);
	uy1 = get<1>(disps[idx1]);
	ux2 = get<0>(disps[idx2]);
	uy2 = get<1>(disps[idx2]);

        proj = (ux2 - ux1)*rx + (uy2 - uy1)*ry;

	energy += .5 * proj * proj;
    }
    energy *= k;

    printf("%12.10le\n", energy);
    
    return 0;
}
