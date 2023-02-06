/*
 * This program computes a non-affine parameter for a displacement field of
 * a strained, filamentous, elastic network. The non-affine parameter 
 * characterizes the departure of a displacement field from its reference
 * value for a well-connected network that exhibits near-affine deformation
 * in response to an externally applied load, and is defined as the mean
 * square distance between the actual and reference value for the displacement
 * field at each vertex in the network, divided by the area of the network
 * and the square of the applied strain.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <tuple>
#include <map>
#include <cfloat>
#include <cstring>
#include <boost/filesystem.hpp>
#include <algorithm>
#include <sstream>
#include <boost/regex.hpp>
#include <stdio.h>

using namespace std;

//Simple structure for representing a Cartesian point or vector in 2D
#define vec2d tuple<double, double>

//Read a set of vertices from a file
void read_points(string file_name, vector<vec2d> &points, double &minx, double &miny, double &maxx, double &maxy){

	double x, y, dummy;
	ifstream file_stream;
	string nextline;
	int num_read, num_points, num_edges;

	minx = miny = FLT_MAX;
	maxx = maxy = FLT_MIN;

	file_stream.open(file_name);
	if(! file_stream.is_open()){
		cerr << "Error: Could not open " << file_name << ".\n";
		return;
	}

	getline(file_stream, nextline);
	num_read = sscanf(nextline.c_str(), "%d %d %lf\n", &num_points, &num_edges, &dummy);
	if(num_read < 3){
		cerr << "Invalid header detected.\n";
		file_stream.close();
		return;
	}

	while(! file_stream.eof() && points.size() < num_points){
		getline(file_stream, nextline);
		if((num_read = sscanf(nextline.c_str(),"%lf %lf",&x,&y)) == 2){
			points.push_back(make_tuple(x, y));
		}
		if(x < minx) minx = x;
		if(x > maxx) maxx = x;
		if(y < miny) miny = y;
		if(y > maxy) maxy = y;
	}

	file_stream.close();

	if(points.size() < num_points){
		cerr << "The indiciated number of points could not be read.\n";
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

//Given a file and directory name, create a string giving the path to a file
string get_full_name(string name, string dir){
	ostringstream full_name;
	full_name << dir << name;
	return full_name.str();
}

//Divide a string into pieces, separated by a specified delimiter
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

//Find the mean square discrepancy between a reference displacement field
//and another displacement field
/*double find_mean_discrepancy(map<vec2d, vec2d> ref, map<vec2d, vec2d> curr){

	//Total discrepancy, and components of reference and current networks'
	//displacement fields
	double discrepancy = 0, urx, ury, ucx, ucy;

	//Iterate over the vertices of the reference network. If a vertex is
	//missing from the current network, the contribution to the total
	//discrepancy should just be the squared magnitude of the displacement
	//field of the reference network at the given point. Otherwise, find
	//the squared Euclidean norm of the difference between the current and
	//reference displacement fields.
	for(auto iter = ref.begin(); iter != ref.end(); iter++){
		auto curr_val = curr.find(iter->first);
		urx = get<0>(iter->second);
		ury = get<1>(iter->second);

		if(curr_val == curr.end()){
			discrepancy += urx*urx + ury*ury;
		}
		else{
			ucx = get<0>(curr_val->second);
			ucy = get<1>(curr_val->second);
			discrepancy += (ucx-urx)*(ucx-urx)+(ucy-ury)*(ucy-ury);
		}
	}

	return discrepancy / ref.size();
}*/

//Find the minimum and maximum x and y values in a network, and use these
//to compute the area of a bounding rectangle enclosing the network.
/*double get_area(vector<vec2d> points){

	double minx, miny, maxx, maxy;
	double x, y;

	minx = miny = FLT_MAX;
	maxx = maxy = FLT_MIN;

	for(vec2d next : points){
		x = get<0>(next);
		y = get<1>(next);

		minx = x < minx ? x : minx;
		maxx = x > maxx ? x : maxx;
		miny = y < miny ? y : miny;
		maxy = y > maxy ? y : maxy;
	}

	return (maxx - minx) * (maxy - miny);
}*/

//Create a map from points to displacement field values
/*bool build_disp_map(map<vec2d, vec2d> &disp_map, vector<vec2d> points, vector<vec2d> disps){

	int iter;

	if(points.size() != disps.size()){
		cerr << "The numbers of points and displacements differ.\n";
		cerr << "Counts are " << points.size() << " and " << disps.size() << "\n";
		return false;
	}

	for(iter = 0; iter < points.size(); iter ++){
		disp_map.insert(make_pair(points[iter], disps[iter]));
	}

	return true;
}*/

//Make sure the end of a file has not been reached, then read the next line
//from that file
bool check_and_read(ifstream &input, string &line){

	if(input.eof()){
		cerr << "EOF was reached before all information was read.\n";
		input.close();
		return false;
	}

	getline(input, line);
	return true;
}

/*
 * Process a batch file with the following information
 * 	-A pattern for network displacement files to be compared to the 
 * 	reference case
 * 	-A pattern for displacement files to be comapared to the reference case
 * 	-Axial strain
 * 	-Shear strain
 */
//bool parse_batch_file(ifstream &input, string &ref_dat, string &ref_disp, string &dat_pattern, string &disp_pattern, double &strain, string &report){
bool parse_batch_file(ifstream &input, string &dat_pattern, string &disp_pattern, double &axial_strain, double &shear_strain, string &report){

	string nextline;

	//Find the name for the reference network description file and
	//displacement file
	/*if( ! check_and_read(input, ref_dat)) return false;
	if( ! check_and_read(input, ref_disp)) return false;*/

	//Obtain the patterns for batch processing of diluted network files
	//and accompanying displacement files
	if( ! check_and_read(input, dat_pattern)) return false;
	if( ! check_and_read(input, disp_pattern)) return false;

	//Obtain the strain applied to each network
	if(! check_and_read(input, nextline))	return false;
	if(sscanf(nextline.c_str(), "%lf %lf", &axial_strain, &shear_strain) < 2) return false;

	//Find the file to which calculations should be reported. If none is
	//given, indicate this with an empty string.
	if(! check_and_read(input, report)) report  = "";

	input.close();

	return true;
}

double calc_na_val(vector<vec2d> points, vector<vec2d> disps, double minx, double miny, double axial_strain, double shear_strain){

	//Coordinates and actual and expected displacement field values
	double x, y, ux_a, uy_a, ux_e, uy_e;
	double na_param  = 0;
	int iter;

	for(iter = 0; iter < points.size(); iter++){
		x = get<0>(points[iter]);
		y = get<1>(points[iter]);
		ux_a = get<0>(disps[iter]);
		uy_a = get<1>(disps[iter]);

		ux_e = (y - miny) * shear_strain;
		uy_e = (y - miny) * axial_strain;
		na_param += (ux_a-ux_e)*(ux_a-ux_e)+(uy_a-uy_e)*(uy_a-uy_e);
	}

	return na_param / points.size();
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
	//string ref_dat, ref_disp, dat_pattern, disp_pattern, report_name;
	string dat_pattern, disp_pattern, report_name;
	//map<vec2d, vec2d> disp_map;
	double axial_strain, shear_strain, minx, miny, maxx, maxy, area, na_val;
	vector<string> dat_files, disp_files;
	ifstream input;
	FILE *report;
	bool report_open;

	if(argc < 2){
		cerr << "Usage: batch file";
		return -1;
	}

	input.open(argv[1]);
	if(! input.is_open()){
		cerr << "The specified file could not be opened.\n";
		return -1;
	}

	//Obtain information about the reference network description and
	//displacement files, and the other files to compare to these reference
	//files
	//parse_batch_file(input, ref_dat, ref_disp, dat_pattern, disp_pattern, strain, report_name);
	parse_batch_file(input, dat_pattern, disp_pattern, axial_strain, shear_strain, report_name);

	//If a valid name for a report file was obtained, open this file for
	//writing. Otherwise, make a note that calculations should reported to
	//the terminal.
	if(! report_name.compare("") == 0){
		report = fopen(report_name.c_str(), "w");
	}

	/*//Extract the reference point set and displacement field
	read_points(ref_dat, points);
	area = get_area(points);
	get_displacements(ref_disp, disps);
	if(! build_disp_map(ref_map, points, disps)){
		cerr << "Building the reference map failed.\n";
		if(report) fclose(report);
		return -1;
	}*/

	dat_files = get_file_list(dat_pattern);
	disp_files = get_file_list(disp_pattern);

	//Read each dat file/displacement file pair in turn, extract the points
	//and displacement field, and compare to the reference case. 
	for(int iter = 0; iter < dat_files.size(); iter ++){
		points.clear();
		disps.clear();
		read_points(dat_files[iter], points, minx, miny, maxx, maxy);
		get_displacements(disp_files[iter], disps);
		area = (maxx-minx)*(maxy-miny);
		//curr_map.clear();
		/*if(! build_disp_map(curr_map, points, disps)){
			cerr << "Point and displacement count mismatch: ";
			cerr<<disp_files[iter]<<" and "<< dat_files[iter]<<"\n";
			continue;
		}*/
		if(points.size() != disps.size()){
			cerr << "Point and displacement count mismatch.\n";
			continue;
		}
		na_val = calc_na_val(points, disps, minx, miny, axial_strain, shear_strain);
		if(report){
			fprintf(report, "%s\t%2.10le\n", disp_files[iter].c_str(), na_val);
		}
	}

	if(report) fclose(report);

	return 0;
}
