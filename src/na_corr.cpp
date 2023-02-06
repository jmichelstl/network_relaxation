/*
 * This program computes correlations in the non-affine displacements of
 * vertices in a strained elastic network. Given a reference network 
 * configuration, a set of displacements, and an imposed boundary strain, the
 * program computes the deviation from affine displacement for each vertex
 * in the network, and takes the dot product of this deviation with the
 * deviation of a set of nieghboring points within a given cutoff radius.
 * Correlations are binned according to some radius step size and azimuthally
 * averaged, and are finally normalized by the square of the mean non-affine
 * deviation. An accompanying power spectrum is finally produced.
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
#include <boost/tokenizer.hpp>
#include <sstream>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <cmath>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include <fftw3.h>
#include <stdlib.h>
#include <unistd.h>

typedef boost::tokenizer<boost::char_separator<char>> tokenizer;

using namespace std;

//Shorthand for namespaces used for spatial indexing
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

//Shorthand for types used in spatial indexing
typedef bg::model::point<double, 2, bg::cs::cartesian> bpoint;
typedef bg::model::box<bpoint> box;
typedef std::pair<bpoint, unsigned> value;

//Simple structure for representing a Cartesian point or vector in 2D
#define vec2d tuple<double, double>

#define FLOAT_TOL 1e-8

//This is a comparator to sort Cartesian points. Points are sorted first by
//x index, and then by y index, with a non-zero floating point tolerance.
/*struct pcomp {
	bool operator() (const vec2d &p1, const vec2d &p2) const {
		double x1, x2, y1, y2;
		x1 = get<0>(p1);
		x2 = get<0>(p2);
		y1 = get<1>(p1);
		y2 = get<1>(p2);

		if(y1 < y2 - FLOAT_TOL) return true;
		else if(y2 < y1 - FLOAT_TOL) return false;
		else if(x1 < x2 - FLOAT_TOL) return true;
		else return false;
	}
};*/

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
bool parse_batch_file(ifstream &input, string &dat_pattern, string &disp_pattern, double &axial_strain, double &shear_strain, double &bin_size, double &cutoff){

	string nextline;

	//Obtain the patterns for batch processing of diluted network files
	//and accompanying displacement files
	if( ! check_and_read(input, dat_pattern)) return false;
	if( ! check_and_read(input, disp_pattern)) return false;

	//Obtain the strain applied to each network
	if(! check_and_read(input, nextline))	return false;
	if(sscanf(nextline.c_str(), "%lf %lf", &axial_strain, &shear_strain) < 2) return false;

	//Obtain the bin size and the maxium separation between points to be
	//considered
	if(! check_and_read(input, nextline))	return false;
	if(sscanf(nextline.c_str(), "%lf %lf", &bin_size, &cutoff) < 2) return false;

	input.close();

	return true;
}

vector<vec2d> calc_na_disps(vector<vec2d> points, vector<vec2d> disps, double minx, double miny, double axial_strain, double shear_strain){

	//Coordinates and actual and expected displacement field values
	double x, y, ux_a, uy_a, ux_e, uy_e;
	double na_param  = 0;
	int iter;
	vector<vec2d> na_disps;

	for(iter = 0; iter < points.size(); iter++){
		x = get<0>(points[iter]);
		y = get<1>(points[iter]);
		ux_a = get<0>(disps[iter]);
		uy_a = get<1>(disps[iter]);

		ux_e = (y - miny) * shear_strain;
		uy_e = (y - miny) * axial_strain;

		na_disps.push_back(make_tuple(ux_a - ux_e, uy_a - uy_e));
	}

	return na_disps;
}

//Given a displacement file, prepare companion files for the non-affine
//correlation function and the accompanying power spectrum
vector<string> get_corr_ps_names(string disp_name){

	ostringstream base_name, corr_oss, pow_oss;
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

	if(disp_name[0] == '.') base_name << ".";
	for(int i = 0; i < tokens.size() - 1; i++){
		base_name << tokens[i] << ".";
	}

	//Create the correlation file name by appending the .corr file
	//extension.
	corr_oss << base_name.str() << "corr";
	file_names.push_back(corr_oss.str());

	//Create the power spectrum file by appending the .pow extension
	pow_oss << base_name.str() << "pow";
	file_names.push_back(pow_oss.str());

	return file_names;
}

//Given a displacement file, prepare a companion file for the non-affine
//correlation map.
string get_nacm_name(string disp_name){

	ostringstream base_name, corr_oss;
	boost::char_separator<char> sep{"."};
	tokenizer tok{disp_name, sep};
	vector<string> tokens;

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

	//Create the non-affine correlation map file name by appending the 
	//.nacm file extension.
	corr_oss << base_name.str() << "nacm";
	return corr_oss.str();
}

//Given a set of displacements, return the indices of a set of points whose
//displacement magnitudes are in a given bottom fraction of all displacement
//magnitudes.
vector<int> least_n_disps(vector<vec2d> disps, double fraction){

	vector<int> keepers;
	double x, y, max;

	if(abs(fraction - 1) < FLOAT_TOL){
		for(int iter = 0; iter < disps.size(); iter++){
			keepers.push_back(iter);
		}
	}

	else{
		vector<double> lengths_sq;
		for(vec2d next_disp : disps){
			x = get<0>(next_disp);
			y = get<1>(next_disp);
			lengths_sq.push_back(x*x + y*y);
		}
		vector<double> duplicates(lengths_sq.begin(), lengths_sq.end());
		sort(duplicates.begin(), duplicates.end());
		max = duplicates[ceil(fraction * lengths_sq.size())];

		for(int iter = 0; iter < disps.size(); iter++){
			if(lengths_sq[iter] <= max){
				keepers.push_back(iter);
			}
		}
	}

	return keepers;
}

//Calculate the radially averaged correlation for non-affine displacement
void calc_na_corr(vector<vec2d> points, vector<vec2d> na_disps, double cutoff, double inc, vector<double> &corrs, vector<int> &counts, double fraction){

	double dx, dy, r_sq, rcut_sq, nax1, nay1, nax2, nay2;
	int idx2, c_idx;
	rcut_sq = cutoff*cutoff;
	vector<int> keepers;

	//Initialize the lists of correlations and bin counts with zeros
	corrs.assign(ceil(cutoff / inc) + 1, 0);
	counts.assign(ceil(cutoff / inc) + 1, 0);

	//Construct an R Tree to enable consideration only of points within
	//the cutoff radius
	keepers = least_n_disps(na_disps, fraction);
	bgi::rtree<value, bgi::quadratic<16>> rtree;
	for(int iter : keepers){
		rtree.insert(make_pair(bpoint(get<0>(points[iter]), get<1>(points[iter])), iter));
	}	

	//cout << "Running here.\n";

	//Iterate over each point, find those neighbors within the cutoff
	//distance, and with a higher index in the array of points, and
	//compute the dot product of the two points' non-affine displacements.
	//Also keep a running sum of the squared magnitude of all non-affine
	//displacements for normalization.
	for(int idx1 : keepers){

		//Select points within a box circumscribing a circle of the
		//cutoff radius.
		vec2d p = points[idx1];
		box query_box(bpoint(get<0>(p)-cutoff, get<1>(p)-cutoff), bpoint(get<0>(p)+cutoff, get<1>(p)+cutoff));
		vector<value> results;
		rtree.query(bgi::within(query_box), back_inserter(results));

		nax1 = get<0>(na_disps[idx1]);
		nay1 = get<1>(na_disps[idx1]);

		for(value v : results){
			idx2 = v.second;

			//Avoid redundant calculations.
			if(idx2 < idx1) continue;

			dx = v.first.get<0>() - get<0>(p);
			dy = v.first.get<1>() - get<1>(p);
			r_sq = dx*dx + dy*dy;

			if(r_sq <= rcut_sq + FLOAT_TOL){
				c_idx = floor((sqrt(r_sq) + inc / 2) / inc);
				nax2 = get<0>(na_disps[idx2]);
				nay2 = get<1>(na_disps[idx2]);
				corrs[c_idx] += nax1*nax2 + nay1*nay2;
				counts[c_idx] ++;
			}
		}
	}
	//cout << "Running here.\n";

	//Take the average of the correlation function for each radius, and
	//normalize all values by the value at zero separation. If a bin
	//has no entries, forgo taking the average.
	corrs[0] /= counts[0];

	for(int iter = 1; iter < counts.size(); iter++){
		if(counts[iter] > 0){
			corrs[iter] = corrs[iter] / counts[iter] / corrs[0];
		}
	}

	corrs[0] = 1;
}

//Calculate the radially averaged correlation for non-affine displacement
void na_corr_no_avg(vector<vec2d> points, vector<vec2d> na_disps, double cutoff,double inc, vector<double> &corrs, vector<vec2d> &separations, double fraction){

	double dx, dy, r_sq, rcut_sq, nax1, nay1, nax2, nay2, corr, xdiff,ydiff;
	int iter, idx2, dcount = 1, disp_index;
	rcut_sq = cutoff*cutoff;
	vector<int> counts, to_duplicate, keepers;
	vector<vec2d> reflections;
	vector<double> duplicates;
	bool match;

	//R Tree to search for displacement vectors that have already been
	//encountered to ensure each unique displacement vector is only
	//considered once.
	bgi::rtree<value, bgi::quadratic<16>> disp_rtree;

	//Construct an R Tree to enable consideration only of points within
	//the cutoff radius
	keepers = least_n_disps(na_disps, fraction);
	bgi::rtree<value, bgi::quadratic<16>> point_rtree;
	for(int next_keeper : keepers){
		point_rtree.insert(make_pair(bpoint(get<0>(points[next_keeper]), get<1>(points[next_keeper])), next_keeper));
	}

	separations.push_back(make_tuple(0, 0));
	corrs.push_back(0);
	counts.push_back(points.size());

	//Iterate over each point, find those neighbors within the cutoff
	//distance, and with a higher index in the array of points, and
	//compute the dot product of the two points' non-affine displacements.
	//Also keep a running sum of the squared magnitude of all non-affine
	//displacements for normalization.
	for(int idx1 : keepers){

		//Select points within a box circumscribing a circle of the
		//cutoff radius.
		vec2d p = points[idx1];
		box query_box(bpoint(get<0>(p)-cutoff, get<1>(p)-cutoff), bpoint(get<0>(p)+cutoff, get<1>(p)+cutoff));
		vector<value> results;
		point_rtree.query(bgi::within(query_box), back_inserter(results));

		nax1 = get<0>(na_disps[idx1]);
		nay1 = get<1>(na_disps[idx1]);

		for(value v : results){
			idx2 = v.second;

			//Only consider each distinct pair of neighbors once.
			if(idx2 < idx1) continue;

			//When computing the autocorrelation for zero
			//displacement, a simplified procedure may be used.
			if(idx1 == idx2){
				corrs[0] += nax1*nax1 + nay1*nay1;
				continue;
			}

			dx = v.first.get<0>() - get<0>(p);
			dy = v.first.get<1>() - get<1>(p);
			r_sq = dx*dx + dy*dy;

			if(r_sq <= rcut_sq + FLOAT_TOL){
				nax2 = get<0>(na_disps[idx2]);
				nay2 = get<1>(na_disps[idx2]);
				corr = nax1*nax2 + nay1*nay2;

				//Add a record of the two-point autocorrelation
				//for the current displacement. If the 
				//displacement vector is above the line y = x,
				//add a record of it directly. Otherwise, add 
				//a record of the displacement's reflection 
				//about the line y = x.
				if(dy < dx - FLOAT_TOL){
					dy = -dy;
					dx = -dx;
				}

				//Query the R Tree storing displacement vectors
			        //for the nearest displacement vector to the
			        //current one. If the Euclidean distance is
			        //greater than the floating point tolerance,
			        //a distinct displacement has been found.
			        match = false;
			        vector<value> closest;
			        disp_rtree.query(bgi::nearest(bpoint(dx, dy), 1), std::back_inserter(closest));
			        if(closest.size() > 0){
		                	xdiff = dx - closest[0].first.get<0>();
	                		ydiff = dy - closest[0].first.get<1>();
        			        if(xdiff*xdiff+ydiff*ydiff < FLOAT_TOL){
			                        match = true;
                        			disp_index = closest[0].second;
			                        corrs[disp_index] += corr;
                        			counts[disp_index] ++;
                			}
        			}

			        //If there is no displacement vector suitably
			        //close to the current one in the R Tree,
			        //the current displacement vector, and make
			        //new records for the non-affine correlation
			        //sum and the number of times the displacement
			        //vector has been encountered.
			        if(! match){
			                disp_rtree.insert(make_pair(bpoint(dx, dy), dcount++));
			                corrs.push_back(corr);
			                counts.push_back(1);
			                separations.push_back(make_tuple(dx,dy));
			        }


			}
		}
	}

	//Take the average of the correlation function for each separation.
	for(iter = 0; iter < corrs.size(); iter++){
		corrs[iter] /= counts[iter];
	}

	//Add a record for the reflection of each displacement other those
	//vectors for which y = x about the line y = x.
	for(iter = 1; iter < corrs.size(); iter++){
		dy = -get<1>(separations[iter]);
		dx = -get<0>(separations[iter]);
		if(abs(dy - dx) < FLOAT_TOL) continue;
		
		reflections.push_back(make_tuple(dx, dy));
		to_duplicate.push_back(iter);
	}
	
	for(int next_idx : to_duplicate){
	       	duplicates.push_back(corrs[next_idx]);
	}

	corrs.insert(corrs.end(), duplicates.begin(), duplicates.end());
	separations.insert(separations.end(), reflections.begin(), reflections.end());

}

//If no values were obtained for a certain range of radii, fill in the gap by
//approximating the correlation function using cubic basis splines.
void fill_in_values(vector<double> &corrs, vector<int> &counts, double rmax, double rinc, int nnz){

	size_t i, j, nsplines, nbreak, next_pos;
	gsl_bspline_workspace *bw;
	gsl_vector *y, *c, *w, *s_vals;
	gsl_matrix *X, *cov;
	gsl_multifit_linear_workspace *mw;
	double yerr, chisq;

	//Since cubic splines are desired, basis splines of order 4 are to
	//be constructed. For nsplines spline functions, nsplines + 2 - k,
	//or nsplines - 2 break points are needed.
	nsplines = ceil(rmax / rinc / 10);
	nbreak = nsplines - 2;

	//Set up workspace
	y = gsl_vector_alloc(nnz);
	X = gsl_matrix_alloc(nnz, nsplines);
	c = gsl_vector_alloc(nsplines);
	w = gsl_vector_alloc(nnz);
	cov = gsl_matrix_alloc(nsplines, nsplines);
	bw = gsl_bspline_alloc(4, nbreak);
	s_vals = gsl_vector_alloc(nsplines);

	//Construct a set of cubic basis splines defined on the interval
	//(0, rmax), with knots spaced by rinc
	gsl_bspline_knots_uniform(0, rmax, bw);

	//Evaluate the basis splines at the radii corresponding to non-empty
	//bins, and use these to construct a matrix to fit for spline
	//coefficients.
	next_pos = 0;
	for(i = 0; i < counts.size(); i++){
		if(counts[i] > 0){
			gsl_vector_set(y, next_pos, corrs[i]);
			gsl_bspline_eval(rinc*next_pos, s_vals, bw);

			for(j = 0; j < nsplines; j++){
				gsl_matrix_set(X,i,j,gsl_vector_get(s_vals, j));
			}

			next_pos ++;
		}
	}
	
	//Use chi squared optimization to find coefficient values that yield
	//the best agreement between the spline approximation and the 
	//correlation function values.
	gsl_multifit_wlinear(X, w, y, c, cov, &chisq, mw);
	
	//Evaluate the approximation using the optimized coeffcients to
	//estimate values for g(r) at empty bins.
	for(i = 0; i < counts.size(); i++){
		if(counts[i] == 0){
			gsl_bspline_eval(i*rinc, s_vals, bw);
			gsl_multifit_linear_est(s_vals, c, cov, &corrs[i],&yerr);
		}
	}

	gsl_bspline_free(bw);
	gsl_vector_free(s_vals);
	gsl_vector_free(y);
	gsl_matrix_free(X);
	gsl_vector_free(c);
	gsl_vector_free(w);
	gsl_matrix_free(cov);
	gsl_multifit_linear_free(mw);
}

//Given a set of correlation data, find the discrete Fourier transform, and
//compute the squared magnitude of each Fourier coefficient to obtain a
//power spectrum of the correlation function.
void calc_power_spectrum(vector<double> corrs, vector<double> &pspec_vals, double *in, fftw_complex *out, fftw_plan *p){

	int i;

	//Populate the input array with correlation values
	for(i = 0; i < corrs.size(); i++){
		in[i] = corrs[i];
	}

	//Perform the transform and write the squared magnitudes of the
	//transform values to the output vector.
	fftw_execute(*p);

	for(i = 0; i < pspec_vals.size(); i++){
		pspec_vals[i] = out[i][0]*out[i][0] + out[i][1]*out[i][1];
	}
}

/*
 * Obtain a set of reference configuations and deformation fields for elastic
 * networks. For each network, calculate the non-affine displacement field
 * at each point, find correlations between points of varying separation
 * distances, and Fourier transform correlation functions to aid in identifying
 * a characteristic length scale.
 */
int main(int argc, char **argv){

	vector<vec2d> points, disps;
	string dat_pattern, disp_pattern, report_base;
	double axial_strain, shear_strain, bin_size, cutoff, fraction = 1;
       	double minx, miny, maxx, maxy, waveno_inc;
	vector<string> dat_files, disp_files, report_files;
	ifstream input;
	FILE *corr_report = NULL, *pow_report = NULL;
	vector<vec2d> na_disps, separations;
	vector<double> corrs, pspec;
	vector<int> counts;
	size_t nbins;
	int nnz, num_read;

	//Data for processing a command line switch
	bool az_avg = true;
	char c;

	//Structures for fast fourier transforms
	double *in;
	fftw_complex *out;
	fftw_plan plan;

	if(argc < 2){
		cerr << "Usage: batch file";
		return -1;
	}

	/*for(int next_arg = 1; next_arg < argc; next_arg ++){
		if(argv[next_arg][0] != '-'){
			input.open(argv[next_arg]);
			if(! input.is_open()){
				cerr << "The specified file could not be opened.\n";
				return -1;
			}
		}
	}*/

	//Process an optional command line switch. The -n flag indicates the
	//program should not perform an azimuthal average.
	while((c = getopt(argc, argv, "b:nf:")) != -1){
		switch(c) {
			case 'b':
				input.open(optarg);
				if(! input.is_open()){
					cerr << "The specified file could not be opened.\n";
					return -1;
				}

				else break;
			case 'f':
				num_read = sscanf(optarg, "%lf", &fraction);

				if(num_read < 1){
					cerr << "The -f option requires a number.\n";
					fraction = 1;
				}
				else if(fraction < 0 || fraction > 1){
					cerr << "The fraction must be between 0 and 1.\n";
					fraction = 1;
				}

				break;
			case 'n':
				az_avg = false;
				break;
			case '?' :
				if(optopt == 'b'){
					cerr << "Option -b requires a file name.\n";
				}
				else if(optopt == 'f'){
					cerr << "Option -f requires a number between 0 and 1.\n";
				}
				else if(isprint(optopt)){
					fprintf(stderr, "Unrecognized option: %c\n", optopt);
				}
				else{
					cerr << "Unknown option character.\n";
				}
				break;
			default:
				break;
		}
	}

	if(! input.is_open()){
		cerr << "No batch file was opened for reading.\n";
		return -1;
	}

	//Obtain information about the reference network description and
	//displacement files, and the other files to compare to these reference
	//files
	//parse_batch_file(input, ref_dat, ref_disp, dat_pattern, disp_pattern, strain, report_name);
	parse_batch_file(input, dat_pattern, disp_pattern, axial_strain, shear_strain, bin_size, cutoff);

	dat_files = get_file_list(dat_pattern);
	disp_files = get_file_list(disp_pattern);

	if(az_avg){
		//Adjust cutoff so that it is an integer multiple of the bin 
		//size
		cutoff = floor(cutoff / bin_size) * bin_size;
		nbins = cutoff / bin_size + 1;
		pspec.assign(nbins / 2 + 1, 0);
		waveno_inc = 2*M_PI / bin_size / nbins;

		//For each set of correlations, the same radius bins will be 
		//used, so a plan and a set of input and output arrays can be 
		//recycled.
		in = (double *) malloc(sizeof(double) * nbins); 
		out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(nbins/2 + 1));
		plan = fftw_plan_dft_r2c_1d(nbins, in, out, FFTW_MEASURE);
	}

	//Read each dat file/displacement file pair in turn, extract the points
	//and displacement field, and compare to the reference case. 
	for(int iter = 0; iter < dat_files.size(); iter ++){
		points.clear();
		disps.clear();
		na_disps.clear();
		corrs.clear();
		separations.clear();
		read_points(dat_files[iter], points, minx, miny, maxx, maxy);
		get_displacements(disp_files[iter], disps);

		if(points.size() != disps.size()){
			cerr << "Point and displacement counts do not match.\n";
			cerr << dat_files[iter] << ",\t" << disp_files[iter] << "\n";
			continue;
		}

		if(az_avg){

			report_files = get_corr_ps_names(disp_files[iter]);

			//Open correlation and power spectrum files for 
			//reporting. If either can't be opened, don't bother
			//doing anything else.
			corr_report = fopen(report_files[0].c_str(), "w");
			if(! corr_report){
				cerr << "The correlation file could not be opened.\n";
				cerr << "Name: " << report_files[0].c_str() << "\n";
				continue;
			}
			pow_report = fopen(report_files[1].c_str(), "w");
			if(! pow_report){
				cerr <<"The power spectrum file could not be opened.\n";
				if(corr_report) fclose(corr_report);
				continue;
			}
		}

		else{
			corr_report = fopen(get_nacm_name(disp_files[iter]).c_str(), "w");
			if(! corr_report){
				cerr << "The correlation file could not be opened.\n";
				continue;
			}	
		}

		//Find the non-affine displacements, then extract correlations
		//and, if azimuthal averaging is to be used, compute the power 
		//spectrum, and report correlations and the power spectrum to 
		//files.
		na_disps = calc_na_disps(points, disps, minx, miny, axial_strain, shear_strain);

		//If azimuthal averaging is to be used, compute the correlation
		//between non-affine displacements of points as a function of
		//a single variable - the separation between points.
		if(az_avg){
			counts.clear();
			calc_na_corr(points, na_disps, cutoff, bin_size, corrs,counts, fraction);
			nnz = 0;

			for(int count : counts) if(count > 0) nnz++;
			if(nnz < counts.size()){
				fill_in_values(corrs, counts, cutoff, bin_size, nnz);
			}

			calc_power_spectrum(corrs, pspec, in, out, &plan);

			for(int list_iter = 0; list_iter < nbins; list_iter++){
				fprintf(corr_report, "%lf\t%lf\n", list_iter*bin_size, corrs[list_iter]);
			}
			fclose(corr_report);

			for(int list_iter = 0; list_iter < pspec.size(); list_iter++){
				fprintf(pow_report, "%lf\t%lf\n", waveno_inc*list_iter, pspec[list_iter]);
			}
			fclose(pow_report);
		}

		//If no azimuthal averaging is to be performed, find the
		//non-affine correlation as a function of both x and y
		//separation.
		else{
			na_corr_no_avg(points, na_disps, cutoff, bin_size, corrs, separations, fraction);
			for(vec2d next_sep : separations){
				fprintf(corr_report, "%lf\t%lf\n", get<0>(next_sep), get<1>(next_sep));
			}
			for(double corr : corrs){
				fprintf(corr_report, "%lf\n", corr);
			}
			fclose(corr_report);
		}
	}

	if(az_avg){
		free(in);
		fftw_free(out);
		fftw_destroy_plan(plan);
	}

	return 0;
}
