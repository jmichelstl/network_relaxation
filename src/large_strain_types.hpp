/*
 * Author: Jonathan Michel
 * 
 * Note - This project makes use of the library libxml2, and uses boilerplate
 * code written by Daniel Veillard, who owns the copyright for libxml2 code.
 *
 * This header defines data types for representing elements of a large-strain
 * relaxation job, and implements a parser that reads such information from a
 * file.
 */

#include <vector>
#include <string>
#include <stdio.h>
#include <libxml/xmlreader.h>
#include <stack>
#include <tuple>

using namespace std;

//Codes denoting that various data types should be reported.
#define EREPORT 0x1
#define LREPORT 0x2
#define DREPORT 0x4

#define ELEM_NUM 8

//Simple data structure for attribute-value pairs
typedef tuple<char *, char *> av_pair;

//Enumerate the permissible minimization strategies. At the present, supported
//strategies are FIRE and the nonlinear Conjugate Gradient technique with
//a Fletcher-Reeve line search.
enum RelaxMethods {FIRE, ConjGrad};

const char **elem_names = {"rbatch", "name", "mech_params", "rscheme", "run", "axial", "shear", "report"};

//List all permissible names for component of a batch description file
enum rbatch_t {rbatch, name, mech_params, rscheme, run, axial, shear, report};

//List all data types that are members a given data type
rbatch_t subtypes[][] = {{name, mech_params, rscheme}, {}, {}, {}, {axial, shear, report}, {}, {}, {}};

int subt_counts[] = {3, 0, 0, 0, 3, 0, 0, 0};

//Lists of attributes and attribute counts for each member of a batch 
//description file
const char **attributes = {{}, {"input", "tag"}, {"stretch", "bend"}, {"type", "pars"}, {}, {"min", "max", "inc"}, {"min", "max", "inc"}, {"edata", "log", "disp"}};

int attr_counts[] = {0, 2, 2, 2, 0, 3, 3, 3};

//Given an element name, attempt to find it in the list of valid names. If it
//can be found, assign a reference to an instance of rbatch_t to that type.
//If the type cannot be found, indicate an error.
bool find_rbatch_type(char *elem_name, rbatch_t &rb_elem){

	int idx;

	for(idx = 0; idx < ELEM_NUM; idx++){
		if(strcomp(elem_name, elem_names[idx]) == 0){
			rb_elem = (rbatch_t) idx;
			return true;
		}
	}

	return false;
} 

//Given an element of an xml tree describing a relaxation job, and the
//next-highest hierarchical level in the tree, determine whether the current
//element is a valid subtype.
bool is_subtype(rbatch_t parent, rbatch_t child){

	int idx;

    for(idx = 0; idx < subt_counts[parent]; idx ++){
		if(subtypes[parent][idx] == child){
			return true;
		}
	}

	return false;
}

//Given an element and an attribute name, ensure the attribute name is valid
bool is_attribute(rbatch_t elem, char *attr_name){

	int idx;

	for(idx = 0; idx < attr_counts[elem]; idx ++){
		if(strcmp(attributes[elem][idx], attr_name) == 0){
			return true;
		}
	}

	return false;
}

//This type defines instructions for relaxing a network. The user must specify
//the method to be used, and parameters to be used in carrying out that method.
struct RelaxScheme {

	RelaxScheme(RelaxMethods myMethod, vector<double> myPars){
		method = myMethod;
		parameters = myPars;
	}

	RelaxMethods method;
	vector<double> parameters;
};

//This type defines a relaxation runs. A job file may prescribe multiple
//relaxation runs in succession.
struct Run {

	Run(){}

	/*Run(double axMin, double axMax, double axInc, double sMin, double sMax, double sInc, char repFlags){

		axial_min = axMin;
		axial_max = axMax;
		axial_inc = axInc;
		shear_min = sMin;
		shear_max = sMax;
		shear_inc = sInc;
		report_flags = repFlags;
	}*/

	//Data types representing minima, maxima, and increments for axial and
	//shear deformations of networks
	double axial_min, axial_max, axial_inc, shear_min, shear_max, shear_inc;

	//Flags indicating whether energy files, log files, and displacement field
	//files should be reported.
	char report_flags;
};

//This type defines an overall job file. It should specify an input file
//pattern, material properties, a relaxation scheme, and one or more relaxation
//runs. The user may (and should) specify a tag to be appended to output file
//names along with the input file pattern.
struct RelaxBatch {

	string input_pattern, tag;
	double mech_pars[2];
	RelaxScheme scheme;
	vector<Run> runs;

	RelaxBatch(){}

	/*RelaxBatch(string pat, string myTag,double mpars[2],RelaxScheme myScheme){
		input_pattern = pat;
		tag = myTag;
		
		runs = vector<Run>();
	}*/

	void set_mech_pars(double stretch, double bend){
		mech_pars[0] = stretch;
		mech_pars[1] = bend;
	}

	void add_run(Run nextRun){
		runs.push_back(nextRun);
	}

};

//Obtain information about the current node. In the process, ensure the current
//node obeys all rules for the formatting of relaxation batch files.
bool processNode(xmlTextReaderPtr reader, stack<rbatch_t> &type_stack) {
    
	const xmlChar *name, *value;
	vector<av_apair> pairs;

    name = xmlTextReaderConstName(reader);
    if (name == NULL){
		cerr << "A nameless tag was encountered.\n";
		return false;
	}

    if(name[0] == '#') return true;

	//If a node has a depth that is less than or equal to the previous depth,
	//the parsing context should be updated. If the node is a closing element,
	//all that should be done is to confirm it matches an opening tag.

	//Obtain any attributes, and make a list of attribute-value pairs
    if( xmlTextReaderHasAttributes(reader)){
	    xmlTextReaderMoveToFirstAttribute(reader);

    	do{
        	name = xmlTextReaderConstName(reader);
	    	if(name == NULL) name = "NULL";

			value = xmlTextReaderConstValue(reader);
			if(value == NULL) value = "NULL";

			pairs.push_back(make_tuple((char *) name, (char *) value));
    	}while(xmlTextReaderMoveToNextAttribute(reader) > 0);
	}

	
}

bool streamRelaxBatchFile(const char *filename) {

    xmlTextReaderPtr reader;
    int ret;
	stack<rbatch_t> type_stack;

	reader = xmlReaderForFile(filename, NULL, 0);
    if (reader != NULL) {

        ret = xmlTextReaderRead(reader);
        while (ret == 1) {
            if(! processNode(reader)){
				xmlFreeTextReader(reader);
				return false;
			}

            ret = xmlTextReaderRead(reader, typestack);
        }
    
	    xmlFreeTextReader(reader);
        if (ret != 0) {
            fprintf(stderr, "%s : failed to parse\n", filename);
			return false;
        }

		return true;
    } 

	else {
        fprintf(stderr, "Unable to open %s\n", filename);
		return false;
    }
}


/* Given an input file, attempt to read the data necessary for specifying a
 * batch of relaxation jobs from an XML-formatted file. The overall type of
 * the document should be relaxation_run. The acceptable fields at the next
 * level are name, mech_params, relax_scheme and run. The fields supported
 * by these categories are as follows:
 *
 * name:
 *	-input
 *	-tag
 *
 * mech_params
 *	-stretch
 *	-bend
 *
 * relax_scheme
 *	-method (FIRE or ConjGrad)
 *	-relax_pars (string of space-delimited doubles)
 *
 * run
 *	-axial:		Arguments are min, max and increment
 *	-shear: 	Arguments are min, max and increment
 *	-report:	Arguments are edata = true/false, log = true/false, 
 *			 	disp = true/false
 *
 * Return value:
 *	-True, if the batch file was parsed successfully
 *	-False, if the batch file was not parsed successfully
 */
bool parse_relaxation_batch(const char *filename, RelaxBatch &rbatch){



	return true;
}
