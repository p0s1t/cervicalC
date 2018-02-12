 // NOTE: see what space time complexity differences there are when 
 // parallelizing floating point calculations or long data structure calculations 

#pragma once

#typedef struct struct px
#typedef struct struct pair
#typedef struct struct node // pixel with affected chromicity
#typedef struct struct axon
#typedef struct struct dendrite // linear mapping of edges between nodes
#typedef struct struct quadrant // quadrant within the original image, the dimensions are natural number factors of its total length
#typedef struct struct lobe // all quadrants with affected nodes above confidence threshold

class: cure
{

	public:

	virtual cure(){};
	virtual ~cure(){};

	protected:

///////////////////////////  BEGIN MEMBER VAR DATA /////////////////////////////////


////// RUN CLOCK VARS ////////////////

/*
	clock_t * begin;
	clock_t * end;
	
	long * runtime;
*/	

////////// FILE ADDRESS VARS //////////////////////


	char * i;// input filepath
	char * o;// output filepath
	char * log; // prolog filepath
	

//// ORIGINAL IMAGE DATA ////////////////////
	

	const char * m; // an original piece of the image data in byte form
	int * n; // number of bytes in total image

	const char * byt; // the image in extended rgb data form 
	int * nq; // number of quadrants in this cure

	px ** pxls; // rgb data of total image

	px * cpx; // pixels of the cervix image
	
	int * cSz; // size of cpx pixel array 

	const bool * class; // classification of image (if known)
	
	const int * nPx; // nPixels in pruned image

	pair * factors; // pairs of factors of original image length (multiply by respective set of factors for hamming string lengths) 

	long * fitness; // respective accuracy of prediction, 0 being perfectly accurate for all test samples
	
	//const int * maxRadius; // maximum radius from center for any cluster in this cure
	//const int * minRadius; // min

	//int * avRadius; // average radius of clusters in this cure


////////////// machine learning variables

	int * nFactors; 
	//px * pruned; // array of rgb data to be pruned from original image
	//char * saved; // original bytes saved from original image
	//int * r; // max radius for nearest neighbor calculations
	int * ne; // nearest neighbors
	int * nNeurons;
	// dendrite * d;
	// axon * a; 
	lobe * cortex;

//////////////////// class FUNCTIONS 


	pair * factorize();
	
	void bits(); // get bit data from image
	
	px * pixel(char * c); // convet byte to chromic data

	//void pro(); // make prolog knowledge base with self-modulating lobed data structure for python sdk

	void lobe(); // calculate classification accuracy once arrays have been run through python deep learning sdk
							// do so by running a command line script that compiles this code, then calls the sdk
	
	void log(); // make a log (legend) for reference -- runtime, classification accuracy for each image, total space complexity
				// HPC utility strengths, etc

	void byte(); // make an array from all of the original image data

	void aggregate();

	px * getPxData(int *, int *, int *, int *); 

	quadrant  * quadrant(quadrant *, int *, int *, int *); 

	node * iNode(px *, int *);

	pair * iPair(int *, int * , int *);



//////////////// global hamming variables	


//	struct ham ** hams; // hamming strings between all nodes in all quadrants, organized by quadrant index 


///////////////// STRUCT DEFINITIONS ////////////////////////////


struct pair()
{
	int * i;
	int * f1; 
	int * f2;
}; // pair



struct dendrite(){

	px * e;
	int * qi;
	vertex * v1;
	vertex * v2;

}; // linear hamming distance representation between nodal data set centers


struct axon(){

	px * v; // pixel data
	vertext * v2; // vertext at other end of edge e 
	edge * e; // edge this vertex is connected to
	int * qi; // quadrant index
	int * iN; // nodal index
	int * ii; // index of v in image

}; // vertex of a hamming distance set of pixels' data


struct ham(int * x1, int * x2, int * y1, int * y2)
{

		int * len = int) );
		len = getLen(x1, x2, y1, y2);

		px * pxData = malloc(len, sizeof(px) );
		pxData = getPxData(x1, x2, y1, y2); 

}; // linear hamming strings within the image 




struct px()
{

		int * iN; // index of node for which this color data is connected
		int * qi; // index of cluster of the node
		int * ii; // index in linear byte string from original image's data structure'

		int * x; // row
		int * y; // column
		
		short * Y; // the luminesence (brightness) of the node
		short * Cr; // the reddish chromesence (interesting per Hubble's theory deduced when observing red shift galactic movement)
		short * Cb; // the blueish chromesence
		
		short * pur; // red - blue
		short * pink; // brightness - red
		short * sky; // brightness - blue

		//dendrite * e; 
		//axon * v; 

		bool * affected; // does it meet the criteria for being effectively classifiable in first, second or third stage
		// 0 is first, 1 is second, 2 is third

}; // pixel struct


struct node(){

	px * center; 
	px * ne;

	short * avY;
	short * avCr;
	short * avCb;

	long * avPink;
	long * avPur;
	long * avSky;

	long * sdY;
	long * sdCr;
	long * sdCb;

	long * sdPink;
	long * sdPur;
	long * sdSky;

	long * sd2Y;
	long * sd2Cr;
	long * sd2Cb;

	long * sd2Pink;
	long * sd2Pur;
	long * sd2Sky;

	long * sd3Y;
	long * sd3Cr;
	long * sd3Cb;

	long * sd3Pink;
	long * sd3Pur;
	long * sd3Sky;

	long * sd4Y;
	long * sd4Cr;
	long * sd4Cb;

	long * sd4Pink;
	long * sd4Pur;
	long * sd4Sky;

	short * minY;
	short * minCr;
	short * minCb;

	long * minPink;
	long * minPur;
	long * minSky;

	short * maxY;
	short * maxCr;
	short * maxCb;

	long * maxPink;
	long * maxPur;
	long * maxSky;

	
}; // node


struct quadrant()
{	

	short * minY;
	short * minCr;
	short * minCb;

	long * minPur;
	long * minPink;
	long * minSky;	

	long * maxY;
	long * maxCr;
	long * maxCb;
	
	long * maxPur;
	long * maxPink;
	long * maxSky;

	long * avY;
	long * avCr;
	long * avCb;
	
	long * avPur;
	long * avPink;
	long * avSky;

	// position relative to center of image

	bool * bRight;
	bool * tRight;

	bool * bLeft;
	bool * tLeft;

	int * n; // n pixels in quadrant
	
	int * area; 	
	int * row;
	int * col;

	px * p; // pixels in quadrant

	bool * isEcto; // classification data
	bool * isEndo;

	bool * isTopVis; // is it edgewise bound?

	double * confidence;

	int * cy1; // in ascending order 
	int * cx1; // the coordinates that contain pixels imaging the cervix proper

	int * cy2; // these are the end points of the cervix, so x coordinates go from cx1 to cx2, y from cy1 to cy2
	int * cx2;

	int * ne;

	long * nnRadius; // max radius when calculating nearest neighbors
	
	ham * inter; // hamming distances between geometric center and its respective pxls
	ham * intra;

	px * node; // nearest neighbors

}; // quadrant



// quadrant data




   // hamming functions
	
	double * getAvProximity(quadrant * q, quadrant * q1); // av prox between affected pixels in one quadrant and another

	void getNearestNeighbors(quadrant * q, quadrant * tQ); // one quadrant, every other quadrant

	ham * tHam(quadrant * q);     // hamming distances between all nodes in every quadrant 

	void iHam(quadrant * q);


// data shared by various clusters

struct shared(int){

int * i; // heart cluster index
int * nNodes; // number of nodes that are shared with other clusters
int * sharedN; //indices of shared nodes
int * sharedC; // indices of clusters, starting with the first shared node in the above array, with whom there are shared nodes 
long * sharedE;   // respective hamming distances between the shared nodes edges and the edges of this cluster 
					// (so, how much edge symmetry is there in the radial distances in edgewise clusters
					// first array is the nodal index (or indices), second is hamming distance of the edges that are shared

}; // shared


// hamming functions

// calculates strings within the image and then calcuates their hamming distance(s), looking for similarities within classification strategy
void hamming(bool * pruned, char * imageData); // has been pruned for uninteresting data (see silver and bright light rgb for references)




// data for python sdk deep learning



//////////////////// BEGIN NEURAL DATA //////////////////////////////


struct lobe(int * _n)
{
	neuron * n = malloc(_n, sizeof(neuron));

}

struct neuron(int * c)
{

	// hamming distances between nodes in quadrants of main image

	int * t1 hDistanceThreshTL;
	int * t2 hDistanceThreshTL;
	int * t3 hDistanceThreshTL;

	int * t1 hDistanceThreshTR;
	int * t2 hDistanceThreshTR;
	int * t3 hDistanceThreshTR;

	int * t1 hDistanceThreshBL;
	int * t2 hDistanceThreshBL;
	int * t3 hDistanceThreshBL;

	int * t1 hDistanceThreshBR;
	int * t2 hDistanceThreshBR;
	int * t3 hDistanceThreshBR;


	// cervical coordinate thresholds

	short * cervixYThresh; // where does the rgb data near the middle of the image change to reflect a cervical modulation in the endometrial tissue
	short * cervixCrThresh;
	short * cervixCbThresh;

	long * cervixPinkThresh;
	long * cervixPurThresh;
	long * cervixSkyThresh;


	// coordinates of x and y within image where cervical rgb data matches the above thresholds
	// is class independent because it reflects merely the location, which should be near the center

	int * xThreshY;
	int * yThreshY;

	int * xThreshCr;
	int * yThreshCr;

	int * xThreshCb;
	int * yThreshCb;

	int * xThreshPink;
	int * yThreshPink;

	int * xThreshPur;
	int * yThreshPur;

	int * xThreshSky;
	int * yThreshSky;



	// endoscopic thresholds

	short * t1EndoCrThresh;
	short * t2EndoCrThresh;
	short * t3EndoCrThresh;
	
	short * t1EndoCbThresh;
	short * t2EndoCbThresh;
	short * t3EndoCbThresh;
	
	short * t1EndoYThresh;
	short * t2EndoYThresh;
	short * t3EndoYThresh;
	
	

	// modulated endoscopic thresholds
	
	long * t1EndoPinkThresh;
	long * t1EndoPurThresh;
	long * t1EndoSkyThresh;
	
	long * t2EndoPinkThresh;
	long * t2EndoPurThresh;
	long * t2EndoSkyThresh;
	
	long * t3EndoPinkThresh;
	long * t3EndoPurThresh;
	long * t3EndoSkyThresh;
	
	
	// ectoscopic thresholds
	
	short * t1EctoCrThresh;
	short * t2EctoCrThresh;
	short * t3EctoCrThresh;
	
	short * t1EctoCbThresh;
	short * t2EctoCbThresh;
	short * t3EctoCbThresh;
	
	short * t1EctoYThresh;
	short * t2EctoYThresh;
	short * t3EctoYThresh;
	
	
	// modulated ectoscopic thresholds
	
	long * t1EctoPinkThresh;
	long * t1EctoPurThresh;
	long * t1EctoSkyThresh;
	
	long * t2EctoPinkThresh;
	long * t2EctoPurThresh;
	long * t2EctoSkyThresh;
	
	long * t3EctoPinkThresh;
	long * t3EctoPurThresh;
	long * t3EctoSkyThresh;	


	// total edges and vertices between nodes in the summed quadrant space

	int * nDendrites; // # of quadrants to which this dendrite is connected
	int * nAxons;

	dendrite * d;
	axon * a;


	// simple boolean data on regional locality of affected nodes and their respective quadrants

	bool * isEcto; // is there a significant amount of affected px definitively in the ectoscopic area 
	bool * isEndo; // is there a significant amount of the affected px definitively in the endoscopic area


	// in what class is the lobe

	bool * c1 = 0;
	bool * c2 = 0;
	bool * c3 = 0;

	if(c == 1) c1 = true;
	if(c == 2) c2 = true;
	if(c == 3) c3 = true;

}; // neuron