//#include<math.h>
//#include<iostream>
//#include<clock.h> 
#include<stdlib.h>
//#include<system.h>

// check to see if the system command is listed as above


virtual cure::cure(){};
virtual cure::~cure(){};

const double * cure::getLen(x1, x2, y1, y2)
{
	const double * len = malloc(1, sizeof(double));
	const len = sqrt(fabs(x2-x1)*fabs(x2-x1) + fabs(y2-1y1)*fabs(y2-y1));
	return len;

} // getLen
 

const px * cure::pixel(char * c, , const int * x, const int * y, const int * i)
{
		const px * p = malloc(1, sizeof(px));
		p = iPx(c, x, y, i);
		return p;
} // px



const lobe * cure::lobe(quadrant * q, int * nq, int * class)
{

	int * nNeurons = malloc( rand() % ( nQ*nQ ) ), sizeof(int) );

	lobe * b = malloc(nNeurons, sizeof(lobe));

	b = iLobe(q, nq);

	int * i = malloc(1, sizeof(int));

	for(int i = 0; i < nNeurons; i++)
	{
		b[i] = aggregate(q, class);

	} // for

} // lobe



// TODO: parallelize getData so that each thread runs a different quadrant
/*
void cure::getData(quadrant * q)
{
	int * temp = malloc(1, sizeof(int));
	int * rX = malloc(1, sizeof(int));
	int * cX = malloc(1, sizeof(int));


	long * confidence = malloc(1, sizeof(long) );
	
	q->confidence = rand();

	for(temp = 0; temp < nq; temp++)
	{
		for(rX =  q[index]->startX; rX < q[index]->enn; rX++)
			{
				q[index]->confidence = confidence;

					for(cX = q[index]->startY; cX < q[index]->endY; cX++)
						{
							if( present(q[index]->px[rX][cX]->data) >= q[index]->confidence ) 
							{
								q[index]->px[rX][cX]->affected = true;
								q[index]->affected++;
							} // if
							else 
							{
								 q[index]->px[rX][cX]->affected = false;
								 q[index]->unaffected++;
							} // else

					} // for all columns in quadrant

			} // for all rows in quadrant

	} // for all quadrants
	

} // getData

*/

// largest data structure by far, must be hyperthreaded 
/*
void cure::findEdges(quadrant * q){

	int * y - malloc(1, sizeof(int));
	y = (q->nNodes*q->nNodes);

	edge * e = malloc(y, sizeof(edge));
	
	int * i = malloc(1, sizeof(int));
	int * j = malloc(1, sizeof(int));

for(j = 1; j < q->nNodes; j++)
{	
	for(i = 0; i < q->nNodes-1; i++)
	{	
		e[y] = findEdge(q->nodes[i]->x, q->nodes[i]->y, q->nodes[j]->x, q->nodes[j]->y);	
		y++;
	} // for i
	
	i = 0;

}// for j

	i = 0;
	j = 0;

	while( ( i <  q->nNodes ) )
	{
		while( ( j <  (q->nNodes-1) ) )
		{
			computeHams(e[i*j], e[i*j+1]);
			i++;
			j++;
		} // while

	} // while

} // find edges

void cure::findEdge(int * x1, int * x2, int*  y1,int * y2){

	double * len = malloc(1, sizeof(int));

	len = getLen(x1, x2, y1, y2);

	float * rem = malloc(1, sizeof(float));
	int * index = malloc(1, sizeof(int));

		if(len % .000000001 == 0){

 		px * edge = malloc((int)len, sizeof(px));
 
 			while( ( x1 < x2) && (y1 < y2) && (index < len) )
			{
				edge[index] = pxls[x1++][y1++];
				index++;
			} // while
		} // if
		else
		{
	rem = (len % .000000001);
    if(rem < .5){

		while( ( x1 < x2) && (y1 < y2) && (index < len) )
			{
				if(index = (len-2) ) edge[index] = pxls[x1++][y1+2];
				else{ edge[index] = pxls[x1++][y1++]; } // else
				index++;
 			} // while

		
	} // inner if
	else{

		while( ( x1 < x2) && (y1 < y2) && (index < len) )
			{
				if(index = (len-2) ) edge[index] = pxls[x1+2][y1++];
				else{ edge[index] = pxls[x1++][y1++]; } // else
				index++;
 			} // while
		} // inner else

	} // else


}; // edge struct



ham * cure::getHam(int * len, edge * one, edge * two)
{
	ham * h = malloc(1, sizeof(ham));

	while(len > 0)
	{
	h[0][index] = fabs(nl(one[index]->Y - two[index]->Y));
	h[1][index] = fabs(nl(one[index]->Cr - two[index]->Cr));
	h[2][index] = fabs(nl(one[index]->Cb - two[index]->Cb));

	h[3][index] = fabs(nl(one[index]->pink * two[index]->pink));
	h[4][index] = fabs(nl(one[index]->pur * two[index]->pur));
	h[5][index] = fabs(nl(one[index]->sky * two[index]->sky));
	} // while
	
	return h;
} // hamming string struct 


void cure::iHam(int * _len){

	int * len = malloc(_len, sizeof(int));
	len = _len;

	short * dY = malloc(_len, sizeof(short));
	short * dCr = malloc(_len, sizeof(short));
	short * dCb = malloc(_len, sizeof(short));

	long * dPink = malloc(_len, sizeof(long));
	long * dPur = malloc(_len, sizeof(long));
	long * dSky = malloc(_len, sizeof(long));

	long ** h = malloc(6, sizeof(malloc(len, sizeof(long));

}; // hamming string struct



int * cure::getHams(edge * one, edge * two){

	int * index = malloc(1, sizeof(int));
	index = 0; 
	if(one->edge->len > two->edge->len) int * ham = malloc(two->edge->len, sizeof(int));
	else {
		int * ham = malloc(one->edge->len, sizeof(int));
	}
	while( ( one->edge->len > index) && (two->edge->len > index) ) {

		ham[index] = one->edge[index] - two->edge[index];
		index++;

	} // while

	return ham;

} // ham


// send global data to file
ham * cure::computeHamming(quadrant * q1, quadrant * q)
{
	int * i = malloc(1, sizeof(int));
	int * j = malloc(1, sizeof(int));
	int * k = malloc(1, sizeof(int));

	computeEdges(q1);
	ham ** h = malloc((nq - 1), sizeof(, malloc(int));
			
			for(j = 0; j < (this->nq-1); j++)
			{
					computeEdges(q[j]);

				for(i = 0; i < q1->nEdges; i++)
					{

					for(j = 0; j < (this->nq-1); j++)
					{
						for(k = 0; k < q[j]->nEdges; k++)
						{
								getHams()


						}

					}
			} // for

	}

		return h;

}// computeInterHamming

ham * cure::computeIntraHamming(quadrant * q)
{
	int * i = malloc(1, sizeof(int));
	int * v = malloc(1, sizeof(int));
	node * temp = malloc(1, sizeof(node));


	for(v = 0; v < q->nNodes; v++)
	{
			temp = q->nodes[v];

		for(i = 0 ; i < (q->nNodes * q->nNodes); i++){

			q->intra[i] = getHams();		

		}
	}

}


} // present



void cure::isEcto(quadrant * q)
{

	int * index = malloc(1, sizeof(int));
	
	for(index = 0; index < nq; index++)
	{
		if( (q[index]->enn > cpx[cSz]->x) && (q[index]->startY < cpx[cSz]->y) )
		{	
			q[index]->isEcto == true;
		}
	}
} // isEcto

void cure::isEndo(quadrant * q)
{
	int * index = malloc(1, sizeof(int));
	
	for(index = 0; index < nq; index++)
	{
		if(q[index]->enn < cpx[cSz]->x) && (q[index]->startY > cpx[cSz]->y)
		{	
			q[index]->isEndo == true;
		}

} // isEndo

}

//does this fall between ecto and endo params

bool cure::middle(quadrant * q, long * tolerance)
{
	bool indeterminate = false;
	return indeterminate;
} // middle
*/

px * cure::iPx(char * b, const int * _x, const int * _y, const * int _i)
{
		px * p = malloc(1, sizeof(px));

		int * x = malloc(1, sizeof(int) );
		p->x = _x;

		int * y  = malloc(1, sizeof(int) );
		p->y = _y;

		int * i  = malloc(1, sizeof(int) );
		p->ii = _i; // image index
		
		int * iN  = malloc(1, sizeof(int) ); // index of node for which this color data is connected
		p->nIndx = iN;

		int ** qi; // index of quadrant of the node

		int ** e; // edge indices of hamming strings to which this pixel is connected
		
		p->Y = malloc(1, sizeof(short) ); // the luminesence (brightness) of the node
		p->Cr = malloc(1, sizeof(short) ); // the reddish chromesence (interesting per Hubble's theory deduced when observing red shift galactic movement)
		p->Cb = malloc(1, sizeof(short) ); // the blueish chromesence
		
		p->pur = malloc(1, sizeof(short) ); //red - blue
		p->pink = malloc(1, sizeof(short) ); //brightness - red
		p->sky = malloc(1, sizeof(short) ); //brightness - blue

		p->Y = ( ( (int)b >> 16 ) & 0xFF ) / 255.0; //  luminesence
		p->Cr = ( ( (int)b >> 8 ) & 0xFF ) / 255.0;  //  red shift 
	    p->Cb =  ( (int)b & 0xFF ) / 255.0;       //  blue shift

		p->pur = fabs( ln(Cr * Cb) );
		p->pink = fabs( ln(Y * Cr) );
		p->sky = fabs( ln(Y * Cb) );

} // pixel data



pair * cure::factors(int * n)
{
	int * t1 = malloc(1, sizeof(int));
	int * t2 = malloc(1, sizeof(int));

	const int * x = fabs(sqrt(n));

	// n's' factor magnitude is no greater than the square root of 
	// its magnitude minus 1, computationally the subtraction is trivial

	pair * p = malloc(x, sizeof(pair));
	
	int * g = malloc(1, sizeof(int));
	int * o = malloc(1, sizeof(int));
	
	o = 0;

	int * index = malloc(1, sizeof(int));
	
	index = 0;
	
	this->nFactors = malloc(1, sizeof(int));

	g = n;
	t = 1;

	while( x != 0 )
	{
				t++;
				g--;
				x--;

			if( ((n % t) == 0) && ((n % g) == 0)) 
			{
				this->nFactors += 2;
				pair[o++] = iPair(index++, t, g);
			} // if
				
	} // while  
	
} // factors

//init rows and columns (arbitrarily taken from middle two factors
void cure::iRC()
{
	
	this->r = this->factors[(int)(this->nFactors/2)]->f1;
	this->c = this->factors[(int)(this->nFactors/2)]->f2;

} // iRC


///////////////// INITIALIZE QUADRANTS  /////////////////////////


void cure::quadrant(quadrant * q, int * i, int * x, int * y)
{

	q->i = malloc(1, sizeof(int));
	q->i = i;

	q->cols = malloc(1, sizeof(int));
	q->rows = malloc(1, sizeof(int));

	q->Y = malloc(1, sizeof(short) );
	q->Cr = malloc(1, sizeof(short) );
	q->Cb = malloc(1, sizeof(short) );
	
	q->pur = malloc(1, sizeof(long) );
	q->pink = malloc(1, sizeof(long) );
	q->sky = malloc(1, sizeof(long) );	

	q->minY = malloc(1, sizeof(short) );
	q->minCr = malloc(1, sizeof(short) );
	q->minCb = malloc(1, sizeof(short) );

	q->minPur = malloc(1, sizeof(long) );
	q->minPink = malloc(1, sizeof(long) );
	q->minSky = malloc(1, sizeof(long) );	

	q->maxY = malloc(1, sizeof(short) );
	q->maxCr = malloc(1, sizeof(short) );
	q->maxCb = malloc(1, sizeof(short) );
	
	q->maxPur = malloc(1, sizeof(long) );
	q->maxPink = malloc(1, sizeof(long) );
	q->maxSky = malloc(1, sizeof(long) );

	q->avY = malloc(1, sizeof(short) );
	q->avCr = malloc(1, sizeof(short) );
	q->avCb = malloc(1, sizeof(short) );
	
	q->avPur = malloc(1, sizeof(long) );
	q->avPink = malloc(1, sizeof(long) );
	q->avSky = malloc(1, sizeof(long) );

	q->bRight = malloc(1, sizeof(bool));
	q->tTight = malloc(1, sizeof(bool));

	q->bLeft = malloc(1, sizeof(bool));
	q->tLeft = malloc(1, sizeof(bool));

	q->cy1 = malloc( 1, sizeof(int) ); // in ascending order 
	q->cx1 = malloc( 1, sizeof(int) ); // the coordinates that contain pixels imaging the cervix proper

	q->cy2 = malloc( 1, sizeof(int) ); // these are the end points of the cervix, so x coordinates go from cx1 to cx2, y from cy1 to cy2
	q->cx2 = malloc( 1, sizeof(int) );

	q->ne = malloc( 1, sizeof(int) );
	
	q->cols = x;
	q->rows = y;

	int * e = malloc(2, sizeof(int)); // x and y distance from edge of image

	// TODO: calculate size of cervix and add buffer to this number to acertain endoscopic or ectoscopic position of pixel data

	e[0] = (this->r/2) - q->rows; // if positive, then in lower half of image, else in upper half
	e[1] = (this->c/2) - q->cols; // if positive, then in left half of image, else in right half


	//TODO: check parameters to ensure these placements are right

	if( (e[0] > 0) && (e[1] > 0) ) q->bRight = true; 
	else if( (e[0] < 0) && (e[1] < 0) ) q->tRight = true;

	else if( (e[0] > 0) && (e[1] < 0) ) q->bLeft = true;
	else if( (e[0] < 0) && (e[1] > 0) ) q->tLeft = true;	
	
	if( (e[0] == (this->r/2)) && (e[1] == (this->c/2)) ) isBRightVis = false;
	else if( (e[0] == (-1)*(this->r/2)) && (e[1] == (-1)*(this->c/2)) ) isTRightVis = false;

	else if( (e[0] == (this->r/2)) && (e[1] == (-1)*(this->c/2)) ) isBLeftVis = false;
	else if( (e[0] == (-1)*(this->r/2)) && (e[1] == (this->c/2)) ) isTLeftVis = false;

	q->area = q->rows * q->cols;

	q->nNodes = ( rand() % ( (int)( pow(q->area, (1/2)) ) ) ) + 1; // square root of area plus one is limit of nearest neighbor 
	
	ne = (int)pow(q->nNodes, (1/2) );	

	int tempX = malloc(1, sizeof(int));
	int tempy = malloc(1, sizeof(int));

	ham * inter = malloc(q->nNodes * q->nNodes, sizeof(ham) ); // hamming distances between center of quadrant and its respective nodes
	ham * intra; // hamming distances between each nodal center and nodal centers of each linearly adjacent quadrant

	px * nodes = malloc(q->nNodes, sizeof(px) ); // central pixel + nearest neighbors
	
	q->nodes = nodes;

} // iQuadrant



// assign values to each variable in the quadrant

quadrant * cure::iQ(int * i, int * x, int * y)
{
	quadrant * q = malloc(1, sizeof(quadrant));

	// x and y are boudary conditions derived from natural number factors of the original image
	// no natural number larger than ~half of the image is used as either x or y, hence the name: quadrant

	quadrant(q, i, x, y); 

	int * n = malloc( 1, sizeof(int) );
	n = 0;

	while(n < q->nNodes)
	{
			// choose nodes from within quadrant
			tempY = ( rand() % q->rows) + q->rows) ); 
			tempX = ( rand() % q->cols) + q->cols) ); 
			
			// last two variables are for island sharing
			q->node[n] = iNode(pxls[tempX][tempY], q->area);
			
			if( q->node[n]->minY < q->minY ) q->minY = q->node[n]->minY;
			if( q->node[n]->maxY > q->maxY ) q->maxY = q->node[n]->maxY;

			if( q->node[n]->minCr < q->minCr ) q->minCr = q->node[n]->minCr;
			if( q->node[n]->maxCr > q->maxCr ) q->maxCr = q->node[n]->maxCr;

			if( q->node[n]->minCb < q->minCb ) q->minCb = q->node[n]->minCb;
			if( q->node[n]->maxCb > q->maxCb ) q->maxCb = q->node[n]->maxCb;

			if( q->node[n]->minPink < q->minPink ) q->minPink = q->node[n]->minPink;
			if( q->node[n]->maxPink > q->maxPink ) q->maxPink = q->node[n]->maxPink;

			if( q->node[n]->minPur < q->minPur ) q->minPur = q->node[n]->minPur;
			if( q->node[n]->maxPur > q->maxPur ) q->maxPur = q->node[n]->maxPur;

			if( q->node[n]->minSky < q->minSky ) q->minSky = q->minSky;
			if( q->node[n]->maxSky > q->maxSky ) q->maxSky = q->maxSky;


			q->avY += q->node[n]->avY;
			q->avCr += q->node[n]->avCr;
			q->avCb += q->node[n]->avCb;

			q->avPink += q->node[n]->avPink;
			q->avPur += q->node[n]->avPur;
			q->avSky += q->node[n]->avSky;
		
			n++;

		} // while

			q->avY /= q->nNodes;
			q->avCr /= q->nNodes;
			q->avCb /= q->nNodes;
			
			q->avPink /= q->nNodes;
			q->avPur /= q->nNodes;
			q->avSky /= q->nNodes;

			int * index = malloc(1, sizeof(int));

			for(index = 0; index < q->nNodes; index++){

				q->sdCr += (q->avCr - q->ne[index]->Cr) * (q->avCr - q->ne[index]->Cr) ;
				q->sdCb +=  (q->avCb - q->ne[index]->Cb) * (q->avCb - q->ne[index]->Cb) ;
				q->sdY +=  (q->avY - q->ne[index]->Y) * (q->avY - q->ne[index]->Y) ;

				q->sdPink +=  (q->avPink - q->ne[index]->Pink) * (q->avPink - q->ne[index]->pink) ;
				q->sdPur +=  (q->avPur - q->ne[index]->Pur) * (q->avPur - q->ne[index]->Pur) ;
				q->sdSky +=  (q->avSky - q->ne[index]->Sky) * (q->avSky - q->ne[index]->Sky) ;
	
				} // for all nearest neighbor variances


				q->sdCr /= q->nNodes;
				q->sdCb /= q->nNodes;
				q->sdY /= q->nNodes;

				q->sdPink /= q->nNodes;
				q->sdPur /= q->nNodes;
				q->sdSky /= q->nNodes;


				// root mean square

				q->sdCr = sqrt(q->sdCr); 
				q->sdCb = sqrt(q->sdCb); 
				q->sdY = sqrt(q->sdY); 
				q->sdPur = sqrt(q->sdPur); 
				q->sdPink = sqrt(q->sdPink); 
				q->sdSky = sqrt(q->sdSky); 

				// second,     third and fourth deviations away from Gaussian center

				q->sd2Cr = q->sdCr + q->avCr * 2; 
				q->sd2Cb = q->sdCb + q->avCr * 2;
				q->sd2Y = q->sdY + q->avCr * 2;
				q->sd2Pur = q->sdPur + q->avCr * 2; 
				q->sd2Pink = q->sdPink + q->avCr * 2;
				q->sd2Sky = q->sdSky + q->avCr * 2;

				q->sd3Cr = q->sdCr + q->avCr * 3;
				q->sd3Cb = q->sdCb + q->avCb * 3; 
				q->s3Y = q->sdY + q->avY * 3;
				q->sd3Pur = q->sdPur + q->avPur * 3;
				q->sd3Pink = q->sdPink + q->avPink * 3;
				q->sd3Sky = q->sdSky + q->avSky * 3;

				q->sd4Cr = q->sdCr + q->avCr * 4; 
				q->sd4Cb = q->sdCb + q->avCb * 4;
				q->sd4Y = q->sdY + q->avY * 4;
				q->sd4Pur = q->sdPur + q->avPur * 4; 
				q->sd4Pink = q->sdPink + q->avPink * 4; 
				q->sd4Sky = q->sdSky + q->avSky * 4; 

} // iQ



node * cure::iNode(px * center)
{
		node * n = malloc(1, sizeof(node) );
		int * index = malloc(1, sizeof(int));

		n->center = center;

		n->ne = malloc(8, sizeof(px)); 
										
		short * Y = malloc(8, sizeof(short) ); // the luminesence (brightness) of the node
		short * Cr = malloc(8, sizeof(short) ); // the reddish chromesence (interesting per Hubble's theory deduced when observing red shift galactic movement)
		short * Cb = malloc(8, sizeof(short) ); // the blueish chromesence
		
		long * pur = malloc(8, sizeof(short) ); //red - blue
		long * pink = malloc(8, sizeof(short) ); //brightness - red
		long * sky = malloc(8, sizeof(short) ); //brightness - blue

		// top left neighbor

		if(NULL != this->pxls[center->x-1][center->y+1]) n->ne[0] = this->pxls[center->x-1][center->y+1];
		if(NULL != this->pxls[center->x][center->y+1])   n->ne[1] = this->pxls[center->x][center->y+1];
		if(NULL != this->pxls[center->x+1][center->y+1]) n->ne[2] = this->pxls[center->x+1][center->y+1];
		if(NULL != this->pxls[center->x-1][center->y])   n->ne[3] = this->pxls[center->x-1][center->y];
		if(NULL != this->pxls[center->x+1][center->y])   n->ne[4] = this->pxls[center->x+1][center->y];
		if(NULL != this->pxls[center->x-1][center->y-1]) n->ne[5] = this->pxls[center->x-1][center->y-1];
		if(NULL != this->pxls[center->x][center->y-1])   n->ne[6] = this->pxls[center->x][center->y-1];
		if(NULL != this->pxls[center->x+1][center->y-1]) n->ne[7] = this->pxls[center->x+1][center->y-1];

		for(index = 0; index < 8; index++){

				n->avCr += n->ne[index]->Cr;
				n->avCb += n->ne[index]->Cb;
				n->avY += n->ne[index]->Y;

				n->avPink += n->ne[index]->pink;
				n->avPur += n->ne[index]->pur;
				n->avSky += n->ne[index]->sky;

				
				if(n->minY > n->ne[index]->Y) n->minY = n->ne[index]->Y;
				if(n->minCr > n->ne[index]->Cr) n->minCr = n->ne[index]->Cr;
				if(n->minCb > n->ne[index]->Cb) n->minCb = n->ne[index]->Cb;
			
				if(n->minPink > n->ne[index]->Pink) n->minPink = n->ne[index]->Pink;
				if(n->minPur > n->ne[index]->Pur) n->minPur = n->ne[index]->Pur;
				if(n->minSky > n->ne[index]->Sky) n->minSky = n->ne[index]->Sky;
			
				if(n->maxY < n->ne[index]->Y) n->maxY = n->ne[index]->Y;
				if(n->maxCr < n->ne[index]->Cr) n->maxCr = n->ne[index]->Cr;
				if(n->maxCb < n->ne[index]->Cb) n->maxCb = n->ne[index]->Cb;
			
				if(n->maxPink < n->ne[index]->Pink) n->maxPink = n->ne[index]->Pink;
				if(n->maxPur < n->ne[index]->Pur) n->maxPur = n->ne[index]->Pur;
				if(n->maxSky < n->ne[index]->Sky) n->maxSky = n->ne[index]->Sky;

	
		} // for all nearest neighbor max/min 

				// average values

				n->avCr /= 8;
				n->avCb /= 8;
				n->avY /= 8;

				n->avPink /= 8;
				n->avPur /= 8;
				n->avSky /= 8;

			// variances

				for(index = 0; index < 8; index++){

				n->sdCr += (n->avCr - n->ne[index]->Cr) * (n->avCr - n->ne[index]->Cr) ;
				n->sdCb +=  (n->avCb - n->ne[index]->Cb) * (n->avCb - n->ne[index]->Cb) ;
				n->sdY +=  (n->avY - n->ne[index]->Y) * (n->avY - n->ne[index]->Y) ;

				n->sdPink +=  (n->avPink - n->ne[index]->Pink) * (n->avPink - n->ne[index]->pink) ;
				n->sdPur +=  (n->avPur - n->ne[index]->Pur) * (n->avPur - n->ne[index]->Pur) ;
				n->sdSky +=  (n->avSky - n->ne[index]->Sky) * (n->avSky - n->ne[index]->Sky) ;
	
		} // for all nearest neighbor variances

				n->sdCr /= 8;
				n->sdCb /= 8;
				n->sdY /= 8;

				n->sdPink /= 8;
				n->sdPur /= 8;
				n->sdSky /= 8;


				// root mean square

				n->sdCr = sqrt(n->sdCr); 
				n->sdCb = sqrt(n->sdCb); 
				n->sdY = sqrt(n->sdY); 
				n->sdPur = sqrt(n->sdPur); 
				n->sdPink = sqrt(n->sdPink); 
				n->sdSky = sqrt(n->sdSky); 

				// second, third and fourth deviations away from Gaussian center

				n->sd2Cr = n->sdCr + n->avCr * 2; 
				n->sd2Cb = n->sdCb + n->avCr * 2;
				n->sd2Y = n->sdY + n->avCr * 2;
				n->sd2Pur = n->sdPur + n->avCr * 2; 
				n->sd2Pink = n->sdPink + n->avCr * 2;
				n->sd2Sky = n->sdSky + n->avCr * 2;

				n->sd3Cr = n->sdCr + n->avCr * 3;
				n->sd3Cb = n->sdCb + n->avCb * 3; 
				n->s3Y = n->sdY + n->avY * 3;
				n->sd3Pur = n->sdPur + n->avPur * 3;
				n->sd3Pink = n->sdPink + n->avPink * 3;
				n->sd3Sky = n->sdSky + n->avSky * 3;

				n->sd4Cr = n->sdCr + n->avCr * 4; 
				n->sd4Cb = n->sdCb + n->avCb * 4;
				n->sd4Y = n->sdY + n->avY * 4;
				n->sd4Pur = n->sdPur + n->avPur * 4; 
				n->sd4Pink = n->sdPink + n->avPink * 4; 
				n->sd4Sky = n->sdSky + n->avSky * 4; 


} // iNode



/*
double * cure::intelligence()
{

	// how to classify?
	
	/*
	 
	write node data, radii from center, and number of nodes for each cluster

	write rules that entail aggregating the most common results for k-nearest neighbors regarding nodal data in local clusters
	
	exchange nodes between clusters (island sharing) for accuracy improvements, do this with messgae passing

	mutate buffer, radius, nNodes, nClusters 

	fitness tested by least quadrants, nodes, and eventually hamming distances between vertices and edges
	necessary to classify, thresholds of rgb data, least change in variances, etc

	smooth noise with neural net sdk, and prune attributes
	
	
	*/
	
} // fit
 
*/

 // TODO: aggregate hamming distances * lnfabs(v1-v2) for each edge pair of vertices
 // 

/////////// MOVE DATA TO PROLOG //////////////////////

// classify in prolog after hamming strings and pruning is done in C 
/*
void cure::log()
{
	
	File * f = fopen(this->o, 'w');
	
	this->begin.stop();
	
	fprint(%c, %d, 'runtime, ' this->t);
	fprint(%c, %i, 'pixel count', this->npx);
	
	fprint(%c, %d, 'fitness, ' this->fitness);
	
	fprint(%c, %i, 'cluster lim, ', this->C);	
	fprint(%c, %i, 'radial lim, ', this->R);
	
	fprint(%c, %d, 'aveRadius, ', this->avRadius);
	fprint(%c, %i, 'node lim, ', this->node);
	
	fprint(%c, %i, 'tMargin, ', margin[0]);  
	fprint(%c, %i, 'rMargin, ', margin[1]); 
	
	
} // log	
 
 
*/


// organize image data from original jpg


char * cure::bits()
{
		
File * m = fopen(this->i, 'r');

fseek( m, 0L, SEEK_END );

this->n = malloc(1, sizeof(int));
this->n = static_cast<int>( ftell(m) );

this->byt = malloc( this->n, 1 );

this->pxls = malloc(this->n, sizeof(px) );

int e = malloc(1, sizeof(int) );
	
this->factors = factors(this->n);

iRC(this->factors);

	for(e = 0; e < this->n; e++)
		{
			fread( this->byt[e], 1, 1, m );
			this->pxls[e] = pixel( this->byt[e], (e % this->r), (e % this->c), e );
		} // for
	
	fclose( m );
 
}  // bits	


///////////////// HAMMING DISTANCE CALCULATIONS /////////////////////////////

/*
char** cure::hamString(px * one, px * two){
 
		int len = malloc(1, sizeof(int) );

		len = (int)sqrt( abs( one->x1 - one->x2 ) * abs( two->x1 - two->x2 ) + abs( one->y1 - one->y2 ) * abs( two->y1 - two->y2 ) ) ;
		
		char * ham = malloc(len, 1);  
		
		int dn;
		
		int n = x1;
		int dy = y1;
		int x = 0;
		int y = 0;

			for(dn = 0; dn < len ; dn++){

				string[x][y] = byt[n + dn][dy + dn];		
				x++;
				y++;
			
			} // for string length

	}// linear hamming string

char** cure::cHamstring(int * x1, int * x2, int * y1, int * y2, int * curveX, int * curveY){

		int len = malloc(1, sizeof(int) );

		len = sqrt(abs(x1-x2)*abs(x1-x2) + abs(y1-y2)*abs(y1-y2) ) ;
		
		char * ham = malloc(len, 1); 
		
		int dn = 0;

		int x = 0;
		int y = 0;
		
		int newX = curveX;
		int newY = curveY;

			while(dn < len) && ( (newX *= dn) < this->rows) ) && (newY *= dn) <= this->cols) ) ){

				newX *= dn;
				newY *= dn;

				string[x][y] = byt[x1 + newX][y1 + newY];
				
				x++;
				y++;
			

			} // for string length

			return ham;

	} // curved hamming string

//write to log file
}

*/


void cure::log()
{
	// cluster metrics
	char * radii = malloc(3, 1); // min/max/av
	
	char * ys = malloc(3, 1); // min/max/av
	char * crs = malloc(3, 1); // min/max/av
	char * cbs = malloc(3, 1); // min/max/av
	
	char * pinks = malloc(3, 1); // same as above
	char * purs = malloc(3, 1);
	char * skys = malloc(3, 1);
	
	char * nodes = malloc(3, 1);
	
	char ** data[8][3];

		
	File * f = fopen(this->pro, 'w');

	int i;
	int n; 
		
	for(i = 0; i < this->nNeurons; i++ )
	{	
			radii[0] = c[i]->minRadius;
			radii[1] = c[i]->maxRadius;
			radii[2] = c[i]->avRadius;
			
			crs[0] = c[i]->minCr;
			crs[1] = c[i]->maxCr;
			crs[2] = c[i]->avCr;
			
			cbs[0] = c[i]->minCb;
			cbs[1] = c[i]->maxCb;
			cbs[2] = c[i]->avCb;
			
			ys[0] = c[i]->minY;
			ys[1] = c[i]->maxY;
			ys[2] = c[i]->avY;
			
			pinks[0] = c[i]->minPink;
			pinks[1] = c[i]->maxPink;
			pinks[2] = c[i]->avPink;
			
			purs[0] = c[i]->minPur;
			purs[1] = c[i]->maxPur;
			purs[2] = c[i]->avPur;
			
			skys[0] = c[i]->minSky;
			skys[1] = c[i]->maxSky;
			skys[2] = c[i]->avSky;
			
			nodes[0] = c[i]->center->y;
			nodes[0] = c[i]->center->x;
			
			nodes[1] = c[i]->nNodes;
			nodes[2] = c[i]->class;
		
					
			// aggregate the data
				
			data[0] = radii;
			data[1] = ys;
			data[2] = crs;
			data[3] = cbs;
			data[4] = pinks;
			data[5] = purs;
			data[6] = skys;
			data[7] = nodes;


/*			
			//lobal data, quadrants' index, nNodes, list of radic data, lists of chromic data, list 
			fprint(%c, %i, %c, %c, %i, 'lobe: ', i, ', class: ', data[7][2], '\n');
			
			fprint(%c, %i, %c, %i, %c, 'quadrant(', i, ,',', data'):-\n');
			
			fprint(%c, %i, %c, %i, %c, 'rowLen(', i, ,',', data'):-', ',');
			fprint(%c, %i, %c, %i, %c, 'colLen(', i, ,',', data'):-', ',');

			fprint(%c, %i, %c, %i, %c, 'nNodes(', i, ',', data[7][1], ',', data[7][2], '),\n'); 
			
			fprint(%c, %i, %c, %i, %c, 'centerIndex(', i, ',', data[7][0], ',', data[7][2], '),\n'); 
			
			fprint(%c, %i, %c, %i, %c, 'maxRadius(', i, ',', data[0][1], ',', data[7][2], '),\n'); 
			fprint(%c, %i, %c, %i, %c, 'minRadius(', i, ',', data[0][0], ',', data[7][2], '),\n'); 
			
			fprint(%c, %i, %c, %i, %c, 'maxY(', i, ',', data[1][1], ',',data[7][2], '),\n'); 
			fprint(%c, %i, %c, %i, %c, 'minY(', i, ',', data[1][0], ',', data[7][2], '),\n'); 
			
			fprint(%c, %i, %c, %i, %c, 'maxCr(', i, ',', data[2][1], ',', data[7][2], '),\n'); 
			fprint(%c, %i, %c, %i, %c, 'minCr(', i, ',',  data[2][0], ',',data[7][2], '),\n'); 
			
			fprint(%c, %i, %c, %i, %c, 'maxCb(', i, ',',  data[3][1], ',', data[7][2], '),\n'); 
			fprint(%c, %i, %c, %i, %c, 'minCb(', i, ',',  data[3][0], ',', data[7][2], '),\n'); 
			
			fprint(%c, %i, %c, %i, %c, 'maxPink(', i, ',', data[4][1], ',', data[7][2], '),\n'); 
			fprint(%c, %i, %c, %i, %c, 'minPink(', i, ',', data[4][0], ',', data[7][2], '),\n'); 
			
			fprint(%c, %i, %c, %i, %c, 'maxPur(', i, ',',  data[5][1], ',', data[7][2], '),\n'); 
			fprint(%c, %i, %c, %i, %c, 'minPur(', i,',',  data[5][0], ',', data[7][2], '),\n'); 
			
			fprint(%c, %i, %c, %i, %c, 'maxSky(', i, ',',  data[6][1], ',', data[7][2], '),\n'); 
			fprint(%c, %i, %c, %i, %c, 'minSky(', i, ',', data[6][0], ',', data[7][2], '),\n'); 
			
			fprint(%c, %i, %c, %i, %c, 'avRadius(', i, ',',  data[0][2], ',', data[7][2], '),\n'); 

			fprint(%c, %i, %c, %i, %c, 'avY(', i, ',',  data[1][2], ',', data[7][2], '),\n'); 
			fprint(%c, %i, %c, %i, %c, 'avCr(', i, ',', data[2][2], ',', data[7][2], '),\n'); 
			fprint(%c, %i, %c, %i, %c, 'avCb(', i, ',', data[3][2], ',', data[7][2], '),\n');
			
			fprint(%c, %i, %c, %i, %c, 'avPink(', i, ',', data[4][2], ',', data[7][2], '),\n'); 
			fprint(%c, %i, %c, %i, %c, 'avPur(', i, ',', data[5][2], ',', data[7][2], '),\n'); 
			fprint(%c, %i, %c, %i, %c, 'avSky(', i, ',',  data[6][2], ',', data[7][2], ').\n'); 
		
			fprint(%c, %i, %c, %i, %c, 'hammingDistances(', i, ',',  data[0][2], ',', data[7][2], '),\n'); // TODO: update this

			fprint(%c, '\n');

		
		for( n = 0; n < this->c[i]->nodes[n]; n++)		
		{

			fprint(%c, %i, %c, '% node data from cluster: ', i, '\n');
	
			fprint(%c, %i, %c, %i, %c, %i, %c, 'node(', n, ',', i, ',' this->class,'):-\n');
		
			fprint(%c, %i, %c, %i, %c, %i, %c, %i, %c, 'Y(', n, ',', i, ',', this->c[i]->nodes[n]->Y, ',', this->class' ),\n'); 
			fprint(%c, %i, %c, %i, %c, %i, %c, %i, %c, 'Cr(', n, ',' , i,',', this->c[i]->nodes[n]->Cr, ',', this->class' ),\n'); 
			fprint(%c, %i, %c, %i, %c, %i, %c, %i, %c, 'Cb(', n, ',' , i, ',', this->c[i]->nodes[n]->Cb, ',', this->class' ),\n'); 
		
			fprint(%c, %i, %c, %i, %c, %i, %c, %i, %c, 'pink(', n,',' , i, ',', this->c[i]->nodes[n]->pink, ',', this->class' ),\n'); 
			fprint(%c, %i, %c, %i, %c, %i, %c, %i, %c, 'pur(', n, ',', i, ',', this->c[i]->nodes[n]->pur,',', this->class' ),\n'); 
			fprint(%c, %i, %c, %i, %c, %i, %c, %i, %c, 'sky(', n, ',', i, ',', this->c[i]->nodes[n]->sky,',', this->class' ),\n'); 
			
			fprint(%c, %i, %c, %i, %c, %i, %c, %i, %c, 'loc(', n, ',', i,',', this->c[i]->nodes[n]->loc,',', this->class' ).\n'); 
			
			fprint(%c, '\n');
			
		} // for all nodes
*/			
	} // for all clusters

	fclose(f);
		
} // pro


void cure::aggregate(quadrant * q, int * class)
{

	int * i = malloc(this->nq, sizeof(int));

	for(i = 0; i < this->nq; i++)
	{
		this->avY += q->Y;

			if(q->Y < this->minY) this->minY = q->minY;
			if(q->Y > this->maxY) this->maxY = q->Y;
	
		this->avCr += q->Cr; // weight more heavily as affected tissue appears to be pinkish red
	
			if(q->Cr < this->minCr) this->minCr = q->Cr;
			if(q->Cr > this->maxCr) this->maxCr = q->Cr;
	
		this->avCb += q->Cb;
	
			if(q->Cb < this->minCb) this->minCb = q->Cb;
			if(q->Cb > this->maxCb) this->maxCb = q->Cb;

		this->avPink += q->pink; // weight more heavily as affected tissue appears to be pinkish red
	
			if(q->pink < this->minPink) this->minPink = q->minPink;
			if(q->pink > this->maxPink) this->maxPink = q->maxPink;
	
		this->avPur += q->pur;
	
			if(q->pur < this->minPur) this->minPur = q->pur;
			if(q->pur > this->minPur) this->maxPur = q->pur;
	
		this->avSky += q->sky;
	
			if(q->sky < this->minSky) this->minSky = q->sky;
			if(q->sky > this->maxSky) this->maxSky = q->sky;

		this->tNodes += q->nNodes;

			if(q->nNodes < this->minNNodes) this->minNNodes = q->nNodes;
			if(q->nNodes > this->maxNNodes) this->maxNNodes = q->nNodes;

	} // for all quadrants

	this->avY /= this->nq;
	this->avCr /= this->nq;
	this->avCb /= this->nq;
	
	this->avPink /= this->nq;
	this->avPur /= this->nq;
	this->avSky /= this->nq;
	
	this->avNodes = this->tNodes/this->nq;

} // aggregate


/*
// largest data structure by far, must be hyperthreaded 
void cure::ham(quadrant * q){

		px * center = malloc(1, sizeof(px) );

		q->center = findMiddle( center, q );

		px * endPoints = malloc( (q->nNodes*2), sizeof(px) ); // approximates the line between two nodes as closely as possible			
				   // given the two chosen row and column variables 

		q->endPoints = getEnds( endPoints );

		edge * rnEdges = malloc( (q->nNodes-1), sizeof(edge) ); // linear mapping of the edges that 
						// theoretically exist between the geometrically central node and the outer nodes

		q->rnEdges = getRNEdges( q->center, q->endPoints, q->nodes );

		edge * neEdges = malloc( (q->nNodes * q->ne), sizeof(edge) ); // linear mapping of the edges that 
						// theoretically exist between the nodes and their respective neighbors

		q->neEdges = getNEEdges( q->nodes, q->nodes );

		edge * nnEdges = malloc( (q->nNodes*q->nNodes), sizeof(edge) ); // linear mapping of all edges between every node 
																		 // and every other node in the given quadrant

		q->nnEdges = getNNEdges( q->nodes );

		edge ** rnED = malloc( (q->nNodes-1), sizeof( malloc( 1, sizeof(edge) ) ) ); // radial node edge hamming distance structure

		edge ** neED = malloc( (q->nNodes * q->ne), sizeof( malloc( 1, sizeof(edge) ) ) ); // nearest neighbor hamming distance structure 

		edge ** nnED = malloc( (q->nNodes * q->nNodes), sizeof( malloc( 1, sizeof(edge) ) ) ); // intra nodal hamming distance structure


		q->rnEdgeHams = getHams( q->rnEdges, q->rnEdges );
		q->neEdgeHams = getHams( q->neHams, q->neHams );
		q->nnEdgeHams = getHams( q->nnHams, q->nnHams );



}; // hamming strings struct

*/

pair * cure::pair(int * x, int * a, int * b)
{
	pair * p = malloc(1, sizeof(pair));

	int * i = malloc(1, sizeof(int));
	i = x;
	int * f1 = malloc(1, sizeof(int) );
	int * f2 = malloc(1, sizeof(int) );

	if(a < b)
	{
		f1 = a;
		f2 = b;
	} // if
	else
	{
		f2 = a;
		f1 = b;
	} // else

} // pair


void cure::com(char * command){

	system(command);

} // com

}; // cure


// cure cancer, computer
int main(int argv, char ** argc)
{
	complete = false;

	else
	{ 
		cure * c = malloc(1, sizeof(cure) ); // allocating memory for the cure!
		c->begin.start(); // begin system clock for runtime comparison with effectiveness of classification

		char * i = argc[1]; // input file (presumed .jpg, as the test data is that way)	
		
		char * o = 'output.';
		char * p = '.csv';

		strcat(o, i);
		strcat(o, p);

		int * class = argc[2];
	
				strcpy( prolog, pro );
				strcat( prolog, '/' );
				strcat( prolog, m );
				strcat( prolog, sfx ); 
				
				c->i = input;
				c->o = output;
			
				// cure the shit out of it, bit!
			 
				c->bits();

				c->factors(c->n);

				int * x = malloc(1, sizeof(int));

				x = 0;

				c->nQ = rand() % (c->n/c->nFactors) ) + c->nFactors;
	
				quadrant * q = malloc( c->nQ, sizeof(quadrant) );

				c->q = q;
				
				while(x < c->nq){

						fact1 = this->factors[sz];
						fact2 = this->factors[sz];

						q[sz] = c->quadrant(fact1->f1, fact2->f1); // first factor from each pair of factors ensures low quadrant size, and there
																// fore a more meticulous coverage of the image's data
						c->iQ(q[x]);
						
					} // while

				c->lobe(c->q, c->nq, class);

				c->log();
	
/*
				c[x]->ham(c[x]->q, c[x]->nq); // init hamming distance data structures

				// TODO: parallelize 
				for(sz = 0; sz < c[x]->nq; sz++)
				{
					c[x]->intraHam(c[x]->q[sz]); // compute hamming distances for all quadrants in this cure
				}

					// TODO: hyperthread
				c[x]->interHam(c[x]->q, c[x]->nq); // compute hamming distances for every node in the lobe
				
				c[x]->log(c[x]->q); // write PROLOG data structures
				
				c[x]->com(); // call python sdk from command line

				c[x]->learn(c[x]->class); // evaluate results of python Deep Learning SDK, and modify PROLOG knowledge base if/when nec.


			*/	

			bool complete = c[x]->log(c[x]->q); // write PROLOG data structures

	} // if argv == argc[1]

	return complete;
	
} // main

	

//////// prolog knowledge base ////////////////////


/*

	C = class, I = index of cluster, D = data
	
	classify(C, cluster(I, D) ):-
	
	


*/