#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <math.h>

#define FIELDSIZE 10
#define CountMeasurePoints 6

typedef struct vectC{
	float x;
	float y;
	float c;
} vectC;

typedef struct Measure{
	unsigned count;
	vectC point[CountMeasurePoints];
} Measure;

void setRandPos(vectC * v){
	if (v){
		v->x = ((double)rand() / RAND_MAX) * FIELDSIZE;
		v->y = ((double)rand() / RAND_MAX) * FIELDSIZE;
	}
}

float calcDist(vectC * v1, vectC * v2){
	float dist = __FLT_MIN__;
	if (v1 && v2)
		dist = sqrt(pow((v1->x-v2->x),2)+pow((v1->y-v2->y),2));
	return dist;
}

float calcConcentration(vectC * src, vectC * MeasurePoint){
	float concenration = 0;
	if (src && MeasurePoint){
		concenration = src->c * exp(- calcDist(src, MeasurePoint));
	}
	return concenration;
}

void init(vectC * RealSource, Measure * msr){
	if (RealSource && msr){
		setRandPos(RealSource);
		RealSource->c = 1;
		for (unsigned i = 0; i < msr->count; i++){
			setRandPos(msr->point+i);
			msr->point[i].c = calcConcentration(RealSource, msr->point+i);
		}
	}
}

vectC * findVectWithMaxConcentration(Measure * msr){
	vectC * out = NULL;
	if (msr){
		if (msr->count > 0) out = msr->point;
		for (unsigned i = 1; i < msr->count; i++){
			if (msr->point[i].c > out->c) out = msr->point + i;
		}
	}
	return out;
}

void printVect(char* name, vectC * v){
	if (v) printf("%10s\t\t%.7f\t%.7f\t%.7f\n", name, v->x, v->y, v->c);
}

void print_Data(vectC * RealSource, Measure * msr, vectC * Vect_MaxC){
	printf("Name\t\t\tx\t\ty\t\tc\n");
	printVect("RealSource", RealSource);
	for (unsigned i = 0; msr && i < msr->count; i++)
		printVect("msr->point:", msr->point +i);
	printVect("Vect_MaxC", Vect_MaxC);
}

int main(){
	srand(time (NULL));
	vectC RealSource;
	Measure datameasure;
	datameasure.count = CountMeasurePoints;
	init(&RealSource, &datameasure);
	vectC * Vect_MaxC = findVectWithMaxConcentration(&datameasure);
	if (Vect_MaxC){
		print_Data(&RealSource, &datameasure, Vect_MaxC);
	}
	return 0;
}