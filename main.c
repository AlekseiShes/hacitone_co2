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

typedef struct orient{
	vectC pos;
	vectC dir;
} orient;

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
		RealSource->x += 10;
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

void print_Data(vectC * RealSource, Measure * msr, vectC * Vect_MaxC, orient *orients){
	printf("Name\t\t\tx\t\ty\t\tc\n");
	printVect("RealSource", RealSource);
	for (unsigned i = 0; msr && i < msr->count; i++)
		printVect("msr->point:", msr->point +i);
	for (unsigned i = 0; msr && i < CountMeasurePoints-1; i++){
		printVect("____dir", &(orients[i].dir));
		printVect("____pos", &(orients[i].pos));
	}
	printVect("Vect_MaxC", Vect_MaxC);
}

void summvect(vectC * res, vectC *a, vectC *b){
	if (res && (a || b)){
		if (a){
			res->x = a->x + (b ? b->x: 0);
			res->y = a->y + (b ? b->y: 0);
			res->c = a->c + (b ? b->c: 0);
		} else {
			summvect(res, b, a);
		}
	}
}

void multvectnum(vectC * res, vectC *a, float num){
	if (res){
		res->x = (a? a->x: res->x) * num;
		res->y = (a? a->y: res->y) * num;
		res->c = (a? a->c: res->c) * num;
	}
}

void sub_vect(vectC * res, vectC *a, vectC *b){
	if (res && a && b){
		summvect(res, b, NULL);
		multvectnum(res, res, -1);
		summvect(res, res, a);
	}
}

void calcpos(vectC ** vecset, vectC * pos){
	if (vecset && *vecset && pos){
		pos->x = pos->y = pos->c = 0;
		summvect(pos, pos, vecset[0]);
		summvect(pos, pos, vecset[1]);
		summvect(pos, pos, vecset[2]);
	}
	multvectnum(pos, pos, 1.0/3.0);
}

void vectmult(vectC * res, vectC *a, vectC *b){
	if (res && a && b){
		res->x = a->y * b->c - a->c * b->y;
		res->y = a->c * b->x - a->x * b->c;
		res->c = a->x * b->y - a->y * b->x;
	}
}

float scalarmult(vectC *a, vectC *b){
	float res=FLT_MIN;
	if (a && b)
		res = a->x*b->y + a->y*b->y + a->c*b->c;
	return res;
}

void projectVect2Surf(vectC *res, vectC *v, vectC *Ns, vectC *Ps){
    multvectnum(res, res, 0);
    summvect(res,res,Ns);
    float num = - scalarmult(Ns,Ns);
    if (fabs(num)<1E-7) num=1/num;
    num*=scalarmult(Ns,v);
    multvectnum(res, res, num);
    summvect(res,res,Ps);
    summvect(res,res,v);
}


void calcdir(vectC ** vecset, orient *orients){
	if (vecset && *vecset && orients){
		orients->dir.x = orients->dir.y = orients->dir.c = 0;
		vectC v1, v2, N;
		sub_vect(&v1, vecset[0], vecset[1]);
		sub_vect(&v2, vecset[0], vecset[2]);
		vectmult(&N, &v1, &v2);
		if (N.c > 0) multvectnum(&N, &N, -1);
	    vectC Ns={0,0,1};
	    vectC res={0,0,0};
		projectVect2Surf(&res,&N, &Ns,  &(orients->pos));
		summvect(&(orients->dir), &(orients->dir), &res);
	}

}

void calcOrient(vectC ** vecset, orient*orients){
	if (vecset && *vecset && orients){
		calcpos(vecset, &(orients->pos));
		calcdir(vecset, orients);
	}
}

void OrientInit(orient*orients, Measure * msr, vectC * Vect_MaxC){
	if (orients && msr && Vect_MaxC){
		vectC * vecset[3] = {Vect_MaxC, NULL, NULL};
		int io = 0;
		for (unsigned i = 0; i < msr->count; i++){
			int needcalc = 1;
			if (msr->point+i != Vect_MaxC){
				vecset[1] = msr->point+i;
				if (i < msr->count-1 && msr->point+i+1 != Vect_MaxC)
					vecset[2] = msr->point+i+1;
				else if (i < msr->count-2)
					vecset[2] = msr->point+i+2;
				else if (msr->point != Vect_MaxC)
					vecset[2] = msr->point;
				else if(msr->count>1) vecset[2] = msr->point+1;
				else needcalc = 0;
				if(needcalc){
					if (io < CountMeasurePoints-1) calcOrient(vecset, orients+io);
					io++;
				}
			}
		}
	}
}

int main(){
	srand(time (NULL));
	vectC RealSource;
	Measure datameasure;
	datameasure.count = CountMeasurePoints;
	init(&RealSource, &datameasure);
	vectC * Vect_MaxC = findVectWithMaxConcentration(&datameasure);
	if (Vect_MaxC){
		orient orients[CountMeasurePoints-1] = {0};
		OrientInit(orients, &datameasure, Vect_MaxC);
		//
		print_Data(&RealSource, &datameasure, Vect_MaxC, orients);
	}
	return 0;
}