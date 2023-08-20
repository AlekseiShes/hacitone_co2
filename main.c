#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <math.h>

#define FIELDSIZE 10
#define CountMeasurePoints 6
#define CountOrient 10
#define CountCross 45

typedef struct vectC{
	double x;
	double y;
	double c;
} vectC;

typedef struct orient{
	vectC pos;
	vectC dir;
} orient;

typedef struct Measure{
	unsigned count;
	vectC point[CountMeasurePoints];
} Measure;

short unsigned IsNotZero(double a){ return fabs(a)>1E-16; }

void setRandPos(vectC * v){
	if (v){
		v->x = ((double)rand() / RAND_MAX) * FIELDSIZE;
		v->y = ((double)rand() / RAND_MAX) * FIELDSIZE;
	}
}

double calcDist(vectC * v1, vectC * v2){
	double dist = __FLT_MAX__;
	if (v1 && v2)
		dist = sqrt(pow((v1->x-v2->x),2)+pow((v1->y-v2->y),2));
	return dist;
}

double calcConcentration(vectC * src, vectC * MeasurePoint){
	double concenration = 0;
	if (src && MeasurePoint){
		concenration = src->c * exp(- calcDist(src, MeasurePoint)/10);
	}
	return concenration;
}

void init(vectC * RealSource, Measure * msr){
	if (RealSource && msr){
		// setRandPos(RealSource);
		// RealSource->x += 10;
		// RealSource->c = 1;
		// for (unsigned i = 0; i < msr->count; i++){
		// 	setRandPos(msr->point+i);
		// 	msr->point[i].c = calcConcentration(RealSource, msr->point+i);
		// }
		
		RealSource->x = 20;
		RealSource->y = 20;
		RealSource->c = 1;
		int i = 0;
		(msr->point+i)->x = 4;	(msr->point+i)->y = 8;	(msr->point+i)->c = 0;	i++;
		(msr->point+i)->x = 3;	(msr->point+i)->y = 6;	(msr->point+i)->c = 0;	i++;
		(msr->point+i)->x = 8;	(msr->point+i)->y = 5;	(msr->point+i)->c = 0;	i++;
		(msr->point+i)->x = 4;	(msr->point+i)->y = 4;	(msr->point+i)->c = 0;	i++;
		(msr->point+i)->x = 6;	(msr->point+i)->y = 3;	(msr->point+i)->c = 0;	i++;
		(msr->point+i)->x = 1;	(msr->point+i)->y = 2;	(msr->point+i)->c = 0;	i++;
		for (unsigned i = 0; i < msr->count; i++){
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
	if (v) printf("%10s\t%10.7f   %10.7f   %10.7f\n", name, v->x, v->y, v->c);
}

void print_Data(vectC * RealSource, Measure * msr, vectC * Vect_MaxC, orient *orients, vectC * Cross, vectC* Arr){
	printf("Name\t\tx\t\ty\t\tc\n");

	for (unsigned i = 0; msr && i < msr->count; i++)
		printVect("msr->point:", msr->point +i);
	for (unsigned i = 0; orients && i < CountOrient; i++){
		printVect("____pos", &(orients[i].pos));
		printVect("____dir", &(orients[i].dir));
		printf("\n");
	}
	for (unsigned i = 0; i < CountCross; i++){
		printVect("Arr    ", &(Arr[i]));
	}
	printVect("Vect_MaxC", Vect_MaxC);
	printVect("Cross     ", Cross);
	printVect("RealSource", RealSource);
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

void multvectnum(vectC * res, vectC *a, double num){
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

double scalarmult(vectC *a, vectC *b){
	double res=__FLT_MAX__;
	if (a && b)
		res = a->x *b->x + a->y*b->y + a->c*b->c;
	return res;
}

void normalize(vectC *a){
	if (a){
		double length = scalarmult(a, a);
		if (IsNotZero(length)){
			length = 1/sqrt(length);			
			multvectnum(a, a, length);
		}
	}
}

void projectVect2Surf(vectC *res, vectC *v, vectC *Ns, vectC *Ps){
    multvectnum(res, res, 0);
    summvect(res,res,Ns);
    double num = - scalarmult(Ns,Ns);
    if (IsNotZero(num)) num=1/num;
    num*=scalarmult(Ns,v);
    multvectnum(res, res, num);
    summvect(res,res,Ps);
    summvect(res,res,v);
}

void calcdir(vectC ** vecset, orient *orients){
	if (vecset && *vecset && orients){
		orients->dir.x = orients->dir.y = orients->dir.c = 0;
		vectC v1, v2, N;
		sub_vect(&v1, vecset[1], vecset[0]);
		sub_vect(&v2, vecset[2], vecset[0]);
		vectmult(&N, &v1, &v2);
		if (N.c > 0) multvectnum(&N, &N, -1);
		N.c = 0;
		//normalize(&N);
	    //vectC Ns={0,0,1}, res={0,0,0};
		//projectVect2Surf(&res, &N, &Ns,  &(orients->pos));
		//summvect(&(orients->dir), &(orients->dir), &res);


		// vectC real = {12, 12, 0};
		// vectC temppos ={0,0,0};
		// summvect(&temppos, &(orients->pos), NULL);
		// temppos.c = 0;
		// sub_vect(&v1, &temppos, &real);
		// summvect(&v2, &N, NULL);
		// normalize(&v1);
		// normalize(&v2);
		// vectC temp;
		// vectmult(&temp, &v1, &v2);
		// temp.c = 0;
		// if (! IsNotZero( scalarmult(&temp, &temp)))
		// 	N.x = N.y =  N.c= __FLT_MAX__;
		
		summvect(&(orients->dir), &(orients->dir), &N);
		
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
					if (io < CountOrient) calcOrient(vecset, orients+io);
					io++;
				}
			}
		}
	}
}

void OrientInit_2(orient*orients, Measure * msr, vectC * Vect_MaxC){
	if (orients && msr && Vect_MaxC){
		vectC * vecset[3] = {Vect_MaxC, NULL, NULL};
		int io = 0;
		for (unsigned i = 0; i < msr->count; i++){
			int needcalc = 1;
			if (msr->point+i != Vect_MaxC){
				vecset[1] = msr->point+i;
				for (unsigned j = i+1; j < msr->count; j++){
					if (msr->point+j != Vect_MaxC){
						vecset[2] = msr->point+j;
						if(needcalc){
							if (io < CountOrient)
								calcOrient(vecset, orients+io);
							io++;
						}
					}
				}
			}
		}
	}
}

void SetCAsZero(orient*orients){
	for (unsigned i = 0; orients && i < CountOrient; i++)
		orients->dir.c = orients->pos.c = 0;
}

void calcIcross(vectC * cross, orient * a, orient * b){
	if (cross && a && b){
		cross->x = cross->y = __FLT_MAX__; cross->c = 0;
		if (IsNotZero(a->dir.x) && IsNotZero(b->dir.x)){
			double ka = a->dir.y/a->dir.x,
						kb = b->dir.y/b->dir.x;
			if (IsNotZero(ka - kb)){
				cross->x = (b->pos.y - b->pos.x * kb - (a->pos.y - a->pos.x * ka )) / (ka - kb);
				cross->y = ka * (cross->x - a->pos.x) + a->pos.y;
			}
		}
		vectC v1, v2, N;
		sub_vect(&v1, &(a->pos), cross);
		sub_vect(&v2, &(a->dir), NULL);
		double num = scalarmult(&v1, &v2);
		if (num < 0)
			cross->x = cross->y = cross->c = __FLT_MAX__;
	}
}

void calcCross(vectC *res, orient*orients,vectC *Arr ){
	if (res && orients && Arr){
		unsigned k = 0;
		for (unsigned i = 0; i < CountOrient; i++){
			for (unsigned j = i+1; j < CountOrient; j++){
				//if (orients->dir.x )
				//без упора в точку
				//без  расходящ векторов
				//больше треугольников
				    calcIcross(Arr+k, orients+i, orients+j);
				k++;
			}
		}
		multvectnum(res, res, 0);
		int cnt = 0;
		for (unsigned i = 0; i < CountCross; i++){
			if ((Arr+i)->c < __FLT_MAX__-1){
				summvect(res, res, Arr+i);
				cnt++;
			}
		}
		if (cnt>0)
			multvectnum(res, res, 1.0/(cnt+0.0));
	}
}

void debug(vectC *res){
	vectC Arr[1];
	orient orient[2];
	//double a = sqrt(2)/2;
	orient->dir.x = 3;
	orient->dir.y = 5;
	orient->dir.c = 0;
	orient->pos.x = -1;
	orient->pos.y = 0;
	orient->pos.c = 0;
	(orient+1)->dir.x = -3;
	(orient+1)->dir.y = 5;
	(orient+1)->dir.c = 0;
	(orient+1)->pos.x = 1;
	(orient+1)->pos.y = 0;
	(orient+1)->pos.c = 0;
	calcIcross(Arr, orient, orient+1);
}

void checkOrients(orient*orients, vectC * RealSource){
	for (unsigned i = 0; i < CountOrient; i++) {
		vectC v1, v2, N;
		sub_vect(&v1, &(orients[1].pos), RealSource);
		sub_vect(&v2, &(orients[1].dir), NULL);
		vectmult(&N, &v1, &v2);
		double num = scalarmult(&N, &N);
		if (fabs(num)<1E-7){
			printf("Is___Zero( scalarmult(&N, &N)) for i = %u, num =%.15lf\n", i, num);
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
		orient orients[CountOrient] = {0};
		//OrientInit(orients, &datameasure, Vect_MaxC);
		OrientInit_2(orients, &datameasure, Vect_MaxC);
		SetCAsZero(orients);
		//checkOrients(orients, &RealSource);
		vectC Cross= {0,0,0};
		vectC Arr[CountCross];
		calcCross(&Cross, orients, Arr);
		print_Data(&RealSource, &datameasure, Vect_MaxC, orients, &Cross, Arr);
		//debug(&Cross);
	}
	return 0;
}