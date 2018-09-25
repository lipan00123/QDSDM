#include <math.h>
#include <stdio.h>
#include <string>
#include <random>
#include <time.h>
#include <unordered_map>
#include "mex.h"
#include "matrix.h"

#define maxR 4096
#define maxN 65536
#define maxe 32768
// QRCDM does not support parallelization now. 
#define maxK 1
#define accr 1e-12

void swap_double(double * a, int p, int q){
	double temp = a[p];
	a[p] = a[q];
	a[q] = temp;
}
void swap_int(int * a, int p, int q){
	int temp = a[p];
	a[p] = a[q];
	a[q] = temp;
}

void MYsort(double * a, double * b, int len, char type){
	// sort a according to rule type, permutate b according to a 
	if(len <= 1) return;
 	std::default_random_engine generator;
	std::uniform_int_distribution<int> distribution(0, len-1);
	int pos = distribution(generator), p = 0;
	double pivot = a[pos];
    swap_double(a, pos, len-1);
    swap_double(b, pos, len-1);
	if(type == 'd'){
		for(int i = 0; i < len-1; i++){
			if(a[i] > pivot){
				swap_double(a, i, p);
				swap_double(b, i, p);
			    p++;
			}
		}
		swap_double(a, p, len-1);
		swap_double(b, p, len-1);
	}
	if(type == 'a'){
		for(int i = 0; i < len-1; i++){
			if(a[i] < pivot){
				swap_double(a, i, p);
				swap_double(b, i, p);
			    p++;
			}
		}
		swap_double(a, p, len-1);
		swap_double(b, p, len-1);
	}
	MYsort(a, b, p, type);
	MYsort(a+p+1, b+p+1, len-p-1, type);
}

void subgradient_homo(double * dx, double * x, int ** incidence_list,
        int * incidence_list_size, double ** parameter_list, int N, int R){
    int i, j, maxindex, minindex;
    double tempmax, tempmin, dis;
    for(i = 0; i < N; i++){
        dx[i] = 0;
    }
    for(i = 0; i < R; i++){
            tempmax = -INFINITY;
            maxindex = -1;
            tempmin = INFINITY;
            minindex = -1;
            for(j = 0; j < incidence_list_size[i]; j++){
                if(tempmax < x[incidence_list[i][j]]){
                    tempmax = x[incidence_list[i][j]];
                    maxindex = incidence_list[i][j];
                }
                if(tempmin > x[incidence_list[i][j]]){
                    tempmin = x[incidence_list[i][j]];
                    minindex = incidence_list[i][j];
                }
            }
            dis = 2*parameter_list[i][0]*(tempmax - tempmin);
            dx[maxindex] += dis;
            dx[minindex] -= dis;
    }
}

double primal_val(double * x, int ** incidence_list, int * incidence_list_size,
        double ** parameter_list, int * parameter_list_size,
        char ** submodular_type, int * submodular_type_size,
        double * bias_vec, double *W, int N, int R){
    double gap= 0, dv = 0, tempmax, tempmin;
    int i,j;
    for(i = 0; i < N; i++){
        gap += W[i]*pow(x[i] - bias_vec[i],2);
    }
    for(i = 0; i < R; i++){
        if(submodular_type[i][0] == 'h'){
            tempmax = -INFINITY;
            tempmin = INFINITY;
            for(j = 0; j < incidence_list_size[i]; j++){
                tempmax = fmax(tempmax, x[incidence_list[i][j]]);
                tempmin = fmin(tempmin, x[incidence_list[i][j]]);
            }
            gap += pow((tempmax - tempmin)*parameter_list[i][0],2);
        }
//         if(submodular_type[i][0] == 'e'){
//             gap += fabs(x[incidence_list[i][0]] - x[incidence_list[i][1]])*parameter_list[i][0];
//         }
//         if(submodular_type[i][0] == 'p'){
//             for(j = 0; j < incidence_list_size[i]/2; j++){
//                 gap += fabs(x[incidence_list[i][2*j]] - x[incidence_list[i][2*j+1]])*parameter_list[i][j];
//             }
//         }
//         if(submodular_type[i][0] == 'c'){
//             double * tempx = new double [incidence_list_size[i]];
//             double * index = new double [incidence_list_size[i]]();
//             for(j = 0; j < incidence_list_size[i]; j++){
//                 tempx[j] = x[incidence_list[i][j]];
//             }
//             MYsort(tempx, index, incidence_list_size[i], 'd');
//             for(j = 0; j < incidence_list_size[i]; j++){
//                 gap += tempx[j]*parameter_list[i][j];
//             }
//             delete[] tempx;
//             delete[] index;
//         }
    }
    return gap;
}

double gap_to_final_dual(double * x, int ** incidence_list, int * incidence_list_size,
        double ** parameter_list, int * parameter_list_size,
        char ** submodular_type, int * submodular_type_size,
        double * bias_vec, double *W, int N, int R, double dual_final){
    double gap= -  dual_final, dv = 0, tempmax, tempmin;
    int i,j;
    for(i = 0; i < N; i++){
        gap += W[i]*pow(x[i] - bias_vec[i],2);
    }
    for(i = 0; i < R; i++){
        if(submodular_type[i][0] == 'h'){
            tempmax = -INFINITY;
            tempmin = INFINITY;
            for(j = 0; j < incidence_list_size[i]; j++){
                tempmax = fmax(tempmax, x[incidence_list[i][j]]);
                tempmin = fmin(tempmin, x[incidence_list[i][j]]);
            }
            gap += pow((tempmax - tempmin)*parameter_list[i][0],2);
        }
//         if(submodular_type[i][0] == 'e'){
//             gap += fabs(x[incidence_list[i][0]] - x[incidence_list[i][1]])*parameter_list[i][0];
//         }
//         if(submodular_type[i][0] == 'p'){
//             for(j = 0; j < incidence_list_size[i]/2; j++){
//                 gap += fabs(x[incidence_list[i][2*j]] - x[incidence_list[i][2*j+1]])*parameter_list[i][j];
//             }
//         }
//         if(submodular_type[i][0] == 'c'){
//             double * tempx = new double [incidence_list_size[i]];
//             double * index = new double [incidence_list_size[i]]();
//             for(j = 0; j < incidence_list_size[i]; j++){
//                 tempx[j] = x[incidence_list[i][j]];
//             }
//             MYsort(tempx, index, incidence_list_size[i], 'd');
//             for(j = 0; j < incidence_list_size[i]; j++){
//                 gap += tempx[j]*parameter_list[i][j];
//             }
//             delete[] tempx;
//             delete[] index;
//         }
    }
    return gap;
}

void main_func(double * x, double * record,
        int ** incidence_list, int * incidence_list_size,
        double ** parameter_list, int * parameter_list_size,
        char ** submodular_type, int * submodular_type_size,
        double * bias_vec, double * W, int N, int R, unsigned int T,
        unsigned int record_dis){
    int i, t;
    double gap;
    double * dx = new double [N];
    double normdx;
    double maxW = 0;
    for(i = 0; i < N; i++){
        x[i] = 0;
        dx[i] = 0;
        maxW = fmax(maxW,W[i]);
    }
    for(t = 0; t <= T; t++){
        if(t%record_dis == 0){
//             gap = gap_to_final_dual(x, incidence_list, incidence_list_size,
//                     parameter_list, parameter_list_size, 
//                     submodular_type, submodular_type_size,
//                     bias_vec, W, N, R, dual_final);
            gap = primal_val(x, incidence_list, incidence_list_size,
                    parameter_list, parameter_list_size, 
                    submodular_type, submodular_type_size,
                    bias_vec, W, N, R);

            record[t/record_dis] = gap;
        }
        subgradient_homo(dx, x, incidence_list, incidence_list_size,
                    parameter_list, N, R);
        normdx = 0;
        for(i = 0; i < N; i++){
            dx[i] += 2*W[i]*(x[i] - bias_vec[i]);
            normdx += dx[i]*dx[i];
        }
        //mexPrintf("%f\n",normdx);
        for(i = 0; i < N; i++){
            x[i] -= 1.05*dx[i]/sqrt(normdx)/(t+1)/maxW;
        }
    }
    delete[] dx;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    // [y, record] = Subgradient_QDSFM_cversion(incidence_list, parameter_list, submodular_type, bias_vec, W, N, R, T, record_dis, dual_final);
    // Compared to RCDM: no parallel parameter K while have weight matrix W
	if (nlhs != 2 || nrhs != 9) {
    mexWarnMsgTxt("Check Parameters");
    return;
    }
    // read incidence_list
    int R = *(mxGetPr(prhs[6]));
    const mxArray * incidence_list_org = prhs[0];
    mxArray * incidence_Element;
    //const int * R_org = mxGetDimensions(incidence_list_org);
    int * * incidence_list = new int * [R];
    int * incidence_list_size = new int [R];
    double * templist;
    int j, k;
    for(j = 0; j < R;j++){
       incidence_Element = mxGetCell(incidence_list_org, j);
       incidence_list_size[j] = (int)mxGetN(incidence_Element);
       incidence_list[j] = new int [incidence_list_size[j]];
       templist = mxGetPr(incidence_Element);
       for(k = 0; k < incidence_list_size[j]; k++){
           incidence_list[j][k] = templist[k] - 1;
           //mexPrintf("%d, %d, %d\n", j, k, incidence_list[j][k]);
       }
    }
    // read parameter_list
    const mxArray * parameter_list_org = prhs[1];
    mxArray * parameter_Element;
    double * * parameter_list = new double * [R];
    int * parameter_list_size = new int [R];
    for(j = 0; j < R;j++){
       parameter_Element = mxGetCell(parameter_list_org, j);
       parameter_list_size[j] = (int)mxGetN(parameter_Element);
       parameter_list[j] = new double[parameter_list_size[j]];
       templist = mxGetPr(parameter_Element);
       for(k = 0; k < parameter_list_size[j]; k++){
           parameter_list[j][k] = templist[k];
           //mexPrintf("%d, %d, %f\n", j, k, parameter_list[j][k]);
       }
    }
    // read submodular_type
    const mxArray * submodular_type_org = prhs[2];
    mxArray * submodular_type_Element;
    char ** submodular_type = new char * [R];
    int * submodular_type_size = new int [R];
    char * temp_type;
    for(j = 0; j < R;j++){
       submodular_type_Element = mxGetCell(submodular_type_org, j);
       temp_type = (char *)mxGetPr(submodular_type_Element);
       submodular_type_size[j] = (int)mxGetN(submodular_type_Element);
       submodular_type[j] = new char[submodular_type_size[j]];
       for(k=0; k < submodular_type_size[j]; k++){
           submodular_type[j][k] = (char)temp_type[2*k];
           //mexPrintf("%d, %d, %c\n", j, k, temp_type[2*k]);
       }
    }
    
    double * bias_vec = mxGetPr(prhs[3]);
    double * W = mxGetPr(prhs[4]);
    int N = *(mxGetPr(prhs[5]));
    unsigned int T = *(mxGetPr(prhs[7]));
    unsigned int record_dis = *(mxGetPr(prhs[8]));
    int record_num = T/record_dis+1;
    
    plhs[0] = mxCreateDoubleMatrix(1, N, mxREAL);
    double * x = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(record_num, 1, mxREAL);
    double * record = mxGetPr(plhs[1]);
// no discrete gap for QDSFM    
//     plhs[2] = mxCreateDoubleMatrix(record_num, 1, mxREAL);
//     double * record_discrete = mxGetPr(plhs[2]);
    
    main_func(x, record,
            incidence_list, incidence_list_size,
            parameter_list, parameter_list_size,
            submodular_type, submodular_type_size,
            bias_vec, W, N, R, T, record_dis);
    delete[] incidence_list;
    delete[] incidence_list_size;
    delete[] parameter_list;
    delete[] parameter_list_size;
    delete[] submodular_type;
    delete[] submodular_type_size;
}







