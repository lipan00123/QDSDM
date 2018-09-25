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

// void projection_edge_weight(double * y, double * a, double * W, const double * para){
//     y[0] = (W[0]*a[0] - W[1]*a[1])/(W[0] + W[1]);
//     if(y[0] < 0) y[0] = fmax(-para[0], y[0]);
//     else y[0] = fmin(para[0], y[0]);
//     y[1] = - y[0];
// }
// 
// void projection_para_edge_weight(double * y, double * a, double * W, const double * para, int len){
//     int n = len/2, i, twoi;
//     for(i = 0; i < n;i++){
//         twoi = 2*i;
//         y[twoi] = (W[twoi]*a[twoi] - W[twoi+1]*a[twoi+1])/(W[twoi] + W[twoi+1]);
//         if(y[twoi] < 0) y[twoi] = fmax(-para[i], y[twoi]);
//         else y[twoi] = fmin(para[i], y[twoi]);
//         y[twoi+1] = - y[twoi];
//     }
// }


void proximal_undirected_hyperedge(double * x, double * a, double * W, const double * para, int len){
    double * sortW = new double [len]();
    double * weighta = new double [len]();
    double * index = new double [len]();
    double r, dr, wr, s, ds, ws, was, war, deltar, deltas;
    int i, ri, si, r_or_s;
    for(i = 0; i < len; i++){
        weighta[i] = a[i];
        index[i] = i;
    }
    MYsort(weighta, index, len, 'd');
    for(i = 0; i < len; i++){
        sortW[i] = W[(int)index[i]]/para[0]/para[0];
        //mexPrintf("%f,%f\n", weighta[i],sortW[i]);
    }
    //mexPrintf("%f\n", weighta[0]);
    //mexPrintf("%f\n", weighta[len-1]);
    r = weighta[0];
    dr = 0;
    ri = 0;
    wr = sortW[0];
    s = weighta[len-1];
    ds = 0;
    si = len-1;
    ws = sortW[len-1];
    //mexPrintf("r-s0: %f\n", r-s);
    while(r - s > accr){
        r_or_s = 1;
        deltar = r - weighta[ri+1];
        //mexPrintf("deltar:%f\n", deltar);
        deltas = deltar*wr/ws;
        if(s + deltas > weighta[si-1]){
            r_or_s = 2;
            deltas = weighta[si-1] - s;
            deltar = deltas*ws/wr;
        }
        r -= deltar;
       // mexPrintf("r:%f\n", r);
        dr += wr*deltar;
        s += deltas;
       // mexPrintf("s:%f\n", s);
        ds -= wr*deltar;
        if(dr > r - s) break;
        if(r_or_s == 1){
            ri ++;
            wr += sortW[ri];
        }else{
            si --;
            ws += sortW[si];
        }    
    }
    //mexPrintf("r-se:%f\n", r-s);
    if(fabs(r - s) < accr){
        for(i = 0; i < len; i++){
            x[i] = weighta[i];
        }
        
    }else{
        wr += 1;
        ws += 1;
        war = 0;
        was = 0;
        for(i=0; i <= ri; i++)
            war += sortW[i]*weighta[i];
        for(i=len-1; i >= si; i--)
            was += sortW[i]*weighta[i];
        r = (ws*war + was)/(ws*wr - 1);
        s = (wr*was + war)/(ws*wr - 1);
        for(i=0; i <= ri; i++)
            x[(int)index[i]] = r;
        for(i=len-1; i >= si; i--)
            x[(int)index[i]] = s;
        for(i = ri+1; i<si; i++)
            x[(int)index[i]] = weighta[i];
    }
    delete[] sortW;
    delete[] index;
    delete[] weighta;
}

double primal_val(int ** incidence_list, int * incidence_list_size,
        double ** parameter_list, int * parameter_list_size,
        char ** submodular_type, int * submodular_type_size,
        double * bias_vec, double * x,
        double *W, int N, int R){
    double gap=0, dv = 0, tempmax, tempmin;
    int i,j;
    for(i = 0; i < N; i++){
        gap +=  pow(x[i] - bias_vec[i],2)*W[i];
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
double gap_eff(int ** incidence_list, int * incidence_list_size,
        double ** parameter_list, int * parameter_list_size,
        char ** submodular_type, int * submodular_type_size,
        double * bias_vec, double * y_sum, double * phi, double * x,
        double *W, int N, int R){
    double gap=0, dv = 0, tempmax, tempmin;
    int i,j;
    for(i = 0; i < N; i++){
        gap +=  pow(x[i] - bias_vec[i], 2)*W[i] + y_sum[i]*y_sum[i]/W[i] - 2*y_sum[i]*bias_vec[i];
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
            gap += pow(phi[i],2);
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

void degree_stat(double * total_degree_vec, int ** incidence_list, int * incidence_list_size, int N, int R){
    int i, j;
    for(i = 0; i < N; i++){
        total_degree_vec[i] = 0;
    }
    for(i = 0; i < R; i++){
        for(j = 0; j < incidence_list_size[i]; j++){
            total_degree_vec[incidence_list[i][j]]++;
        }
    }
}


void main_func(double * x, double * record, double * final_gap,
        int ** incidence_list, int * incidence_list_size,
        double ** parameter_list, int * parameter_list_size,
        char ** submodular_type, int * submodular_type_size,
        double * bias_vec, double * W, int N, int R, unsigned int T,
        unsigned int record_dis){
    double ** y = new double * [R];
    double * a = new double [maxe];
    double * proximal_weight = new double [maxe];
    double * total_degree_vec = new double [N];
    double * y_sum = new double [N]();
    double * localx = new double [maxe];
    double * barx = new double [N]();
    int i, j, t;
    double gap;
    double * phi = new double [R]();
    double sigma, tau, theta = 1, tempysum, tempx, max_degree;
    
    degree_stat(total_degree_vec, incidence_list, incidence_list_size, N, R);
    max_degree = 0;
    for(i = 0; i < N; i++){
        max_degree = fmax(max_degree, total_degree_vec[i]);
    }
    sigma = 1/sqrt(max_degree + 1);
    tau = 1/sqrt(max_degree + 1);
    //mexPrintf("max_degree: %f\n", max_degree);
    
    for(i = 0; i < R; i++){
        y[i] = new double [incidence_list_size[i]]();
        phi[i] = 0;
    }
    for(i = 0; i < N; i++){
        x[i] = 0;
        barx[i] = 0;
    }
    
    for(i = 0; i < maxe; i++){
        proximal_weight[i] = sigma;
    }
    
    for(t = 0; t <= T; t++){
        if(t%record_dis == 0){
            gap = gap_eff(incidence_list, incidence_list_size,
                parameter_list, parameter_list_size,
                submodular_type, submodular_type_size,
                bias_vec, y_sum, phi, x, W, N, R);
//             gap = primal_val(incidence_list, incidence_list_size,
//                 parameter_list, parameter_list_size,
//                 submodular_type, submodular_type_size,
//                 bias_vec, x, W, N, R);  
            //mexPrintf("%d, gap: %f\n", t, gap);
            record[t/record_dis] = gap;
        }
        for(i = 0; i < N; i++){
            y_sum[i] = 0;
        }
        for(i = 0; i < R; i++){
            for(j = 0; j < incidence_list_size[i]; j++){
                a[j] = y[i][j] + sigma * barx[incidence_list[i][j]];
                //mexPrintf("a: %f\n", a[j]);
            }
            if(submodular_type[i][0] == 'h'){
                proximal_undirected_hyperedge(localx, a, proximal_weight, parameter_list[i], incidence_list_size[i]);
            }
            tempysum = 0;
            for(j = 0; j < incidence_list_size[i]; j++){
                //mexPrintf("localx: %f\n", localx[j]);
                y[i][j] = a[j] - localx[j];
                //mexPrintf("y: %f\n", y[i][j]);
                tempysum += fabs(y[i][j]);
                y_sum[incidence_list[i][j]] += y[i][j];
            }
            phi[i] = tempysum/2/parameter_list[i][0];
            //mexPrintf("phi: %f\n", phi[i]);
        }
        for(i = 0; i < N; i++){
            tempx = (tau*(W[i]*bias_vec[i] - y_sum[i]) + x[i])/(tau*W[i]+1);
            barx[i] = tempx + theta*(tempx - x[i]);
            x[i] = tempx;
            //mexPrintf("x: %f\n", x[i]);
            //mexPrintf("barx: %f\n", barx[i]);
        }
    }
    *final_gap = gap_eff(incidence_list, incidence_list_size,
            parameter_list, parameter_list_size,
            submodular_type, submodular_type_size,
            bias_vec, y_sum, phi, x, W, N, R);
    delete[] y;
    delete[] y_sum;
    delete[] localx;
    delete[] barx;
    delete[] proximal_weight;
    delete[] total_degree_vec;
    delete[] a;
    delete[] phi;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    // [y, record] = PDHG_QDSFM_cversion(incidence_list, parameter_list, submodular_type, bias_vec, W, N, R, T, record_dis);
    // Compared to RCDM: no parallel parameter K while have weight matrix W
	if (nlhs != 3 || nrhs != 9) {
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
    plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double * final_gap = mxGetPr(plhs[2]);
// no discrete gap for QDSFM    
//     plhs[2] = mxCreateDoubleMatrix(record_num, 1, mxREAL);
//     double * record_discrete = mxGetPr(plhs[2]);
    
    main_func(x, record, final_gap,
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







