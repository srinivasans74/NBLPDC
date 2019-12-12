
#include "stdio.h"
#include "stdlib.h"
#include <math.h>
#include <time.h>
#include <string.h>
//#include <ap_fixed.h>
#include "nbldpcq16.h"
#ifdef VS
#define _CRT_SECURE_NO_DEPRECATE
#endif

//===================================
// Allocate 2D array of double
//===================================
double **malloc2Ddouble(int arraySizeX, int arraySizeY) // allocates array[a][b]
{
	double** theArray;  
	theArray = (double**) malloc(arraySizeX*sizeof(double*));  
	for (int i = 0; i < arraySizeX; i++)  
		theArray[i] = (double*) malloc(arraySizeY*sizeof(double));  
	return theArray; 
}

//===================================
// Galios Field Multiplication
//===================================
unsigned char mul_gf(unsigned char a, unsigned char b)
{
#pragma HLS RESOURCE variable=logq core=ROM_1P_BRAM
#pragma HLS RESOURCE variable=expq core=ROM_1P_BRAM
	if (a == 0 || b == 0) return 0;
	if (a == 1) return b;
	if (b == 1) return a;
	if ((logq[a] + logq[b])>(Q-1))
			{
			return	expq[(logq[a] + logq[b])- (Q-1)];
			}
			else
			return expq[(logq[a] + logq[b])];	//table look-up multiplication
}

//===================================
// Galios Field Inverse
//===================================
unsigned char inv_gf(unsigned char a)
{
#pragma HLS RESOURCE variable=logq core=ROM_1P_BRAM
#pragma HLS RESOURCE variable=expq core=ROM_1P_BRAM
	//if(a == 0) exit (EXIT_FAILURE);
	if(a == 0) return 0;
	if(a == 1) return 1;
	return expq[(Q - 1 - logq[a])];
	//return expq[(Q - 1 - logq[a]) % (Q-1)];

}

//===================================
// Galios Field Multiplication
//===================================
int div_gf(int a, int b)
{

	if (a == 0) return 0;
	if (b == 1) return a;
	return expq[(logq[a] - logq[b] + (Q-1))% (Q-1)];	//table look-up multiplication
}

//===================================
// Forward function
//===================================
#ifdef TREE
double forward(int symbol, int row, int col, int index)
{

	int a,b;
	double max_value, min_value;
	int index_;

	index_ = index - 1;
	col = col_row[row][index_];
	max_value = 0.0;
	min_value = 100.0;

	if(index_ == 0){	// reach the first non-zero element in current row
		return ALPHAmn[row][col][add_gf((H_nb[row][col]),symbol)];
	}
	else{
		for(b = 0; b < Q; b++){
			a = sub_gf(symbol,add_gf(H_nb[row][col],b));
			max_value = max(forward(a, row, col, index_),ALPHAmn[row][col][b]);
			min_value = min(min_value, max_value);
		}
		return min_value;
	}
}
#else
void forward(int row, int col, int index, double* F)
{
	int index_, col_;
	int a,b,c,i;
	double min,max;
	double* F_;

	F_ = (double*)malloc(sizeof(double) * Q);
	max = min = F[0];

	for(i = 0; i < Q; i++){
		F_[i] = F[i];
	}

	for(index_ = 1, col_ = col_row[row][0]; index_ <= index ; index_++, col_ = col_row[row][index_] ){
		for(c = 0; c < Q; c++){
			for(b = 0; b < Q; b++){
				a = sub_gf(c,add_gf(H_nb[row][col],b));
				max = max(F[a],ALPHAmn[row][col_][b]);
				min = min(min,max);
			}
			F[c] = min;
		}
	}
}
#endif

//===================================
// Backward function
//===================================
#ifdef TREE
double backward(int symbol, int row, int col, int index)
{
	int a,b;
	double max_value, min_value;
	int index_;

	index_ = index + 1;
	col = col_row[row][index_];
	max_value = 0.0;
	min_value = 100.0;

	if(index_ == row_weight[row]-1){	// reach the last non-zero element in current row
		return ALPHAmn[row][col][add_gf((H_nb[row][col]),symbol)];
	}
	else{
		for(b = 0; b < Q; b++){
			a = sub_gf(symbol,add_gf(H_nb[row][col],b));
			max_value = max(backward(a, row, col, index_),ALPHAmn[row][col][b]);
			min_value = min(min_value, max_value);
		}
		return min_value;
	}
}
#else
void backward(int row, int col, int index, double* B)
{
	int index_, col_;
	int a,b,c,i;
	double min,max;
	double* B_;

	B_ = (double*)malloc(sizeof(double) * Q);
	max = min = B[0];

	for(i = 0; i < Q; i++){
		B_[i] = B[i];
	}

	for(index_ = row_weight[row]-1, col_ = col_row[row][index_]; index_ >= index ; index_--, col_ = col_row[row][index_] ){
		for(c = 0; c < Q; c++){
			for(b = 0; b < Q; b++){
				a = sub_gf(c,add_gf(H_nb[row][col],b));
				max = max(B[a],ALPHAmn[row][col_][b]);
				min = min(min,max);
			}
			B[c] = min;
		}
	}
}
#endif



//===================================
// Root function
//===================================
#ifdef TREE
double root(int symbol, int row, int col, int index)
{
	int a,b;
	double max_value,min_value;

	max_value = 0.0;
	min_value = 100.0;

	for(a = 0; a < Q; a++)
	{
		b = add_gf(add_gf(H_nb[row][col],symbol),a);
		max_value = max(forward(a, row, col, index), backward(b, row, col, index));
		min_value = min(min_value, max_value);
	}
	return min_value;
}
#else
void root(int row, int col, int index, double* temp)
{
	int a,b,c;
	double min,max;
	int index_f,index_b;
	int col_f,col_b;
	double* F;
	double* B;

	index_f = index - 1;
	index_b = index + 1;

	col_f = col_row[row][index_f];
	col_b = col_row[row][index_b];

	F = (double *)malloc(sizeof(double) * Q);
	B = (double *)malloc(sizeof(double) * Q);

	forward(row, col_f, index, F);
	backward(row, col_b, index, B);

	max = min = F[0];

	for(c = 0; c < Q; c++){
		for(b = 0; b < Q; b++){
			a = sub_gf(c,add_gf(H_nb[row][col],b));
			max = max(F[a],B[b]);
			min = min(min,max);
		}
		temp[c] = min; 
	}
}
#endif

//===================================
// Minmax nonbinary decode function
//===================================
// Refer to "Min-max decoding for nonbinary LDPC codes"
int minmax(unsigned char BETAmn_[U][V][Q], unsigned char ALPHAmn_[U][V][Q], unsigned char GAMMAn_[Q][V], int* iteration, int* decoded_bit)
{


	int row,col,index,index_c,index_v,col_v,row_v,error;
	int index_, row_;
	int last_col;
	int index_B, col_B;
	int index_p, row_p, col_p;
	int row_weight_;
	int i,j;
	int a,b,c;
	int iter = 0;
	int min_index;
	//int syndrome;
	double min_F,min_B;
	double max_F,max_B;
	double temp;
	double min_value,max_value;
	double ALPHA_t[Q];
	int codeword_sym[V];
	double B[Q][Y];
	double F[Q][Y];
#ifdef DEBUG
	FILE *fp = fopen("C:/Dropbox/codes/non_binary_LDPC/vs/minmax/codeword.txt","w");
	FILE *fp1 = fopen("C:/Dropbox/codes/non_binary_LDPC/vs/minmax/alpha.txt","w");
	FILE *fp2 = fopen("C:/Dropbox/codes/non_binary_LDPC/vs/minmax/gamma_post.txt","w");
	FILE *fp3 = fopen("C:/Dropbox/codes/non_binary_LDPC/vs/minmax/gamma.txt","w");
	FILE *fp4 = fopen("C:/Dropbox/codes/non_binary_LDPC/vs/minmax/beta.txt","w");
	FILE *fp5 = fopen("C:/Dropbox/codes/non_binary_LDPC/vs/minmax/F.txt","w");
	FILE *fp6 = fopen("C:/Dropbox/codes/non_binary_LDPC/vs/minmax/B.txt","w");
	FILE *fp7 = fopen("C:/Dropbox/codes/non_binary_LDPC/vs/minmax/alpha_t.txt","w");
	FILE *fp8 = fopen("C:/Dropbox/codes/non_binary_LDPC/vs/minmax/beta_t.txt","w");
#endif





#ifdef DEBUG
	for(col = 0; col < N; col++){
		for(a = 0; a < Q; a++){
			fprintf(fp3, "%3f ", GAMMAn_[a][col]);
		}
		fprintf(fp3, "\n");
	}
#endif

	for(i = 0; i < Q; i++){
		for(j = 0; j < Y; j++){
			F[i][j] = B[i][j] = 0.0;
		}
	}

	for(i = 0; i < Q; i++){
		for(j = 0; j < V; j++){

			GAMMAn_post[i][j] = 0.0;
		}
	}





	minmax_label1:for(row = 0; row < U; row++)// initialize variable node message ALPHA with channel info GAMMA
	{
		col = col_row[row][0];
		minmax_label26:for(index = 0; index < row_weight[row] ; index++)
		{
#pragma HLS LOOP_TRIPCOUNT min=310 max=620

			col = col_row[row][index];
			for(a = 0; a < Q; a++){
				ALPHAmn_[row][col][a] = GAMMAn_[a][col];
			}
		}
	}

#ifdef PRINT_A_B
	for(row = 0; row < M; row++){
		for(col = 0; col < N; col++){
			for(a = 0; a < Q; a++){
				fprintf(fp1, "%3f ", ALPHAmn_[row][col][a]);
			}
			fprintf(fp1, "\n");
		}
	}
#endif

#ifdef DEBUG
	for(col = 0; col < N; col ++){

		fprintf(fp,"%3d ",encoded_sym[col]);
	}
	fprintf(fp, "\n");
#endif

	criticalpath1:
	for(col = 0; col < V; col ++){ // tentative decoding
		a = 0;
		for(i = 0; i < Q; i++){


#ifdef DEBUG
			fprintf(fp3,"%3f ",GAMMAn_[i][col] );
#endif
			if(GAMMAn_[i][col] < GAMMAn_[a][col])
			{

				a = i;
			}

		}

		codeword_sym[col] = a;
#ifdef DEBUG
		fprintf(fp,"%3d ",codeword_sym[col]);
		fprintf(fp3, "\n");
#endif
	}

	error = 0;

	for(i = 0;i < (U);i ++){
		if(encoded_sym[i + U] != codeword_sym[i + U])
			error++;
	}

#ifdef DEBUG
#else
	if(error == 0) {
#ifdef _NOSYNTH_
		printf("This received codeword has no error.\n");
		for(i = 0; i < Q; i++){
			free(B[i]);
		}
		free(B);

		for(i = 0; i < Q; i++){
			free(F[i]);
		}
		free(F);

		free(ALPHA_t);
		free(codeword_sym);

		for(i = 0; i < Q; i++){
			free(GAMMAn_post[i]);
		}
		free(GAMMAn_post);
#endif
		return 0;
	}
#endif

#ifdef DEBUG
	fprintf(fp,"error iss %d\n", error);

	for(row = 0; row < M; row ++){
		syndrome = 0;
		for(col = 0; col < N; col ++){
			syndrome = add_gf(syndrome,mul_gf(H_nb[row][col],codeword_sym[col]));
		}
		printf("row = %d syndrome = %d\n",row, syndrome);
	}
#endif	

	minmax_label2:while(iter < MAX_ITERATION){
	// iteration starts


		minmax_label3:for(row = 0; row <U; row++){
	// check node processing: update check node message BETA


#ifdef PRINTMSG			
			printf("row = %d\n",row);
#endif
			row_weight_ = row_weight[row];
			minmax_label4:for(index_c = 0, col = col_row[row][0]; index_c < row_weight_; index_c++, col = col_row[row][index_c]){
#pragma HLS LOOP_TRIPCOUNT min=310  max=620
				// 1. Initialize F and B value for each row
				if(index_c == 0){
					minmax_label5:for(a = 0; a < Q; a++){

						//F[a][0] = ALPHAmn_[row][col][mul_gf(div_gf(1,H_nb[row][col]),a)];
						F[a][0] = ALPHAmn_[row][col][mul_gf(inv_gf(H_nb[row][col]),a)];
						//F[a][0] = ALPHAmn_[row][col][a];
#ifdef DEBUG
						if(row == 58){
							fprintf(fp5,"F[%d][0] = %3.3f b = %d\n",a,F[a][0],mul_gf(inv_gf(H_nb[row][col]),a));
						}
#ifdef PRINTMSG
						printf("F[%d][0] = %3.3f b = %d\n",a,F[a][0],mul_gf(inv_gf(H_nb[row][col]),a));
#endif
#endif
						last_col = col_row[row][row_weight_ - 1];
						//B[a][row_weight_ - 1] = ALPHAmn_[row][last_col][mul_gf(div_gf(1,H_nb[row][last_col]),a)];
						B[a][row_weight_ - 1] = ALPHAmn_[row][last_col][mul_gf(inv_gf(H_nb[row][last_col]),a)];
						//B[a][row_weight_ - 1] = ALPHAmn_[row][last_col][a];
#ifdef DEBUGl
						if(row == 58){
							fprintf(fp6,"B[%d][%d] = %3.3f\n",a,(row_weight_ - 1),B[a][row_weight_ - 1],mul_gf(inv_gf(H_nb[row][last_col]),a));
						}
#ifdef PRINTMSG
						printf("B[%d][%d] = %3.3f H[%d][%d] = %d , inv(H[%d][%d]) = %d inv(H)* %d = %d \n",a,(row_weight_ - 1),B[a][row_weight_ - 1],row, last_col, H_nb[row][last_col],row, last_col, inv_gf(H_nb[row][last_col]), a, mul_gf(inv_gf(H_nb[row][last_col]),a));
#endif
#endif
					}
				}
				else{
#ifdef PRINTMSG
					printf("row = %d col = %d index_c = %d\n",row,col,index_c);
#endif
					index_p = index_c - 1;
					col_p = col_row[row][index_p];
					minmax_label6:for(c = 0; c < Q; c++){

#ifdef PRINTMSG
						printf("c = %d\n",c);
#endif
						min_F = max(F[c][index_c - 1],ALPHAmn_[row][col][0]);
#ifdef PRINTMSG
						printf("min_F = %3.3f\n",min_F);
#endif
						minmax_label7:for(b = 0; b < Q; b++){



							a = sub_gf(c,mul_gf(H_nb[row][col],b));
							//a = mul_gf(sub_gf(c,mul_gf(H_nb[row][col],b)),inv_gf(H_nb[row][col_p]));
#ifdef PRINTMSG
							printf("a = %d b = %d\n",a,b);
							printf("F[%d][%d] = %3.3f ALPHAmn_[%d][%d][%d] = %3.3f\n",a,index_c - 1, F[a][index_c - 1],row, col, b, ALPHAmn_[row][col][b]);
#endif
							//max_F = max(F[mul_gf(H_nb[row][col],b)][index_c - 1], ALPHAmn_[row][col][a]);
							max_F = max(F[a][index_c - 1], ALPHAmn_[row][col][b]);
							min_F = min(min_F, max_F);
#ifdef PRINTMSG
							printf("max_F = %3.3f min_F = %3.3f\n",max_F,min_F);
#endif
						}
						F[c][index_c] = min_F;
#ifdef DEBUG
						if(row == 58){
							fprintf(fp5,"F[%d][%d] = %3.3f\n",c,index_c,F[c][index_c]);
						}
#endif

#ifdef PRINTMSG
						printf("F[%d][%d] = %3.3f\n",c,index_c,F[c][index_c]);
#endif
					}

					index_B = row_weight_ - index_c - 1;
					col_B = col_row[row][index_B];
#ifdef PRINTMSG
					printf("index_B = %d col_B = %d\n",index_B,col_B);
#endif
					index_p = index_B + 1;
					col_p = col_row[row][index_p];
					minmax_label8:for(c = 0; c < Q; c++){


						min_B = max(ALPHAmn_[row][col_B][0],B[c][index_B + 1]);
#ifdef PRINTMSG
						printf("c = %d\n",c);
						printf("min_B = %3.3f\n",min_B);
#endif
						minmax_label9:for(b = 0; b < Q; b++){
							a = sub_gf(c,mul_gf(H_nb[row][col_B],b));
							//a = mul_gf(sub_gf(c,mul_gf(H_nb[row][col],b)),inv_gf(H_nb[row][col_p]));
#ifdef PRINTMSG
							printf("a = %d b = %d\n",a,b);
							printf("B[%d][%d] = %3.3f ALPHAmn_[%d][%d][%d] = %3.3f\n",a,index_B + 1,B[a][index_B + 1],row,col_B,b,ALPHAmn_[row][col_B][b]);
#endif
							max_B = max(B[a][index_B + 1], ALPHAmn_[row][col_B][b]);
							//max_B = max(B[mul_gf(H_nb[row][col_B],b)][index_B + 1], ALPHAmn_[row][col_B][a]);

							min_B = min(min_B, max_B);
#ifdef PRINTMSG
							printf("max_B = %3.3f min_B = %3.3f\n",max_B,min_B);
#endif
						}
						B[c][index_B] = min_B;
#ifdef DEBUG
						if(row == 58){
							fprintf(fp6,"B[%d][%d] = %3.3f\n",c,index_B,B[c][index_B]);
						}
#endif

#ifdef PRINTMSG
						printf("B[%d][%d] = %3.3f\n",c,index_B,B[c][index_B]);
#endif
					}
				}
			}

			minmax_label10:for(index_c = 0, col = col_row[row][0]; index_c < row_weight_; index_c++, col = col_row[row][index_c]){
	//2. update BETA values from B and F values
#pragma HLS LOOP_TRIPCOUNT min=620 max=620


				if(index_c == 0){
					minmax_label11:for(a = 0; a < Q; a++){


						BETAmn_[row][col][a] = B[mul_gf(H_nb[row][col],a)][index_c + 1];
#ifdef PRINTMSG
						printf("BETAmn_[%d][%d][%d] = B[%d][%d] = %3.3f\n",row,col,a,a,index_c + 1,BETAmn_[row][col][a]);
#endif

#ifdef DEBUG
						fprintf(fp8,"BETAmn_[%d][%d][%d] = B[%d][%d] = %3.3f\n",row,col,a,a,index_c + 1,BETAmn_[row][col][a]);
#endif
					}
				}
				else if(index_c == (row_weight[row] - 1)){
					minmax_label12:for(a = 0; a < Q; a++){

						BETAmn_[row][col][a] = F[mul_gf(H_nb[row][col],a)][index_c - 1];
#ifdef PRINTMSG
						printf("BETAmn_[%d][%d][%d] = F[%d][%d] = %3.3f\n",row,col,a,a,index_c - 1,BETAmn_[row][col][a]);
#endif

#ifdef DEBUG
						fprintf(fp8,"BETAmn_[%d][%d][%d] = F[%d][%d] = %3.3f\n",row,col,a,a,index_c - 1,BETAmn_[row][col][a]);
#endif
					}
				}
				else{

					minmax_label13:for(c = 0; c < Q; c++){


#ifdef PRINTMSG
						printf("c = %d\n",c);
#endif
						min_value = max(F[mul_gf(H_nb[row][col],c)][index_c - 1], B[0][index_c + 1]);
#ifdef PRINTMSG
						printf("F[%d][%d]  = %3.3f B[%d][%d] = %3.3f\n",mul_gf(H_nb[row][col],c),index_c - 1,F[mul_gf(H_nb[row][col],c)][index_c - 1],0,index_c + 1,B[0][index_c + 1]);
						printf("min_value = %3.3f\n",min_value);
#endif

						minmax_label14:for(b = 0; b < Q; b++){



							a = add_gf(b,mul_gf(H_nb[row][col],c));
#ifdef PRINTMSG
							printf("a = %d b = %d\n",a,b);
							printf("F[%d][%d] = %3.3f  B[%d][%d] = %3.3f\n",a,index_c - 1,F[a][index_c - 1],b,index_c + 1,B[b][index_c + 1]);
#endif
							max_value = max(F[a][index_c - 1],B[b][index_c + 1]);
#ifdef PRINTMSG
							printf("max_value = %3.3f\n",max_value);
#endif
							min_value = min(max_value, min_value);
#ifdef PRINTMSG
							printf("min_value = %3.3f\n",min_value);
#endif
						}
						BETAmn_[row][col][c] = min_value;
#ifdef PRINTMSG
						printf("BETAmn_[%d][%d][%d] = B[%d][%d] = %3.3f\n",row,col,a,a,index_c + 1,BETAmn_[row][col][c]);
#endif
					}
				}
			}
		}
		//printf("check node processing finished\n");

#ifdef PRINT_A_B
		for(row = 0; row < M; row++){
			for(col = 0; col < N; col++){
				for(a = 0; a < Q; a++){
					fprintf(fp4, "%3f ", BETAmn_[row][col][a]);
				}
				fprintf(fp4, "\n");
			}
		}
#endif

		minmax_label15:for(col_v = 0; col_v < V; col_v ++){
	// 	 processing: update variable node message ALPHA
#ifdef PRINTMSG_V
			printf("col_v = %d\n",col_v);
#endif
			index_v = 0; 
			row_v = row_col[0][col_v];
			minmax_label16:while(index_v < (col_weight[col_v])){
#pragma HLS LOOP_TRIPCOUNT min=620 max=620 avg=620
				row_v = row_col[index_v][col_v];
#ifdef PRINTMSG_V
				printf("index_v = %d row_v = %d\n",index_v, row_v);
#endif
				minmax_label17:for(a = 0; a < Q; a++){
#ifdef PRINTMSG_V
					printf("a = %d\n",a);
#endif
					temp = 0;
					index_ = 0;
					row_ = row_col[0][col_v];
					minmax_label18:while(index_ < col_weight[col_v]){
#pragma HLS LOOP_TRIPCOUNT min=620 max=620 avg=620
						row_ = row_col[index_][col_v];
						if(index_ != index_v){
#ifdef PRINTMSG_V
							printf("index_ = %d row_ = %d\n",index_, row_);
#endif
#ifdef PRINT_A_B
							fprintf(fp7, "%3f ", BETAmn_[row_][col_v][a]);
#endif
							temp += BETAmn_[row_][col_v][a];
#ifdef PRINTMSG_V
							printf("BETAmn_[%d][%d][%d] = %3f temp = %3f\n",row_, col_v, a, BETAmn_[row_][col_v][a], temp);
#endif
						}
						index_++;
					}
#ifdef PRINT_A_B
					fprintf(fp7, "  ");
#endif
					ALPHA_t[a] = temp + GAMMAn_[a][col_v];
#ifdef PRINTMSG_V
					printf("temp = %3f GAMMAn_[%d][%d] = %3f ALPHA_t[%d] = %3f\n",temp, a, col_v, GAMMAn_[a][col_v] , a, ALPHA_t[a]);
#endif
				}
#ifdef PRINT_A_B
				fprintf(fp7, "\n");
#endif
				min_index = 0;
				minmax_label19:for(a = 1; a < Q; a++){
					if(ALPHA_t[a] < ALPHA_t[min_index]){
						min_index = a;
					}
				}

				minmax_label20:for(a = 0; a < Q; a++){
					ALPHAmn_[row_v][col_v][a] = ALPHA_t[a] - ALPHA_t[min_index];
				}
#ifdef PRINTMSG_V
				printf("ALPHA_t[min] = %3f ALPHA_t[%d] = %3f ALPHAmn_[%d][%d] = %3f\n", ALPHA_t[min_index], a, ALPHA_t[a], row_v, col_v, ALPHAmn_[row_v][col_v][a]);
#endif
				index_v++;
			}
		}
		//printf("variable node processing finished\n");
#ifdef PRINT_A_B
		for(row = 0; row < M; row++){
			for(col = 0; col < N; col++){
				for(a = 0; a < Q; a++){
					fprintf(fp1, "%3f ", ALPHAmn_[row][col][a]);
				}
				fprintf(fp1, "\n");
			}
		}
#endif

		minmax_label21:for(col = 0; col < V; col ++){
	// post processing

			minmax_label22:for( a = 0; a < Q; a++){

				temp = 0;
				index_ = 0;
				row_ = row_col[0][col];
				minmax_label23:while(index_ <col_weight[col] ){
#pragma HLS LOOP_TRIPCOUNT min=620 max=620

					row_ = row_col[index_][col];
#ifdef PRINTMSG_P
					printf("index_ = %d row_ = %d\n",index_, row_);
#endif
					temp += BETAmn_[row_][col][a];
#ifdef PRINTMSG_P
					printf("BETAmn_[%d][%d][%d] = %3.3f temp = %3.3f\n", row_, col, a, BETAmn_[row_][col][a], temp);
#endif
					index_++;
				}
				GAMMAn_post[a][col] = temp + GAMMAn[a][col];
#ifdef DEBUG
				printf("GAMMAn[%d][%d] = %3.3f GAMMAn_post[%d][%d] = %3.3f\n", a, col, GAMMAn_post[a][col], a, col, GAMMAn_[a][col]);
#endif
			}
		}
		//printf("post processing finished\n");
		minmax_label24:for(col = 0; col < V; col ++){	// tentative decoding
			a = 0;
			minmax_label25:for(i = 0; i < Q; i++){
#ifdef DEBUG
				fprintf(fp2,"%3f ",GAMMAn_post[i][col] );
#endif
				if(GAMMAn_post[i][col] < GAMMAn_post[a][col]){
					a = i;
				}

			}
			codeword_sym[col] = a;
#ifdef DEBUG
			fprintf(fp,"%3d ",codeword_sym[col]);
			fprintf(fp2, "\n");
#endif
		}

		error = 0;

		for(i = 0;i < (U);i ++){
			if(encoded_sym[i + U] != codeword_sym[i + U])
				error++;
		}
		if(error == 0){
#ifdef _NOSYNTH_
			for(i = 0; i < Q; i++){
				free(B[i]);
			}
			free(B);

			for(i = 0; i < Q; i++){
				free(F[i]);
			}
			free(F);

			free(ALPHA_t);
			free(codeword_sym);

			for(i = 0; i < Q; i++){
				free(GAMMAn_post[i]);
			}
			free(GAMMAn_post);
#endif
			return 0;

		}

#ifdef DEBUG
		fprintf(fp,"error = %d ",error);
		fprintf(fp, "\n");
#endif		

#ifdef DEBUG
		for(row = 0; row < M; row ++){
			syndrome = 0;
			for(col = 0; col < N; col ++){
				syndrome = add_gf(syndrome,mul_gf(H_nb[row][col],codeword_sym[col]));
			}
			//printf("row = %d syndrome = %d\n",row, syndrome);
		}
#endif		
		iter ++;
	}

#ifdef DEBUG
	fclose(fp);
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);
	fclose(fp5);
	fclose(fp6);
#endif

	for(i = 0; i < Q; i++){
#ifdef _NOSYNTH_
		free(B[i]);
#endif
		}
#ifdef _NOSYNTH_

	free(B);
#endif
	for(i = 0; i < Q; i++){
#ifdef _NOSYNTH_

		free(F[i]);
#endif
	}
#ifdef _NOSYNTH_

	free(F);
#endif
#ifdef _NOSYNTH_

	free(ALPHA_t);
	free(codeword_sym);
#endif
	for(i = 0; i < Q; i++){
#ifdef _NOSYNTH_
		free(GAMMAn_post[i]);
#endif
	}
#ifdef _NOSYNTH_

	free(GAMMAn_post);
#endif

	return error;

}
