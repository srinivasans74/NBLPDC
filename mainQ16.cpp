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

	/*if (a == 0 || b == 0) return 0;
		if (a == 1) return b;
		if (b == 1) return a;
		return expq[(logq[a] + logq[b]) % (Q-1)];
*/
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
// Minmax nonbinary decode function
//===================================
// Refer to "Min-max decoding for nonbinary LDPC codes"
int minmax(unsigned char BETAmn_[U][V][Q], unsigned char ALPHAmn_[U][V][Q], unsigned char GAMMAn_[Q][V])
{
#pragma HLS ARRAY_PARTITION variable=GAMMAn_ complete dim=1
#pragma HLS ARRAY_PARTITION variable=BETAmn_ complete dim=1

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
	double min_F,min_B;
	double max_F,max_B;
	double temp;
	double min_value,max_value;
	double ALPHA_t[Q];
#pragma HLS ARRAY_RESHAPE variable=ALPHA_t complete factor=1
	int codeword_sym[V];
	double B[Q][Y];
	double F[Q][Y];


	for(i = 0; i < Q; i++){
#pragma HLS unroll
		for(j = 0; j < Y; j++){

#pragma HLS unroll

#pragma HLS ARRAY_PARTITION variable=F complete dim=1
#pragma HLS ARRAY_PARTITION variable=B complete dim=1

			F[i][j] = B[i][j] = 0.0;
		}
	}

	for(i = 0; i < Q; i++){
#pragma HLS unroll
		for(j = 0; j < V; j++){
#pragma HLS unroll
#pragma HLS ARRAY_PARTITION variable=GAMMAn_post complete dim=1
			GAMMAn_post[i][j] = 0.0;
		}
	}

	minmax_label1:for(row = 0; row < U; row++)// initialize variable node message ALPHA with channel info GAMMA
	{
#pragma HLS unroll
		col = col_row[row][0];
		minmax_label23:for(index = 0; index <2 ; index++)
		{
#pragma HLS unroll
			col = col_row[row][index];
			for(a = 0; a < Q; a++){
#pragma HLS unroll
				ALPHAmn_[row][col][a] = GAMMAn_[a][col];
				printf(" ALPHAMN[%d][%d][%d] %hhu \n ",row, col,a,ALPHAmn_[row][col][a]);
			}
		}
	}

	criticalpath1:
	for(col = 0; col < V; col ++){// tentative decoding
		a = 0;
#pragma HLS unroll
		for(i = 0; i < Q; i++){
#pragma HLS unroll
			if(GAMMAn_[i][col] < GAMMAn_[a][col])
			{
				a = i;
			}
		}
		codeword_sym[col] = a;
		printf("%dcodeword symbol \n",codeword_sym[col]);

	}

for (i=0;i<V;i++)
{
	printf("encoded symbol %d \n",encoded_sym[V]);

}

	error = 0;
	for(i = 0;i < (U);i ++){
#pragma HLS unroll
#pragma HLS ARRAY_PARTITION variable=GAMMAn_ complete dim=1
		if(encoded_sym[i + U] != codeword_sym[i + U])
			error++;
	printf("updated error %d",error);
	}


	if(error == 0) {

	//return 0;
	}

	minmax_label2:while(iter < MAX_ITERATION){
#pragma HLS pipeline II=20
		minmax_label3:for(row = 0; row <U; row++){	// check node processing: update check node message BETA
			row_weight_ = 2;
			minmax_label4:for(index_c = 0, col = col_row[row][0]; index_c < row_weight_; index_c++, col = col_row[row][index_c]){
				// 1. Initialize F and B value for each row
				if(index_c == 0){
					minmax_label5:for(a = 0; a < Q; a++){
#pragma HLS unroll
						F[a][0] = ALPHAmn_[row][col][mul_gf(inv_gf(H_nb[row][col]),a)];
						last_col = col_row[row][row_weight_ - 1];
						B[a][row_weight_ - 1] = ALPHAmn_[row][last_col][mul_gf(inv_gf(H_nb[row][last_col]),a)];
					}
				}
				else{

					index_p = index_c - 1;
					col_p = col_row[row][index_p];
					minmax_label323:for(c = 0; c < Q; c++){
#pragma HLS unroll
						min_F = max(F[c][index_c - 1],ALPHAmn_[row][col][0]);
						minmax_label7:for(b = 0; b < Q; b++){
#pragma HLS unroll
							a = sub_gf(c,mul_gf(H_nb[row][col],b));
							max_F = max(F[a][index_c - 1], ALPHAmn_[row][col][b]);
							min_F = min(min_F, max_F);
				}
						F[c][index_c] = min_F;
					}
					index_B = row_weight_ - index_c - 1;
					col_B = col_row[row][index_B];
					index_p = index_B + 1;
					col_p = col_row[row][index_p];
					minmax_label8:for(c = 0; c < Q; c++){

#pragma HLS unroll
						min_B = max(ALPHAmn_[row][col_B][0],B[c][index_B + 1]);
						minmax_label9:for(b = 0; b < Q; b++){
#pragma HLS unroll
							a = sub_gf(c,mul_gf(H_nb[row][col_B],b));
							max_B = max(B[a][index_B + 1], ALPHAmn_[row][col_B][b]);
							min_B = min(min_B, max_B);
						}
						B[c][index_B] = min_B;

					}
				}
			}
			minmax_label10:for(index_c = 0, col = col_row[row][0]; index_c < row_weight_; index_c++, col = col_row[row][index_c]){
	//2. update BETA values from B and F values
				if(index_c == 0){
					minmax_label11:for(a = 0; a < Q; a++){
#pragma HLS unroll
						BETAmn_[row][col][a] = B[mul_gf(H_nb[row][col],a)][index_c + 1];
					}
				}
				else if(index_c == (row_weight[row] - 1)){
					minmax_label12:for(a = 0; a < Q; a++){
#pragma HLS unroll
						BETAmn_[row][col][a] = F[mul_gf(H_nb[row][col],a)][index_c - 1];

					}
				}
				else{

					minmax_label13:for(c = 0; c < Q; c++){
#pragma HLS unroll
						min_value = max(F[mul_gf(H_nb[row][col],c)][index_c - 1], B[0][index_c + 1]);
						minmax_label14:for(b = 0; b < Q; b++){
#pragma HLS unroll
							a = add_gf(b,mul_gf(H_nb[row][col],c));
							max_value = max(F[a][index_c - 1],B[b][index_c + 1]);
							min_value = min(max_value, min_value);
						}
						BETAmn_[row][col][c] = min_value;

					}
				}
			}
		}
		//printf("check node processing finished\n");
		minmax_label15:for(col_v = 0; col_v < V; col_v ++){
	// 	 processing: update variable node message ALPHA
#pragma HLS ARRAY_PARTITION variable=row_col complete dim=1
			index_v = 0; 
			row_v = row_col[0][col_v];
			minmax_label132:while(index_v < (1)){
#pragma HLS unroll
				row_v = row_col[index_v][col_v];
		minmax_label17:for(a = 0; a < Q; a++){
		temp = 0;
					index_ = 0;
					row_ = row_col[0][col_v];
					minmax_label18:while(index_<1){
						row_ = row_col[index_][col_v];
						if(index_ != index_v){
							temp += BETAmn_[row_][col_v][a];
						}
						index_++;
					}
					ALPHA_t[a] = temp + GAMMAn_[a][col_v];
				}
				min_index = 0;
				minmax_label19:for(a = 1; a < Q; a++){
#pragma HLS unroll
					if(ALPHA_t[a] < ALPHA_t[min_index]){
						min_index = a;
					}
				}
				minmax_label20:for(a = 0; a < Q; a++){
#pragma HLS unroll
					ALPHAmn_[row_v][col_v][a] = ALPHA_t[a] - ALPHA_t[min_index];
	printf("ALPHAMAN[%d][%d][%d]=%d\n",row_v,col_v,a,ALPHAmn_[row_v][col_v][a]);
				}
				index_v++;
			}
		}



		//printf("variable node processing finished\n");
		//FILE *fp66 = fopen("C:/Users/mtech/Desktop/Q16/data/Gamman_post", "w+");
		/*minmax_label21:for(col = 0; col < V; col ++){ // post processing
			minmax_label22:for( a = 0; a < Q; a++){
				temp = 0;
				index_ = 0;
				row_ = row_col[0][col];
				minmax_label723:while(index_ <1 ){
					row_ = row_col[index_][col];
					temp += BETAmn_[row_][col][a];
					index_++;
				}
				GAMMAn_post[a][col] = temp + GAMMAn[a][col];
			}
			//fprintf(fp66, "%d ",GAMMAn_post[a][col]);
		}
		//fclose(fp66);
		//printf("post processing finished\n");
		minmax_label24:for(col = 0; col < V; col ++){// tentative decoding
			a = 0;
#pragma HLS unroll
			minmax_label25:for(i = 0; i < Q; i++){
				if(GAMMAn_post[i][col] < GAMMAn_post[a][col]){
#pragma HLS unroll
					a = i;
				}
			}
			codeword_sym[col] = a;
		}
		error = 0;
		for(i = 0;i < (U);i ++){
			if(encoded_sym[i + U] != codeword_sym[i + U])
#pragma HLS unroll
				error++;
		}
		if(error == 0){
			return 0;
		}*/
		iter ++;
	}
	return error;
}
