//***************************************************************************
//	Nonbinary LDPC decoder simulation
//	Minmax soft decoding
//	Hao Shen
//	Rice University
//***************************************************************************
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "nbldpcq16.h"

#include <time.h>        // for Windows APIs
clock_t t1, t2;           // ticks


int				N;
int				M;
int				dc, dv;
int				**H_qc;
int				**G_qc;
int				H_nb[U][V];
int				G_nb[U][V];
int				row_weight[U];
int				col_weight[V];
int				row_col[X][V];
int				col_row[U][Y];
unsigned char	BETAmn[U][V][Q];
unsigned char	ALPHAmn[U][V][Q];
unsigned char	GAMMAn[Q][V];
int				encoded_sym[V];
unsigned char	GAMMAn_post[Q][V];

//===================================
// Allocate 2D array of int
//===================================
int **malloc2Dint(int a, int b) // allocates array[a][b]
{
	int i;
	int **pp = (int **)malloc(sizeof(int*) * a);
	int *p = (int *)malloc(sizeof(int) * a * b);
	if (pp == NULL || p == NULL) exit(-1);
	for (i = 0; i < a; i++) {
		pp[i] = p + b * i;
	}
	return pp;
}



double ***malloc3Ddouble(int a, int b, int c)
{
	int i, j;
	double*** ppp = (double ***)malloc(sizeof(double **) * M);
	if (ppp == NULL) exit(-1);
	for (i = 0; i < a; i++)
	{
		ppp[i] = (double **)malloc(sizeof(double *) * b);
		for (j = 0; j < b; j++)
			ppp[i][j] = (double *)malloc(sizeof(double) * c);
	}
	return ppp;
}

unsigned char ***malloc3Dunsigned_char(int a, int b, int c)
{
	int i, j;
	unsigned char*** ppp = (unsigned char ***)malloc(sizeof(unsigned char **) * M);
	if (ppp == NULL) exit(-1);
	for (i = 0; i < a; i++)
	{
		ppp[i] = (unsigned char **)malloc(sizeof(unsigned char *) * b);
		for (j = 0; j < b; j++)
			ppp[i][j] = (unsigned char *)malloc(sizeof(unsigned char) * c);
	}
	return ppp;
}

//===================================
// Generate Parity Check Matrix
// Non-binary H = (310,620) degree = (3,6)
// Refer to "Construction of Non-Binary Quasi-cyclic LDPC Codes by Arrays and Array Dispersions", Bo Zhou etc.
//===================================
void qc_gen(void)
{
	int i, j, k, x, x_index, m, n;
	int** W;
	int** H_disp;
	int* H_disp_loc_orig;
	int* H_disp_loc_current;
	int** Hj_up;
	int** Hj_low;
	int* H_disp_mask_loc_orig;
	int* H_disp_mask_loc_current;
	int** H_tmp;
	int* one_row;
	int** col_ext_matrix;
	int** H;
	int K, L, S, T;
	int row_w;
	int col_w;
	int** Hj;
	FILE *fp = fopen("C:/Users/mtech/Desktop/Q16/data/qc.txt", "w+");
	FILE *fp1 = fopen("C:/Users/mtech/Desktop/Q16/data/h_nb.txt", "w+");
	FILE *fp_row_col = fopen("C:/Users/mtech/Desktop/Q16/data/row_col.txt", "w+");
	FILE *fp_col_row = fopen("C:/Users/mtech/Desktop/Q16/data/col_row.txt", "w+");
	FILE *fprow_col1 = fopen("C:/Users/mtech/Desktop/Q16/data/row_col1.txt", "w+");
	FILE *fpcol_row1 = fopen("C:/Users/mtech/Desktop/Q16/data/col_row1.txt", "w+");
	FILE *fprow_weight = fopen("C:/Users/mtech/Desktop/Q16/data/row_weight.txt", "w+");




	if (Q > 8) {
		K = 2;
		L = 2;
		S = 2;
		T = 5;
#ifdef ERROR_CHECK
		printf("check Q>8\n");
#endif
	}

	else {
		K = 2;
		L = 1;
		S = 0;
		T = 1;
#ifdef ERROR_CHECK
		printf("check Q<=8\n");
#endif
	}

	M = (Q - 1) * T * L;
	N = (Q - 1) * T * L * K;

	dc = T - S;
	dv = K * (T - S);

	W = malloc2Dint((Q - 1), (Q - 1));
	Hj = malloc2Dint(T * (Q - 1),T * (Q - 1));
	H_tmp = malloc2Dint((Q - 1) * (Q - 1),(Q - 1) * (Q - 1));
	Hj_up = malloc2Dint(T * (Q - 1),T * (Q - 1));
	Hj_low = malloc2Dint(T * (Q - 1),T * (Q - 1));
	one_row = (int*)malloc(sizeof(int) * (Q - 1));
	col_ext_matrix = malloc2Dint((Q - 1) * (Q - 1), (Q - 1));
	H_disp = malloc2Dint(L * T * (Q - 1), L * T * (Q - 1));
	H_disp_loc_orig = (int*)malloc(sizeof(int) * L);
	H_disp_loc_current = (int*)malloc(sizeof(int) * L);
	H_disp_mask_loc_orig = (int*)malloc(sizeof(int) * L * T);
	H_disp_mask_loc_current = (int*)malloc(sizeof(int) * L * T);
	H = malloc2Dint((Q - 1) * T, K * (Q - 1) * T);



#ifdef ERROR_CHECK
	printf("Qc_gen() alloc check\n");
#endif

	for (i = 0; i < M; i++) {
		row_weight[i] = 0;
	}

#ifdef ERROR_CHECK
	printf("Qc_gen() row_weight init check\n");
#endif

	for (i = 0; i < N; i++) {
		col_weight[i] = 0;
	}

#ifdef ERROR_CHECK
	printf("Qc_gen() col_weight init check\n");
#endif

	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++) {
			for (k = 0; k < Q; k++) {
				BETAmn[i][j][k] = 0;
				ALPHAmn[i][j][k] = 0;
			}
		}
	}

#ifdef ERROR_CHECK
	printf("Qc_gen() alpha beta init check\n");
#endif

	for (i = 0; i < Q; i++) {
		for (j = 0; j < N; j++) {
			GAMMAn[i][j] = 0;
		}
	}

#ifdef ERROR_CHECK
	printf("Qc_gen() gamma init check\n");
#endif

	for (i = 0; i < dc; i++) {
		for (j = 0; j < N; j++) {
			row_col[i][j] = 0;
		}
	}

#ifdef ERROR_CHECK
	printf("Qc_gen() row_col init check\n");
#endif

	for (i = 0; i < M; i++) {
		for (j = 0; j < dv; j++) {
			col_row[i][j] = 0;
		}
	}

#ifdef ERROR_CHECK
	printf("Qc_gen() col_row init check\n");
#endif

	for (i = 0; i < (Q - 1) * T; i++) {
		for (j = 0; j < K * (Q - 1) * T; j++) {
			H[i][j] = 0;
		}
	}

#ifdef ERROR_CHECK
	printf("Qc_gen() H init check\n");
#endif

	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++) {
			H_nb[i][j] = 0;
		}
	}

#ifdef ERROR_CHECK
	printf("Qc_gen() H_nb init check\n");
#endif

	for (i = 0; i < N - M; i++) {
		for (j = 0; j < N; j++) {
			G_nb[i][j] = 0;
		}
	}

#ifdef ERROR_CHECK
	printf("Qc_gen() G_nb init check\n");
#endif

	for (i = 0; i < (Q - 1); i++) {
		for (j = 0; j < (Q - 1); j++) {
			W[i][j] = 0;
		}
	}

#ifdef ERROR_CHECK
	printf("Qc_gen() W init check\n");
#endif

	for (i = 0; i < (Q - 1) * (Q - 1); i++) {
		for (j = 0; j < (Q - 1) * (Q - 1); j++) {
			H_tmp[i][j] = 0;
		}
	}

#ifdef ERROR_CHECK
	printf("Qc_gen() H_tmp init check\n");
#endif

	for (i = 0; i < T * (Q - 1); i++) {
		for (j = 0; j < T * (Q - 1); j++) {
			Hj_up[i][j] = Hj_low[i][j] = Hj[i][j] = 0;
		}
	}

#ifdef ERROR_CHECK
	printf("Qc_gen() Hj init check\n");
#endif

	for (i = 0; i < L * T * (Q - 1); i++) {
		for (j = 0; j < L * T * (Q - 1); j++) {
			H_disp[i][j] = 0;
		}
	}

#ifdef ERROR_CHECK
	printf("Qc_gen() H_disp init check\n");
#endif

	for (i = 0; i < L; i++) {
		if (i == 0) {
			H_disp_loc_orig[i] = 1;
		}
		else if (i == (L - 1)) {
			H_disp_loc_orig[i] = -1;
		}
		else {
			H_disp_loc_orig[i] = 0;
		}
	}

#ifdef ERROR_CHECK
	printf("Qc_gen() H_disp_loc init check\n");
#endif

	for (i = 0; i < L * T; i++) {
#ifdef ERROR_CHECK
		printf("i=%d, L=%d, T=%d\n", i, L, T);
#endif
		if (i == 0) {
			H_disp_mask_loc_orig[i] = 1;
		}
		else if (i >= (T + S + 1)) {
			H_disp_mask_loc_orig[i] = 1;
		}
		else {
			H_disp_mask_loc_orig[i] = 0;
		}
		fprintf(fp, "%2d ", H_disp_mask_loc_orig[i]);
	}

#ifdef ERROR_CHECK
	printf("Qc_gen() H_disp_mask init check\n");
#endif

#ifdef ERROR_CHECK
	printf("Qc_gen() variable init check\n");
#endif

	fprintf(fp, "\n\n\n\n");

	for (i = 0; i < (Q - 1); i++) {
		one_row[i] = expq[i];
	}

	for (i = 0; i < (Q - 1); i++) {
		for (j = 0; j < (Q - 1); j++) {
			W[i][j] = one_row[((j - i) + (Q - 1)) % (Q - 1)] - 1;
			fprintf(fp, "%2d ", W[i][j]);
		}
		fprintf(fp, "\n");
	}

	fprintf(fp, "\n\n\n\n");


	for (i = 0; i < (Q - 1) * (Q - 1); i++) {
		x = i;
		x_index = 0;
		while (x >= (Q - 1)) {
			x -= (Q - 1);
			x_index++;
		}
		for (j = 0; j < (Q - 1); j++) {
			col_ext_matrix[i][j] = mul_gf(W[x_index][j], one_row[x]);
			fprintf(fp, "%2d ", col_ext_matrix[i][j]);
		}
		fprintf(fp, "\n");
	}

	fprintf(fp, "\n");
	fprintf(fp, "\n");
	fprintf(fp, "\n");
	fprintf(fp, "\n");

	fprintf(fp, "hello\n");

	for (i = 0; i < (Q - 1) * (Q - 1); i++) {
		for (j = 0; j < (Q - 1); j++) {
			k = 0;
			while (one_row[k] != col_ext_matrix[i][j]) {
				k++;
			}
			H_tmp[i][j * (Q - 1) + k] = col_ext_matrix[i][j];
		}
	}

	for (i = 0; i < (Q - 1) * (Q - 1); i++) {
		for (j = 0; j < (Q - 1) * (Q - 1); j++) {
			fprintf(fp, "%2d ", H_tmp[i][j]);
		}
		fprintf(fp, "\n");
	}

	fprintf(fp, "\n\n\n\n");

	for (i = 0; i < (Q - 1) * T; i++) {
		for (j = 0; j < K * (Q - 1) * T; j++) {
			H[i][j] = H_tmp[(Q - 1) * (Q - 1) - (Q - 1) * T + i][j];
			fprintf(fp, "%2d ", H[i][j]);
		}
		fprintf(fp, "\n");
	}

	fprintf(fp, "\n\n\n\n");

	for (k = 0; k < K; k++) {
		fprintf(fp, "Hj\n");
		for (i = 0; i < T * (Q - 1); i++) {
			for (j = 0; j < T * (Q - 1); j++) {
				Hj[i][j] = H[i][j + (Q - 1) * k * T];
				fprintf(fp, "%2d ", Hj[i][j]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n\n\n\n");

		for (i = 0; i < T * (Q - 1); i++) {
			for (j = 0; j < T * (Q - 1); j++) {
				if (j > (i + (Q - 1))) {
					Hj_up[i][j] = Hj[i][j];
					Hj_low[i][j] = 0;
				}
				else {
					Hj_low[i][j] = Hj[i][j];
					Hj_up[i][j] = 0;
				}
			}
		}

		fprintf(fp, "Hj_low\n");
		for (i = 0; i < T * (Q - 1); i++) {
			for (j = 0; j < T * (Q - 1); j++) {
				fprintf(fp, "%2d ", Hj_low[i][j]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n\n\n\n");
		fprintf(fp, "Hj_up\n");
		for (i = 0; i < T * (Q - 1); i++) {
			for (j = 0; j < T * (Q - 1); j++) {
				fprintf(fp, "%2d ", Hj_up[i][j]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n\n\n\n");


		fprintf(fp, "\n\n\n\n");


		for (m = 0; m < L; m++) {
			for (n = 0; n < L; n++) {
				H_disp_loc_current[n] = H_disp_loc_orig[((n - m) + L) % L];
				fprintf(fp, "%2d ", H_disp_loc_current[n]);
			}
			for (n = 0; n < L; n++) {
				if (H_disp_loc_current[n] == 1) {
					for (i = 0 + m * T * (Q - 1); i < T * (Q - 1) + m * T * (Q - 1); i++) {
						for (j = 0 + n * T * (Q - 1); j < T * (Q - 1) + n * T * (Q - 1); j++) {
							H_disp[i][j] = Hj_low[i - m * T * (Q - 1)][j - n * T * (Q - 1)];
						}
					}
				}
				else if (H_disp_loc_current[n] == -1) {
					for (i = 0 + m * T * (Q - 1); i < T * (Q - 1) + m * T * (Q - 1); i++) {
						for (j = 0 + n * T * (Q - 1); j < T * (Q - 1) + n * T * (Q - 1); j++) {
							H_disp[i][j] = Hj_up[i - m * T * (Q - 1)][j - n * T * (Q - 1)];
						}
					}
				}
				else {
					for (i = 0 + m * T * (Q - 1); i < T * (Q - 1) + m * T * (Q - 1); i++) {
						for (j = 0 + n * T * (Q - 1); j < T * (Q - 1) + n * T * (Q - 1); j++) {
							H_disp[i][j] = 0;
						}
					}
				}
			}
		}

		for (i = 0; i < L * T * (Q - 1); i++) {
			for (j = 0; j < L * T * (Q - 1); j++) {
				fprintf(fp, "%2d ", H_disp[i][j]);
			}
			fprintf(fp, "\n");
		}

		for (m = 0; m < L * T; m++) {
			for (n = 0; n < L * T; n++) {
				H_disp_mask_loc_current[n] = H_disp_mask_loc_orig[((n - m) + L * T) % (L * T)];
				if (H_disp_mask_loc_current[n] == 0) {
					for (i = m * (Q - 1); i < (m + 1) * (Q - 1); i++) {
						for (j = n * (Q - 1); j < (n + 1) * (Q - 1); j++) {
							H_disp[i][j] = 0;
						}
					}
				}
			}
		}

		for (i = 0; i < L * T * (Q - 1); i++) {
			for (j = 0; j < L * T * (Q - 1); j++) {
				H_nb[i][j + k * L * T * (Q - 1)] = H_disp[i][j];
			}
		}
	}

	fprintf(fp, "\n\n\n\n");
	fprintf(fp, "\n\n\n\n");

	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++) {
			fprintf(fp1, "%d", H_nb[i][j]);
		}
		fprintf(fp1, "\n");
	}

	for (i = 0; i < M; i++) {
		row_w = 0;
		for (j = 0; j < N; j++) {
			if (H_nb[i][j] != 0) {
				row_w++;
			}
		}
		fprintf(fp, "row_w[%d] = %d\n", i, row_w);
	}

	for (j = 0; j < N; j++) {
		col_w = 0;
		for (i = 0; i < M; i++) {
			if (H_nb[i][j] != 0) {
				col_w++;
			}
		}
		fprintf(fp, "col_w[%d] = %d\n", j, col_w);
	}

	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++) {
			if (H_nb[i][j] != 0) {
				row_weight[i]++;
			}
		}
		fprintf(fp, "row_weight[%d] = %d\n", i, row_weight[i]);
		fprintf(fprow_weight, "row_weight[%d] = %d\n", i, row_weight[i]);
	}

	for (j = 0; j < N; j++) {
		for (i = 0; i < M; i++) {
			if (H_nb[i][j] != 0) {
				col_weight[j]++;
			}
		}
		fprintf(fp, "col_weight[%d] = %d\n", j, col_weight[j]);
	}

	for (i = 0; i < M; i++) {
		for (j = 0; j < row_weight[i]; j++) {
			for (k = 0; k < N; k++) {
				if (H_nb[i][k] != 0) {
					col_row[i][j] = k;
					fprintf(fp_col_row, "%d ", col_row[i][j]);
					fprintf(fpcol_row1, "%d ", col_row[i][j]);

					j++;
				}
			}
		}
		fprintf(fp, "\n");
		fprintf(fpcol_row1, "\n");

	}

	for (j = 0; j < N; j++) {
		for (i = 0; i < col_weight[j]; i++) {
			for (k = 0; k < M; k++) {
				if (H_nb[k][j] != 0) {
					row_col[i][j] = k;
					fprintf(fp_row_col, "%d ", row_col[i][j]);
					fprintf(fprow_col1, "%d ", row_col[i][j]);

					i++;
				}
			}
		}
		fprintf(fp_row_col, "\n");
		fprintf(fprow_col1, "\n");

	}



	fprintf(fp, "\n\n\n\n");
	fclose(fprow_weight);
	fclose(fp_row_col);
	fclose(fp_col_row);
	fclose(fprow_col1);
	fclose(fpcol_row1);
	fclose(fp);
	fclose(fp1);
}

//===================================
// Generate Random Parity Check Matrix
// Refer to http://www.inference.phy.cam.ac.uk/mackay/codes/data.html
//===================================
//#endif

//===================================
// H to G transform
// Using Gauss Jordan Elimination
// binary
//===================================
void h2g_nb(int H[U][V], int n)
{
	int i, j, k, p, q, maxrow, scale;
	int tmp;
	int H_[U][V];
	FILE *fp, *fpp;

#ifdef QC
	fp = fopen("data/qc.H'.nb.txt", "w+");
	fpp = fopen("data/qc.G'.nb.txt", "w+");
#else
	if (Q == 16)
	{
		fp = fopen("matrix/204.H'.nb.txt", "w+");
		fpp = fopen("matrix/204.G'.nb.txt", "w+");
	}
	else if (Q == 4)
	{
		fp = fopen("matrix/408.H'.nb.txt", "w+");
		fpp = fopen("matrix/408.G'.nb.txt", "w+");
	}
	else if (Q == 2)
	{
		fp = fopen("matrix/816.H'.nb.txt", "w+");
		fpp = fopen("matrix/816.G'.nb.txt", "w+");
	}
#endif


	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++) {
			H_[i][j] = H[i][j];
		}
	}

	for (j = 0; j < M; j++) {
		maxrow = j;
		for (i = j; i < M; i++) {	//find max entry in one colum
			if ((H_[i][j]) > (H_[maxrow][j]))
				maxrow = i;
		}

		for (k = j; k < n; k++) {	//swap the row with jth row
			tmp = H_[j][k];
			H_[j][k] = H_[maxrow][k];
			H_[maxrow][k] = tmp;
		}

		scale = inv_gf(H_[j][j]);
		//printf("H_[j][j] is %d and scale is %d\n",H_[j][j],scale);
		for (k = j; k < N; k++) {
			H_[j][k] = mul_gf(H_[j][k], scale);
		}

		for (q = 0; q < M; q++) {		//eliminate other entry within this colum
			if ((q != j) && (H_[q][j] != 0)) {
				tmp = H_[q][j];
				for (p = j; p < n; p++) {
					H_[q][p] = sub_gf(H_[q][p], mul_gf(H_[j][p], tmp));
				}
			}
		}
	}
	printf("eliminate finishes\n");

	for (i = 0; i < M; i++) {
		for (j = 0; j < n; j++) {
			fprintf(fp, "%4d ", H_[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	printf("H' written to file\n");

	for (i = 0; i < M; i++) //build identity matrix in G
		for (j = M; j < n; j++)
			if (j == (i + M)) G_nb[i][j] = 1;

	for (i = 0; i < M; i++)
		for (j = 0; j < M; j++)
			G_nb[i][j] = H_[j][i + M];

	for (i = 0; i < M; i++) {
		for (j = 0; j < n; j++) {
			fprintf(fpp, "%4d ", G_nb[i][j]);
		}
		fprintf(fpp, "\n");
	}
	fclose(fpp);

	printf("G written to file\n");
}

//===================================
// Ramdom info data generation
//===================================
void info_gen(int info_bin[])// radom number generation
{
	/*#ifdef DEBUG
		FILE *fpp = fopen("data/info.txt","w");
	#endif
		srand(time(0));
		int i ;
		for (i = 0 ; i <(N-M)*LOG2Q ; i++)
		{
			info_bin [i] = (rand())%2;
			//info_bin[i] = 1;
	#ifdef DEBUG
			fprintf(fpp, "%d\n",info_bin[i]);
	#endif
			printf("info_bin[%d]=%d\n", i, info_bin[i]);
		}
	#ifdef DEBUG
		fclose(fpp);
	#endif*/


	//======================================
	//read pre generated file

	FILE *file;

	if (Q == 4)
		file = fopen("C:/Users/mtech/Desktop/Q16/data/random_data/Q4/info.txt", "r");
	else if (Q == 8)
		file = fopen("C:/Users/mtech/Desktop/Q16/data/random_data/Q8/info.txt", "r");
	else if (Q == 16)
		file = fopen("C:/Users/mtech/Desktop/Q16/data/random_data/Q16/info.txt", "r");
	else if (Q == 32)
		file = fopen("C:/Users/mtech/Desktop/Q16/data/random_data/Q32/info.txt", "r");
	else if (Q == 64)
		file = fopen("C:/Users/mtech/Desktop/Q16/data/random_data/Q64/info.txt", "r");
	else if (Q == 128)
		file = fopen("C:/Users/mtech/Desktop/Q16/data/random_data/Q128/info.txt", "r");
	else if (Q == 256)
		file = fopen("C:/Users/mtech/Desktop/Q16/data/random_data/Q256/info.txt", "r");

	char line[256];
	int d = 0;
	int i;

	while (fgets(line, sizeof(line), file)) {
		sscanf(line, "%d", &i);
		//printf("%d\n", i);
		info_bin[d] = i;
		d++;
	}
	/*for (i = 0 ; i <(N-M)*LOG2Q ; i++){
		printf("info_bin[%d]=%d\n", i, info_bin[i]);
	}*/

	fclose(file);
}

//===================================
// Group bits into symbols
//===================================
void grouping(int info_bin[], int info_sym[])
{
	int i, j;
	int result, power;
#ifdef DEBUG
	FILE *fp = fopen("data/info_sym.nb.txt", "w+");
#endif
	for (i = 0; i < (N - M); i++)
	{
		result = 0;
		power = 1;
		for (j = LOG2Q - 1; j >= 0; j--)
		{
			if (info_bin[i  * LOG2Q + j] == 1) { result += power; }
			power *= 2;
		}
		info_sym[i] = result;
#ifdef DEBUG	
		fprintf(fp, "%d\n", info_sym[i]);
#endif
	}
#ifdef DEBUG
	fclose(fp);
#endif
}

//===================================
// Encoder nonbinary
//===================================
void encoder_nb(int info_sym[], int encoded_sym[])
{
	int i, j, xHmn, sum;
#ifdef DEBUG
	FILE *fp = fopen("data/encoded_sym.nb.txt", "w+");
#endif

	for (i = 0; i < N; i++) {
		sum = 0;
		for (j = 0; j < (N - M); j++) {
			xHmn = mul_gf(info_sym[j], G_nb[j][i]);
			sum = add_gf(sum, xHmn);
		}
		encoded_sym[i] = sum;
#ifdef DEBUG
		fprintf(fp, "%d\n", encoded_sym[i]);
#endif  
	}

#ifdef DEBUG
	fclose(fp);
#endif 
}

//===================================
// Ungroup symbols into bits
//===================================
void ungrouping(int info_sym[], int info_bin[])
{
	int i, j;
	int decimal, remain;
#ifdef DEBUG
	FILE *fp = fopen("data/info_bin_coded.txt", "w+");
#endif
	for (i = 0; i < N; i++)
	{
		decimal = info_sym[i];
		for (j = (LOG2Q - 1); j >= 0; j--)
		{
			remain = decimal % 2;
			decimal = decimal / 2;
			info_bin[i*LOG2Q + j] = remain;
		}
	}
#ifdef DEBUG
	for (i = 0; i < LOG2Q*N; i++) {
		fprintf(fp, "%d \n", info_bin[i]);
	}
	fclose(fp);
#endif
}

//===================================
// BPSK modulation
//===================================
void bpsk_modulation(int code[], double trans[])
{
	int i;
#ifdef DEBUG
	FILE *fp = fopen("data/trans.txt", "w");
#endif	
	for (i = 0; i < LOG2Q * N; i++) {
		if (code[i] == 1)
			trans[i] = 1;
		else
			trans[i] = -1;
#ifdef DEBUG
		fprintf(fp, "%lf\n", trans[i]);
#endif
	}
#ifdef DEBUG
	fclose(fp);
#endif
}

//===================================
// AWGN Channel
//===================================
void awgn(double *trans, double *rec, double sigma)
{
	/*double u1,u2,s,noise,randmum;
	int i;
	#ifdef DEBUG
		FILE *fp = fopen("data/recv.txt","w");
		printf("sigma is %f\n",sigma);
	#endif
	for (i=0; i< LOG2Q * N; i++)
	{
		do
		{
			randmum = (double)(rand())/RAND_MAX;
			u1 = randmum * 2.0 - 1.0;
			randmum = (double)(rand())/RAND_MAX;
			u2 = randmum * 2.0 - 1.0;
			s = u1 * u1 + u2 * u2;
		} while( s >= 1) ;
		noise = u1 * sqrt( (-2.0 * log(s))/s );
		#ifdef NONOISE
			rec[i] = trans[i];
		#else
			rec[i] = trans[i] + noise * sigma;
		#endif
		#ifdef DEBUG
			fprintf(fp,"%f \n",rec[i]);
		#endif
	}
	#ifdef DEBUG
		printf("awgn finishes\n");
		fclose(fp);
	#endif*/


	//======================================
	//read pre generated file

	FILE *file;

	if (Q == 4)
		file = fopen("C:/Users/mtech/Desktop/Q16/data/random_data/Q4/recv.txt", "r");
	else if (Q == 8)
		file = fopen("C:/Users/mtech/Desktop/Q16/data/random_data/Q8_extra/recv.txt", "r");
	else if (Q == 16)
		file = fopen("C:/Users/mtech/Desktop/Q16/data/random_data/Q16/recv.txt", "r");
	else if (Q == 32)
		file = fopen("C:/Users/mtech/Desktop/Q16/data/random_data/Q32/recv.txt", "r");
	else if (Q == 64)
		file = fopen("C:/Users/mtech/Desktop/Q16/data/random_data/Q64/recv.txt", "r");
	else if (Q == 128)
		file = fopen("C:/Users/mtech/Desktop/Q16/data/random_data/Q128/recv.txt", "r");
	else if (Q == 256)
		file = fopen("C:/Users/mtech/Desktop/Q16/data/random_data/Q256/recv.txt", "r");

	char line[256];
	int d = 0;
	double i;

	while (fgets(line, sizeof(line), file)) {
		sscanf(line, "%lf", &i);
		//printf("%lf\n", i);
		rec[d] = i;
		d++;
	}
	/* for (int x = 0 ; x <LOG2Q * N ; x++){
		printf("rec[%d]=%lf\n", x, rec[x]);
	} */

	fclose(file);
}

//===================================
// Find the most likely symbols
//===================================
void hard_decision(double recv[], double sigma, int hard[])
{
	int i;
#ifdef DEBUG
	FILE * fp = fopen("data/llr_fp.txt", "w");
#endif
	for (i = 0; i < LOG2Q * N; i++)
	{
		if (recv[i] > 0)	hard[i] = 1;
		else hard[i] = 0;
#ifdef DEBUG
		fprintf(fp, "recv[%d] = %lf, hard = %d\n", i, recv[i], hard[i]);
#endif
	}
#ifdef DEBUG
	fclose(fp);
#endif
}

//===================================
// calculate symbol LLRs
//===================================
void sym_llr(double* llr_, unsigned char BETAmn_[U][V][Q], unsigned char GAMMAn_[Q][V], int* hard, double* recv, double sigma)
{
	int i, j, k, l, m;
	int dev, remain;
	double temp;

#ifdef DEBUG
	//FILE *fp = fopen("data/symbol_llr.txt", "w");
#endif

	// initialize channel information
	for (i = 0; i < M; i++) {
		for (m = 0, j = col_row[i][0]; m < row_weight[i]; m++, j = col_row[i][m]) {
			for (k = 0; k < Q; k++) {
				dev = k;
				temp = 0;
				for (l = (LOG2Q - 1); l >= 0; l--) {
					remain = dev % 2;
					dev = dev / 2;

					int Log2Q = LOG2Q;
					if (remain == 1 && hard[j * Log2Q + l] == 0) {
						//BETAmn_[i][j][k] -= 2*recv[j * LOG2Q + l]/(sigma * sigm    a);
						temp -= 2 * recv[j * Log2Q + l] / (sigma * sigma);
					}
					else if (remain == 0 && hard[j*Log2Q + l] == 1) {
						//BETAmn_[i][j][k] += 2 * recv[j * LOG2Q + l] / (sigma * sigma);
						temp += 2 * recv[j * Log2Q + l] / (sigma * sigma);
					}
				}
				//GAMMAn[k][j] = BETAmn_[i][j][k];

				GAMMAn[k][j] = (unsigned char)(temp + 0.5);
				//printf("GAMMAn[%d][%d]=%d\n", k, j, GAMMAn[k][j]);

			}
		}
		//printf("row=%d\n", i);
	}

#ifdef DEBUG
	/*for (i = 0; i < N; i++) {
		for (j = 0; j < Q; j++) {
			fprintf(fp, "%d\n", GAMMAn[j][i]);
		}
	}
	fclose(fp);*/
#endif
}

//===================================
// Read H and gamma
//===================================
void read_data_Q256() {

	FILE *fp_params;
	FILE *fp_hnb;
	FILE *fp_gamma;
	FILE *fp_encoded;
	FILE *fp_row_weight;
	FILE *fp_col_weight;

	fp_params = fopen("random_data/Q256/params.txt", "r");
	fp_hnb = fopen("random_data/Q256/hnb_vertical.txt", "r");
	fp_gamma = fopen("random_data/Q256/symbol_llr.txt", "r");
	fp_encoded = fopen("random_data/Q256/encoded_sym.nb.txt", "r");
	fp_row_weight = fopen("random_data/Q256/row_weight.txt", "r");
	fp_col_weight = fopen("random_data/Q256/col_weight.txt", "r");



	char line[5100*4];
	int d = 0;
	int u;
	unsigned int i;

	while (fgets(line, sizeof(line), fp_params)) {
		sscanf(line, "%d", &u);
		//printf("params %d\n", u);
		if(d==0){
			M = u;
		}
		else if(d==1){
			N = u;
		}
		else if(d==2){
			dc = u;
		}
		else if(d==3){
			dv = u;
		}
		d++;
	}

	printf("params check\n");




	int p=0;
	d=0;
	while (fgets(line, sizeof(line), fp_hnb)) {
		sscanf(line, "%d", &i);
		H_nb[d][p] = i;
		p++;
		if(p>=N){
			p=0;
			d++;
		}
	}

		/* for (int j = 0; j < N; j++) {
			printf("H_nb[0][j]=%d\n",H_nb[0][j]);
		} */


	printf("H_nb check\n");

	d=0;
	while (fgets(line, sizeof(line), fp_row_weight)) {
		sscanf(line, "%d", &i);
		//printf("%d\n", i);
		row_weight[d] = i;
		d++;
	}

	printf("row_weight check\n");

	d=0;
	while (fgets(line, sizeof(line), fp_col_weight)) {
		sscanf(line, "%d", &i);
		//printf("%d\n", i);
		col_weight[d] = i;
		d++;
	}

	printf("col_weight check\n");

	/* for (int x = 0 ; x < nnz+1 ; x++){
		printf("ptr_to_val[%d]=%d\n", x, col_ind[x]);
	} */




	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < Q; k++) {
				BETAmn[i][j][k] = 0;
				ALPHAmn[i][j][k] = 0;
			}
		}
	}

	printf("alpha and beta init\n");


	for (int i = 0; i < Q; i++) {
		for (int j = 0; j < N; j++) {
			GAMMAn[i][j] = 0;
		}
	}

	printf("gamma init\n");

	d=0;
	int z=0;
	while (fgets(line, sizeof(line), fp_gamma)) {
		sscanf(line, "%d", &i);
		GAMMAn[d][z] = i;
		//printf("GAMMAn[%d][%d] = %d\n", d, z, GAMMAn[d][z]);
		d++;
		if (d>=Q){
			d = 0;
			z++;
		}
	}

	printf("gamma check\n");

	d=0;
	while (fgets(line, sizeof(line), fp_encoded)) {
		sscanf(line, "%d", &i);
		//printf("%d\n", i);
		encoded_sym[d] = i;
		d++;
	}

	printf("codeword check\n");

	for (int i = 0; i < M; i++) {
		for (int j = 0; j < row_weight[i]; j++) {
			for (int k = 0; k < N; k++) {
				if (H_nb[i][k] != 0) {
					col_row[i][j] = k;
					//fprintf(fp, "%d ", col_row[i][j]);
					j++;
				}
			}
		}
		//fprintf(fp, "\n");
	}

	printf("col_row check\n");

	for (int j = 0; j < N; j++) {
		for (int i = 0; i < col_weight[j]; i++) {
			for (int k = 0; k < M; k++) {
				if (H_nb[k][j] != 0) {
					row_col[i][j] = k;
					//fprintf(fp, "%d ", row_col[i][j]);
					i++;
				}
			}
		}
		//fprintf(fp, "\n");
	}

	printf("row_col check\n");

	fclose(fp_params);
	fclose(fp_row_weight);
	fclose(fp_col_weight);
	fclose(fp_hnb);
	fclose(fp_gamma);
	fclose(fp_encoded);

	printf("read exit\n");



}
//===================================
// Read CSR-CSC data format
//===================================
void read_data() {



	FILE *fp_gamma;


	if (Q == 4){

		fp_gamma = fopen("C:/Users/mtech/Desktop/Oscar ferraz/NB_LDPC-master/nbldpc-cpu_ubuntu/random_data/Q4/symbol_llr.txt", "r");


	}

	else if (Q == 8){

		fp_gamma = fopen("C:/Users/mtech/Desktop/Oscar_ferraz/NB_LDPC-master/nbldpc-cpu_ubuntu/random_data/Q8/symbol_llr.txt", "r");

	}
	else if (Q == 16){

		fp_gamma = fopen("C:/Users/mtech/Desktop/Oscar_ferraz/NB_LDPC-master/nbldpc-cpu_ubuntu/random_data/Q16/symbol_llr.txt", "r");

	}
	else if (Q == 32){

			fp_gamma = fopen("C:/Users/mtech/Desktop/Oscar_ferraz/NB_LDPC-master/nbldpc-cpu_ubuntu/random_data/Q32/symbol_llr.txt", "r");

		}





	/*
	else if (Q == 16){
		fp_params = fopen("random_data/Q16/params.txt", "r");
		fp_row_weight = fopen("random_data/Q16/row_weight.txt", "r");
		fp_col_weight = fopen("random_data/Q16/col_weight.txt", "r");
		fp_val = fopen("random_data/Q16/val.txt", "r");
		fp_col_ind = fopen("random_data/Q16/col_ind.txt", "r");
		fp_row_ptr = fopen("random_data/Q16/row_ptr.txt", "r");
		fp_col_ptr = fopen("random_data/Q16/col_ptr.txt", "r");
		fp_ptr_to_val = fopen("random_data/Q16/ptr_to_val.txt", "r");
		fp_gamma = fopen("random_data/Q16/symbol_llr.txt", "r");
		fp_encoded = fopen("random_data/Q16/encoded_sym.nb.txt", "r");
		fp_row_ind = fopen("random_data/Q16/row_ind.txt", "r");
	}
	else if (Q == 32){
		fp_params = fopen("random_data/Q32/params.txt", "r");
		fp_row_weight = fopen("random_data/Q32/row_weight.txt", "r");
		fp_col_weight = fopen("random_data/Q32/col_weight.txt", "r");
		fp_val = fopen("random_data/Q32/val.txt", "r");
		fp_col_ind = fopen("random_data/Q32/col_ind.txt", "r");
		fp_row_ptr = fopen("random_data/Q32/row_ptr.txt", "r");
		fp_col_ptr = fopen("random_data/Q32/col_ptr.txt", "r");
		fp_ptr_to_val = fopen("random_data/Q32/ptr_to_val.txt", "r");
		fp_gamma = fopen("random_data/Q32/symbol_llr.txt", "r");
		fp_encoded = fopen("random_data/Q32/encoded_sym.nb.txt", "r");
		fp_row_ind = fopen("random_data/Q32/row_ind.txt", "r");
	}
	else if (Q == 64){
		fp_params = fopen("random_data/Q64/params.txt", "r");
		fp_row_weight = fopen("random_data/Q64/row_weight.txt", "r");
		fp_col_weight = fopen("random_data/Q64/col_weight.txt", "r");
		fp_val = fopen("random_data/Q64/val.txt", "r");
		fp_col_ind = fopen("random_data/Q64/col_ind.txt", "r");
		fp_row_ptr = fopen("random_data/Q64/row_ptr.txt", "r");
		fp_col_ptr = fopen("random_data/Q64/col_ptr.txt", "r");
		fp_ptr_to_val = fopen("random_data/Q64/ptr_to_val.txt", "r");
		fp_gamma = fopen("random_data/Q64/symbol_llr.txt", "r");
		fp_encoded = fopen("random_data/Q64/encoded_sym.nb.txt", "r");
		fp_row_ind = fopen("random_data/Q64/row_ind.txt", "r");
	}
	else if (Q == 128){
		fp_params = fopen("random_data/Q128/params.txt", "r");
		fp_row_weight = fopen("random_data/Q128/row_weight.txt", "r");
		fp_col_weight = fopen("random_data/Q128/col_weight.txt", "r");
		fp_val = fopen("random_data/Q128/val.txt", "r");
		fp_col_ind = fopen("random_data/Q128/col_ind.txt", "r");
		fp_row_ptr = fopen("random_data/Q128/row_ptr.txt", "r");
		fp_col_ptr = fopen("random_data/Q128/col_ptr.txt", "r");
		fp_ptr_to_val = fopen("random_data/Q128/ptr_to_val.txt", "r");
		fp_gamma = fopen("random_data/Q128/symbol_llr.txt", "r");
		fp_encoded = fopen("random_data/Q128/encoded_sym.nb.txt", "r");
		fp_row_ind = fopen("random_data/Q128/row_ind.txt", "r");
	}
	else if (Q == 256){
		fp_params = fopen("random_data/Q256/params.txt", "r");
		fp_row_weight = fopen("random_data/Q256/row_weight.txt", "r");
		fp_col_weight = fopen("random_data/Q256/col_weight.txt", "r");
		fp_val = fopen("random_data/Q256/val.txt", "r");
		fp_col_ind = fopen("random_data/Q256/col_ind.txt", "r");
		fp_row_ptr = fopen("random_data/Q256/row_ptr.txt", "r");
		fp_col_ptr = fopen("random_data/Q256/col_ptr.txt", "r");
		fp_ptr_to_val = fopen("random_data/Q256/ptr_to_val.txt", "r");
		fp_gamma = fopen("random_data/Q256/symbol_llr.txt", "r");
		fp_encoded = fopen("random_data/Q256/encoded_sym.nb.txt", "r");
		fp_row_ind = fopen("random_data/Q256/row_ind.txt", "r");
	}

*/
	char line[256];
	int d = 0;
	int u;
	unsigned int i;












	//GAMMAn[Q][N];

	for (int z = 0; z < Q; z++){
		for (int j = 0; j < N; j++){
			GAMMAn[z][j] = 0;
		}
	}

	d=0;
	int z=0;
	while (fgets(line, sizeof(line), fp_gamma)) {
		sscanf(line, "%d", &i);
		GAMMAn[d][z] = i;
		//printf("GAMMAn[%d][%d] = %d\n", d, z, GAMMAn[d][z]);
		d++;
		if (d>=Q){
			d = 0;
			z++;
		}
	}


	fclose(fp_gamma);

}


//===================================
// Main Function
//===================================

int main(int argc, char*argv[])
{int		*info_bin;
int		*decoded_bit;
int		*info_sym;
int		*info_bin_coded;
double 	*trans;
double	*recv;
int		*hard_decision_bit;
double 	*llr;

double 	rate = 0.5;
double 	snr;
double 	sigma;
int 	iteration;
//int 	seed = 69012;
int 	error = 0;
int 	error_num = 0;
int 	fer_threshold = 100;
int 	counter = 0;
int 	fe = 0;

#if (RANDOM_BITS == 1)
srand((unsigned int)time(NULL));
#else
srand((unsigned int)1010);
#endif

if(Q==256){
	read_data_Q256();

	for (snr = SNR_LOW; snr <= SNR_HIGH; snr += SNR_DELTA){

		double total_time = 0.0;


		t1 = clock();
		error = minmax(BETAmn, ALPHAmn, GAMMAn);
		t2 = clock();
		double time_diff = 1000 * double(t2 - t1) / CLOCKS_PER_SEC;
		total_time += time_diff;
		fprintf(stderr, "throughput is %f Kbps\n", counter * 310 * 5 / total_time);

		error_num += error;
		if (error != 0) fe++;
		printf("frame = %d frame error = %d, bit_error=%d\n", counter, fe, error);
	}


	printf("SER is %f at SNR = %f\n", (double)error_num / (counter*(N - M)), snr);
	printf("FER is %f at SNR = %f\n", (double)fe / (counter), snr);
	printf("SNR = %f is finished\n", snr);


}
else{

	FILE *fp = fopen("data/data_new.txt", "w+");
	//FILE *fpp = fopen("matrix/H_nb_main.txt", "w+");
	//srand (seed);

	qc_gen();	// call quicy-cyclic generate function
	trans = (double *)malloc(LOG2Q * N * sizeof(double));
	recv = (double *)malloc(LOG2Q * N * sizeof(double));
	llr = (double *)malloc(LOG2Q * N * sizeof(double));
	info_bin = (int *)malloc(LOG2Q * (N - M) * sizeof(int));
	decoded_bit = (int *)malloc(LOG2Q * (N - M) * sizeof(int));
	info_sym = (int *)malloc((N - M) * sizeof(int));
	info_bin_coded = (int *)malloc(LOG2Q * N * sizeof(int));
	hard_decision_bit = (int*)malloc(LOG2Q * N * sizeof(int));

	printf("Q is %d\n", Q);
	h2g_nb(H_nb, N);
	fprintf(fp, "%d\n", Q);
	printf("check\n");

	for (snr = SNR_LOW; snr <= SNR_HIGH; snr += SNR_DELTA)
	{
	#ifdef NONOISE
			sigma = 0.01;
	#else
			sigma = 1 / sqrt(2.0*rate*pow(10.0, (snr / 10.0)));
	#endif
			counter = 0;
			error_num = 0;
			fe = 0;
double total_time = 0.0;

			counter++;
			info_gen(info_bin);
			grouping(info_bin, info_sym);
			encoder_nb(info_sym, encoded_sym);
			ungrouping(encoded_sym, info_bin_coded);
			bpsk_modulation(info_bin_coded, trans);
			awgn(trans, recv, sigma);

			/* for(int i = 0;i < M;i++){
				for(int j = 0;j < N;j++){
					if( H_nb[i][j]!=0)
					printf("hnb[%d][%d]=%hu\n", i, j, H_nb[i][j] );
				}
			}  */

			t1 = clock();
			hard_decision(recv, sigma, hard_decision_bit);
			sym_llr(llr, BETAmn, GAMMAn, hard_decision_bit, recv, sigma);
			read_data();
			/*for (int j=0;j<N;j++)
						{
							for(int k=0;k<Q;k++)
							{
					printf("%hhu \n",GAMMAn[k][j]);
							}
						}
*/
			error = minmax(BETAmn, ALPHAmn, GAMMAn);

	/*		for(int i=0;i<U;i++)
			{
				for (int j=0;j<V;j++)
				{
					for(int k=0;k<Q;k++)
					{
			printf("%hhu \n",BETAmn[i][j][k]);
					}}}

			*/

					t2 = clock();
			double time_diff = 1000 * double(t2 - t1) / CLOCKS_PER_SEC;
			total_time += time_diff;
			fprintf(stderr, "throughput is %f Kbps\n", counter * 6 * 3 / total_time);

			error_num += error;
			if (error != 0) fe++;
			printf("frame = %d frame error = %d, bit_error=%d\n", counter, fe, error);
		}

		fprintf(fp, "SER is %f at SNR = %f\n", (double)error_num / (counter*(N - M)), snr);
		fprintf(fp, "FER is %f at SNR = %f\n", (double)fe / (counter), snr);
		printf("SER is %f at SNR = %f\n", (double)error_num / (counter*(N - M)), snr);
		printf("FER is %f at SNR = %f\n", (double)fe / (counter), snr);
		printf("SNR = %f is finished\n", snr);

return 0;

}}
