#include <iostream>
#include <string>
#include<fstream>
#define U 150
#define V 300
using namespace std;
int main()
{
unsigned int i,j,k;
FILE *fp1 =fopen("/home/srinivasan/NB_LDPC/required_files/Q16/hnb.txt", "r");
FILE *fp2 =fopen ("/home/srinivasan/GF16/row_weight.txt","w");
FILE *fp3 =fopen("/home/srinivasan/GF16/col_weight.txt","w");
FILE *fp4 =fopen("/home/srinivasan/GF16/col_row.txt","w");
FILE *fp5 =fopen("/home/srinivasan/GF16/row_col.txt","w");
FILE*fp6 =fopen("/home/srinivasan/GF16/Hmatrix.txt","w");
int H[U][V];
int row_weight[U];
int col_weight[V];
int col_row[150][2];
int row_col[1][300];
if(fp1==NULL) exit(1);
	for (i=0;i<150;i++)
	{
		for (j=0;j<300;j++)
		{
			fscanf(fp1,"%d",&H[i][j]);
		}
	}
	
for (i=0;i<150;i++)
{
	for(j=0;j<300;j++)
	{
		if (H[i][j]!=0)
				{
		fprintf(fp6,"H[%d][%d]= %d \n",i,j,H[i][j]);

	}
	}}
fclose(fp6);
	for (i = 0; i < U; i++) {
		for (j = 0; j < V; j++) {
			if (H[i][j] != 0) {
				row_weight[i]++;
			}
		}
		fprintf(fp2, "row_weight[%d] = %d\n", i, row_weight[i]);
	}
fclose(fp2);

	for (j = 0; j < V; j++) {
		for (i = 0; i < U; i++) {
			if (H[i][j] != 0) {
				col_weight[j]++;
			}
		}
		fprintf(fp3, "col_weight[%d] = %d\n", j, col_weight[j]);
	}
fclose(fp3);
	for (i = 0; i < U; i++) {
		for (j = 0; j < row_weight[i]; j++) {
			for (k = 0; k < V; k++) {
				if (H[i][k] != 0) {
					col_row[i][j] = k;
					fprintf(fp4, "col_row[%d][%d]= %d\n",i,j,col_row[i][j]);
					j++;
				}
			}
		}
		
	}
fclose(fp4);

	for (j = 0; j < V; j++) {
		for (i = 0; i < col_weight[j]; i++) {
			for (k = 0; k < U; k++) {
				if (H[k][j] != 0) {
					row_col[i][j] = k;
					fprintf(fp5, "row_col[%d][%d] =%d \n", i,j,row_col[i][j]);
					i++;
				}
			}
		}
		
	}

fclose(fp5);
return 0;
}
