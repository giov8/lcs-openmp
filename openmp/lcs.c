#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

typedef unsigned short mtype;

/* Read sequence from a file to a char vector.
 Filename is passed as parameter */

char* read_seq(char *fname) {
	//file pointer
	FILE *fseq = NULL;
	//sequence size
	long size = 0;
	//sequence pointer
	char *seq = NULL;
	//sequence index
	int i = 0;

	//open file
	fseq = fopen(fname, "rt");
	if (fseq == NULL ) {
		printf("Error reading file %s\n", fname);
		exit(1);
	}

	//find out sequence size to allocate memory afterwards
	fseek(fseq, 0L, SEEK_END);
	size = ftell(fseq);
	rewind(fseq);

	//allocate memory (sequence)
	seq = (char *) calloc(size + 1, sizeof(char));
	if (seq == NULL ) {
		printf("Erro allocating memory for sequence %s.\n", fname);
		exit(1);
	}

	//read sequence from file
	while (!feof(fseq)) {
		seq[i] = fgetc(fseq);
		if ((seq[i] != '\n') && (seq[i] != EOF))
			i++;
	}
	//insert string terminator
	seq[i] = '\0';

	//close file
	fclose(fseq);

	//return sequence pointer
	return seq;
}

mtype ** allocateMatrix(int sizeA, int sizeB) {
	int i;
	//Allocate memory for LCS score matrix
	mtype ** scoreMatrix = (mtype **) calloc((sizeB + 1), sizeof(mtype *));
	for (i = 0; i < (sizeB + 1); i++)
		scoreMatrix[i] = (mtype *) calloc((sizeA + 1), sizeof(mtype));
	return scoreMatrix;
}

void initScoreMatrix(mtype ** scoreMatrix, int sizeA, int sizeB) {
	int i, j;
	//Fill first line of LCS score matrix with zeroes
	for (j = 0; j < (sizeA + 1); j++)
		scoreMatrix[0][j] = 0;

	//Do the same for the first collumn
	for (i = 1; i < (sizeB + 1); i++)
		scoreMatrix[i][0] = 0;
}

char* findUniqueChars(int sizeA, int sizeB, int *sizeUnique, char *seqA, char *seqB) {
	// this function allocates and returns an array of unique chars
	
	char *seqC = (char *) malloc((sizeA+sizeB)*sizeof(char));

	seqC[0] = seqA[0];
	int sizeC = 1;
	int i, c;

	for (i = 1; i < sizeA; i++) {
		for (c = 0; c < sizeC; c++) {
			if (seqC[c] == seqA[i]) {
				break;
			}
		}
		if (c == sizeC) {
			seqC[sizeC] = seqA[i];
			sizeC++;
		}
	}

	for (i = 0; i < sizeB; i++) {
		for (c = 0; c < sizeC; c++) {
			if (seqC[c] == seqB[i]) {
				break;
			}
		}
		if (c == sizeC) {
			seqC[sizeC] = seqB[i];
			sizeC++;
		}
	}

	seqC = (char *) realloc(seqC, sizeC*sizeof(char));
	*sizeUnique = sizeC;

	return (char *) realloc(seqC, (sizeC+1)*sizeof(char));
}

//old
// void fillPMatrix(mtype **P, char *c, int len_c, char *b, int len_b)
// {
//     //#pragma omp parallel for
//     for(int i = 0; i < len_c; i++)
//     {
//         for(int j = 2; j < len_b+1; j++)
//         {
//             if(b[j-2]==c[i]) //j-2 as b we assume here that b has a empty character in the beginning
//             {
//                 P[i][j] = j-1;
//             }
//             else
//             {
//                 P[i][j] = P[i][j-1];
//             }
//         }
//     }
// }

void fillPMatrix(mtype **P, char *c, int len_c, char *b, int len_b)
{
	printf("\nFIL ");
    //#pragma omp parallel for
    for(int i = 0; i < len_c; i++)
    {
		printf("\ni na fill: %d \n", i);
        for(int j = 1; j < len_b+1; j++)
        {
			printf("j na fill: %d \n", j);
			printf("b[j-1]: %c  c[i] %c \n", b[j-1], c[i]);
			if(b[j-1] == c[i])
            {
                P[i][j] = j - 1;
            }

            else
            {
                P[i][j] = P[i][j-1];
            }
        }
    }
}

int LCS(mtype ** scoreMatrix, int sizeA, int sizeB, char * seqA, char *seqB) {
	int i, j;
	for (i = 1; i < sizeB + 1; i++) {
		for (j = 1; j < sizeA + 1; j++) {
			if (seqA[j - 1] == seqB[i - 1]) {
				/* if elements in both sequences match,
				 the corresponding score will be the score from
				 previous elements + 1*/
				scoreMatrix[i][j] = scoreMatrix[i - 1][j - 1] + 1;
			} else {
				/* else, pick the maximum value (score) from left and upper elements*/
				scoreMatrix[i][j] =
						max(scoreMatrix[i-1][j], scoreMatrix[i][j-1]);
			}
		}
	}
	return scoreMatrix[sizeB][sizeA];
}

int getCharIndex(char *str, int len, char x)
{
    for(int i = 0; i < len; i++)
    {
        if(str[i] == x)
        {
            return i;
        }
    }
    return -1; // caractere não encontrado
}

int LcsParallel(mtype ** scoreMatrix, int sizeA, int sizeB, char * seqA, char *seqB, mtype ** pMatrix, int sizeUniqueChars, char *seqUniqueChars) {
	
	int i, j;
	for (i = 1; i < sizeB + 1; i++) {
		int c_i = getCharIndex(seqUniqueChars, sizeUniqueChars, seqA[i-1]); // the characther index in P Matrix
		printf("i: %d  CI: %d ", i, c_i);
		for (j = 1; j < sizeA + 1; j++) {
			if (seqA[j - 1] == seqB[i - 1])
			{
				scoreMatrix[i][j] = scoreMatrix[i - 1][j - 1] + 1;
			}
			else if (pMatrix[c_i][j] == 0)
			{
				scoreMatrix[i][j] =
						max(scoreMatrix[i-1][j], 0);
			}
			else
			{
				scoreMatrix[i][j] = max(scoreMatrix[i-1][j], scoreMatrix[i-1][pMatrix[c_i][j]-1] + 1);
			}
		}
	}
	return scoreMatrix[sizeB][sizeA];
}

// short lcs_yang_v1(short **DP, short **P, char *A, char *B, char *C, int m, int n, int u)
// {
//     for(int i=1;i<m+1;i++)
//     {
//         int c_i = get_index_of_character(C,A[i-1],u);
//         //#pragma omp parallel for schedule(static)
// 	for(int j=0;j<n+1;j++)
//         {
//             if(A[i-1]==B[j-1])
//             {
//                 DP[i][j] = DP[i-1][j-1] + 1;
//             }
//             else if(P[c_i][j]==0)
//             {
//                 DP[i][j] = max(DP[i-1][j], 0);
//             }
//             else
//             {
//                 DP[i][j] = max(DP[i-1][j], DP[i-1][P[c_i][j]-1] + 1);
//             }
//         }
//     }
//     return DP[m][n];
// }

void printMatrix(char * seqA, char * seqB, mtype ** scoreMatrix, int sizeA,
		int sizeB) {
	int i, j;

	//print header
	printf("Score Matrix:\n");
	printf("========================================\n");

	//print LCS score matrix allong with sequences

	printf("    ");
	printf("%5c   ", ' ');

	for (j = 0; j < sizeA; j++)
		printf("%5c   ", seqA[j]);
	printf("\n");
	for (i = 0; i < sizeB + 1; i++) {
		if (i == 0)
			printf("    ");
		else
			printf("%c   ", seqB[i - 1]);
		for (j = 0; j < sizeA + 1; j++) {
			printf("%5d(%d,%d)   ", scoreMatrix[i][j], i,j);
		}
		printf("\n");
	}
	printf("========================================\n");
}

void freeScoreMatrix(mtype **scoreMatrix, int sizeB) {
	int i;
	for (i = 0; i < (sizeB + 1); i++)
		free(scoreMatrix[i]);
	free(scoreMatrix);
}

int main(int argc, char ** argv) {
	// sequence pointers for both sequences
	char *seqA, *seqB, *uniqueChars;

	// sizes of both sequences
	int sizeA, sizeB, sizeUniqChars;

	//read both sequences
	seqA = read_seq("fileA.in");
	seqB = read_seq("fileB.in");

	//find out sizes
	sizeA = strlen(seqA);
	sizeB = strlen(seqB);

	// allocate LCS score matrix
	mtype ** scoreMatrix = allocateMatrix(sizeA, sizeB);

	// allocate LCS score matrix
	mtype ** scoreMatrix2 = allocateMatrix(sizeA, sizeB);
	initScoreMatrix(scoreMatrix2, sizeA, sizeB);

	// initialize LCS score matrix
	initScoreMatrix(scoreMatrix, sizeA, sizeB);

	// allocate unique chars' array
	uniqueChars = findUniqueChars(sizeA, sizeB, &sizeUniqChars, seqA, seqB);


	printf("size uniq: %d\n", sizeUniqChars);
	printf("%s ", uniqueChars);
	printf("\n");
	
	// allocate P matrix
	mtype ** pMatrix = allocateMatrix(sizeB-1, sizeUniqChars-1);

	// fill P matrix
	fillPMatrix(pMatrix, uniqueChars, sizeUniqChars, seqB, sizeB);

	for (int i = 0; i < sizeUniqChars; i++) {
		for (int j = 0; j < sizeB; j++) {
			printf("%d(%d,%d)   ", pMatrix[i][j],i,j);
		}
		printf("\n");
	}

	//printMatrix(uniqueChars, seqB, pMatrix, sizeUniqChars, sizeB);

	//fill up the rest of the matrix and return final score (element locate at the last line and collumn)
	mtype score = LCS(scoreMatrix, sizeA, sizeB, seqA, seqB);

	mtype score2 = LcsParallel(scoreMatrix2, sizeA, sizeB, seqA, seqB, pMatrix, sizeUniqChars, uniqueChars);

	/* if you wish to see the entire score matrix,
	 for debug purposes, define DEBUGMATRIX. */
#ifdef DEBUGMATRIX
	printf("SEQUENCIAL\n");
	printMatrix(seqA, seqB, scoreMatrix, sizeA, sizeB);
	printf("PARALELO\n");
	printMatrix(seqA, seqB, scoreMatrix2, sizeA, sizeB);
#endif

	//print score
	printf("\nScore: %d\n", score);
	//print score
	printf("\nScore2: %d\n", score2);

	//free score matrix
	freeScoreMatrix(scoreMatrix, sizeB);

	return EXIT_SUCCESS;
}
