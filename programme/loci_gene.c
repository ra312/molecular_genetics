#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

// defining the initial structure of the probability space
double w1=1.0; /* the probability measure for the genes set*/
double w2=1.0; /* the probability measure for the loci*/

int gene_values_num=0;
int SENSITIVITY = 3;
typedef struct
{
	int chr_num; 			//chromosome number in the genome
	char chr[15];			//chromosome number in string to include the case of X and Y chromosome
	long long start;		//start position in chromosome
	long long end;			// end position in chromosome
	char name[30]; //gene name
	double value;
	double min;
	double max;
	double mean;
	double mode;
	int omega;
	long long matched; // the number of values
	double* values;
	double dispersion;
	int sign;
	int* signs;
} tLocus;

typedef struct
{
	int chr_num; 			//chromosome number in the genome
	char chr[15];			//chromosome number in string to include the case of X and Y chromosome
	long long start;		//start position in chromosome
	long long end;			// end position in chromosome
	char name[30]; //gene name
	double value;
	double min;
	double max;
	double mean;
	double mode;
	int omega;
	long long matched; // the number of values
	double* values;
	double dispersion;
} tGene;

tLocus* AllocateMoreMemToLoci(tLocus* ReadLoci, int NUMBER_OF_LOCI, int N)
{
	tLocus* Loci = (tLocus*) malloc(sizeof(tLocus)*(N	));
	int i;

	for (i=0; i<NUMBER_OF_LOCI; i++)
	{
		Loci[i].chr_num=ReadLoci[i].chr_num;
		strcpy(Loci[i].chr, ReadLoci[i].chr);
		Loci[i].start=ReadLoci[i].start;		
		Loci[i].end=ReadLoci[i].end;		
		strcpy(Loci[i].name,ReadLoci[i].name); 
		Loci[i].value=ReadLoci[i].value;
		Loci[i].min=ReadLoci[i].min;
		Loci[i].max=ReadLoci[i].max;
	}
	free(ReadLoci);
	return Loci;
}
//char gene_name[300]={0};
int* AllocateMoreMemToIntArray(int* array, int old_size, int new_size)
{
	//printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	//getchar();
	int* new_array = (int*)malloc(sizeof(int)*new_size);
	int i;
	if (new_size<old_size)
	{
		printf("Please enter bigger array size\n");
		return NULL;
	}
	for (i=0; i<old_size; i++)
	{
		new_array[i]=array[i];
	}
	array = NULL;
	//printf("Inside  int* AllocateMoreMemToIntArray(int* array, int old_size, int new_size)\n");
	// for (i=0; i<old_size; i++)
// 	{
// 		printf("new_array[%i]=%i\n",i, new_array[i]);
// 	}
	return new_array;
}

tLocus* ReadLocusValues(char* genes_file, int* NUMBER_OF_LOCI)
{
	int NUM = 50000;
	tLocus* Loci = (tLocus*) malloc(sizeof(tLocus)*NUM);

	char* chr;
	char* chr_num;

	chr = (char*) malloc(sizeof(char)*5);

	long long start;
	long long end;
	double value;

	int chrom_num;
	FILE *file;
	
	if ( !(file=fopen(genes_file,"r")))
	{
		printf("Cannot open file %s\n",genes_file);
		exit(-1);
	}
	
	int k;
	k=0;

	while ( (fscanf(file,"%s %lld %lld %lf",chr,&start,&end,&value))!=EOF)
	{
		chrom_num=atoi(chr+3);
		if (strstr(chr,"Y")) 
		{
			chrom_num=-2;
		}
		else if (strstr(chr,"X"))
		{
			chrom_num=-1;
		}
		else 
		{
		
			chrom_num=atoi(chr+3);
		
		}
		strcpy(Loci[k].chr,chr);
		Loci[k].chr_num=chrom_num;
		
		Loci[k].start=start;
		Loci[k].end=end;
		
		Loci[k].value=value;
		Loci[k].min=0.0;
		Loci[k].min=0.0;
		k++;
		if (k==NUM)
		{
			printf("Allocating more memory ... \n");
			Loci=AllocateMoreMemToLoci(Loci, NUM,2*NUM);
			NUM=2*NUM;
		}
		
	}

	*NUMBER_OF_LOCI=k;
	fclose(file);

	return Loci;
}


void PrintLocuss(tLocus* Loci, int N, char* parameter)
{
	int k;
	int j;
	
	if (strcmp(parameter,"Genes")==0)
	{
		for (k=0; k<N; k++)
		{		
				
			printf("[%3i]:gene=%10s, chr=%5s, start=%12lld, end=%12lld\n",k,Loci[k].name, Loci[k].chr, Loci[k].start, Loci[k].end);
				
			
			for (j=0; j<Loci[k].matched; j++)
			{
				printf("\t\t  %1.8lf\n",Loci[k].values[j]);
			} 		
		}// for k

	}// if
	else if (strcmp(parameter,"Loci")==0)
	{
		
		for (k=0; k<N; k++)
		{		
printf("[%3i]:gene=%10s, chr=%5s, start=%12lld, end=%12lld, value = %1.8lf\n",k,Loci[k].name, Loci[k].chr, Loci[k].start, Loci[k].end,Loci[k].value);
		}// for k

	}// else if 
	else
	{
		printf("Wrong parameter\n");
	}
	
//printf("[%3i]:gene=%10s, chr=%5s, start=%12lld, end=%12lld, value=%1.8lf\n",k,Loci[k].name, Loci[k].chr, Loci[k].start, Loci[k].end, Loci[k].mode);
		
	return;	
}
tLocus* ReadGenes(char* selected_genes_file, int *NUMBER_OF_SELECTED_GENES)
{
	FILE *file;

	if ( !(file=fopen(selected_genes_file,"r")))
	{
		printf("Cannot open file %s\n",selected_genes_file);
		exit(-1);
	}

	int NUM	=	200;

	tLocus* SelectedGenes = (tLocus*)malloc(sizeof(tLocus)*NUM);
	
	int k;
	k=0;
	int j;
	char chr[20];
	char name[20];

	double start;
	double end;
	
	while ( (fscanf(file,"%s %lf %lf %s",chr,&start,&end,name))!=EOF)
	{
		strcpy(SelectedGenes[k].chr,chr);
		
		strcpy(SelectedGenes[k].name,name);

		SelectedGenes[k].start=(long long) start;

		SelectedGenes[k].end=(long long ) end;

		SelectedGenes[k].min=0.0;
		SelectedGenes[k].min=0.0;
		k++;
		if (k==NUM)
		{
			printf("Allocating more memory ... \n");
			SelectedGenes=AllocateMoreMemToLoci(SelectedGenes, NUM,2*NUM);
			NUM=2*NUM;
		}
		//printf("k=%i\n",k);
	}
	
	*NUMBER_OF_SELECTED_GENES=k;
	return SelectedGenes;
}
// FIND MIN VALUE AMONG ALL GENE VALUES WITHIN EACH CHROMOSOME
void MinMax(tLocus* ReadLoci, int NUMBER_OF_LOCI, tLocus* SelectedLoci, int NUMBER_OF_SELECTED_LOCI)
{
	//MinMax computes MIN and MAX per each gene in SelectedLoci
	int k,j;
	double swap;
	double min;
	double max;
	for (k=0; k<NUMBER_OF_SELECTED_LOCI; k++)
	{
		min=+200.0;
		max=-200.0;
		for (j=0; j<NUMBER_OF_LOCI; j++)
		{
			if (!strcmp(SelectedLoci[k].name,ReadLoci[j].name))
			{
				
				if (ReadLoci[j].value<min) 
				{
//					printf("Min\n");
					min=ReadLoci[j].value;
				}
				if (ReadLoci[j].value>max)
				{
							max=ReadLoci[j].value;
				} 
			}
		}//for j
		for (j=0; j<NUMBER_OF_LOCI; j++)
			if (!strcmp(SelectedLoci[k].name,ReadLoci[j].name))
			{
				ReadLoci[j].min=min;
				ReadLoci[j].max=max;
			}
	
	}//for k
			
	

	return;
}
double MeanValue(tLocus* Loci, int NUMBER_OF_LOCI, double* weights, char* parameter)
{
	double sum=0.0;
	int k, counter=0;
	if (strcmp(parameter,"Loci")==0)
	{
		for (k=0; k<NUMBER_OF_LOCI; k++)
		{
			sum += Loci[k].value;
			counter++;
		};
		sum = sum / ((double)NUMBER_OF_LOCI);

	}
	else if (strcmp(parameter,"Genes")==0)
	{
		for (k=0; k<NUMBER_OF_LOCI; k++)
		{
			sum += Loci[k].value*Loci[k].omega;
			counter++;
		}
		sum  = sum / ((double) gene_values_num);
		
	}
	else
	{
	}
	return sum;
}

double Dispersion(tLocus* Loci, int NUMBER_OF_LOCI, double* weights, char* parameter)
{
	double quadr_mean, arithm_mean;
	int k;
	int added;
	double volume;
	
	added=0;
	
	if (strcmp(parameter,"Genes")==0)
	{
			
		quadr_mean=arithm_mean=0.0;
		volume = 0.0;
		for (k=0; k<NUMBER_OF_LOCI; k++)
		{
			added++;
			arithm_mean	+= Loci[k].sign*Loci[k].omega;
			volume += weights[k]*Loci[k].omega;
		}
		arithm_mean = arithm_mean / volume; // conditional expectation.
// 		printf("Genes: arithm_mean = %f\n",arithm_mean);
		for (k=0; k<NUMBER_OF_LOCI; k++)
		{
			// printf("Loci[%i].value=%f\n", k, Loci[k].value);
// 			printf("Loci[%i].omega=%i\n",k,Loci[k].omega);
// 			printf("weights[%i]=%f\n",k, weights[k]);
			quadr_mean 	+= pow(0.0+Loci[k].sign-arithm_mean, 2.0)*Loci[k].omega*weights[k];
// 			printf("Genes:	[%i]:quadr_mean = %f\n",k,quadr_mean);

		}
		quadr_mean = quadr_mean / ((float)volume);
// 		printf("volume - %f\n",volume);
//    		printf("Genes:	quadr_mean = %f\n",quadr_mean);
	
	}
	else if (strcmp(parameter,"Loci")==0)
	{
		quadr_mean=arithm_mean=0.0;
		volume = 0.0;
		for (k=0; k<NUMBER_OF_LOCI; k++)
		{	
			arithm_mean	+= Loci[k].sign*weights[k]*(1-0*Loci[k].omega);
			volume		+= weights[k]*(1-0*Loci[k].omega);
		}
		arithm_mean = arithm_mean / volume;
		for (k=0; k<NUMBER_OF_LOCI; k++)
		{
			quadr_mean +=pow(0.0+Loci[k].sign-arithm_mean, 2.0)*weights[k]*(1-0*Loci[k].omega);
		}
		quadr_mean = quadr_mean / volume;
   			//printf("Loci:	quadr_mean = %f\n",quadr_mean);

	}
	else
	{
	// 	printf("Else executed\n");
	}
	
	return quadr_mean;

}
int LociOverLap(tLocus ReadGene, tLocus SelectedGene);
void MatchValues(tLocus* ReadLoci, int NUMBER_OF_LOCI,tLocus* SelectedLoci, int NUMBER_OF_SELECTED_LOCI)
{
	double* values_buffer = (double*) malloc (sizeof(double)*NUMBER_OF_LOCI);
	int* signs_buffer = (int*) malloc(sizeof(int)*NUMBER_OF_LOCI);
	int k,i;
	for (k=0; k<NUMBER_OF_LOCI; k++) values_buffer[k]=0;
	int num_matched=0, matched=0;
	for (k=0; k<NUMBER_OF_SELECTED_LOCI;k++)
	{
		ReadLoci[k].omega=0;
	}
	double* weights = (double*) malloc(sizeof(double)*NUMBER_OF_LOCI);
	for (k=0; k<NUMBER_OF_LOCI; k++) weights[k]=0;


	for (k=0;k<NUMBER_OF_SELECTED_LOCI; k++)
	{
		matched=0;
		for (i=0; i<NUMBER_OF_LOCI; i++)
		{

			if (!strcmp(SelectedLoci[k].chr,ReadLoci[i].chr))
			{
			
				if (LociOverLap(ReadLoci[i], SelectedLoci[k]))
				{
						values_buffer[matched] = ReadLoci[i].value;
						signs_buffer[matched]=ReadLoci[i].sign;
						matched++;
						//printf("LociOverlap: i = %i\n",i);

						ReadLoci[i].omega=1;			
						strcpy(ReadLoci[i].name, SelectedLoci[k].name);
						//printf("Gene name =%s\n",ReadLoci[i].name);
						
				}//if (LociOverLap(ReadLoci[i], SelectedLoci[k]))

			}//if (!strcmp(SelectedLoci[k].chr,ReadLoci[i].chr))
		}//for i
		//printf("i = %i, matched = %i\n",i,matched);
		num_matched=num_matched+matched;
		//printf("[%i]: mathched = %i\n",k,matched);
		SelectedLoci[k].values = (double*) malloc(sizeof(double)*matched);
		SelectedLoci[k].signs = (int*) malloc(sizeof(double)*matched);
//		printf("(SelectedLoci[k].values == NULL) = %i\n",SelectedLoci[k].values==NULL);
		SelectedLoci[k].matched=matched;
		for (i=0; i<matched; i++) 		
		{
			SelectedLoci[k].values[i]=values_buffer[i];
			SelectedLoci[k].signs[i]=signs_buffer[i];
		}
		matched=0;
		//free(values_buffer);
	}//for k

	printf("I have compared %i genes\n",NUMBER_OF_SELECTED_LOCI);
	printf("There are %i values identified in %i read  genes\n",num_matched, NUMBER_OF_SELECTED_LOCI);
	gene_values_num = num_matched;
	printf("There are %i values throughout the bedfile\n",NUMBER_OF_LOCI);
	return;
}
int LociOverLap(tLocus ReadGene, tLocus SelectedGene)
{
	int NoOverLap1,NoOverLap2;
	int OverLap1, OverLap2,OverLap3,OverLap4;
	int OverLap;
	int NoOverLap;
	long long a1, b1, a2, b2;
	
	a1=(long long) ReadGene.start;
	b1=(long long)ReadGene.end;
	a2=(long long)SelectedGene.start;
	b2=(long long)SelectedGene.end;
	
	NoOverLap1 = (a1<=b1) && (b1 <= a2) && (a2<= b2);
	if (NoOverLap1)
	{
	//	printf("NoOverLap1\n");
	}
	

	NoOverLap2 = (a2<=b2) && (b2 <= a1) && (a1<= b1);
	if (NoOverLap2)
	{
	//	printf("NoOverLap2\n");
	}

	NoOverLap = NoOverLap1 && NoOverLap2;

	OverLap1 = (a1<=a2) && (a2 <= b1) && (b1 <= b2);
	if (OverLap1)
	{
	//	printf("OverLap1\n");
	}

	OverLap2 = (a2<=a1) && (a1 <= b2) && (b2 <= b1);	
	if (OverLap2)
	{
	//	printf("OverLap2\n");
	}

	OverLap3 = (a1<=a2) && (a2 <= b2) && (b2 <= b1);
	if (OverLap3)
	{
	//	printf("OverLap3\n");
	}

	OverLap4 = (a2<=a1) && (a1 <= b1) && (b1 <= b2);
	if (OverLap4)
	{
	//	printf("OverLap4\n");
	}

	OverLap =0;
	OverLap = OverLap1 || OverLap2 || OverLap3 || OverLap4;
	if (OverLap)
	{
	//	printf("OverLap\n");
	}
	if (NoOverLap)
	{
	//	printf("NoOverLap\n");
	}

	
	return OverLap;
}
int cmpfunc (const void * a, const void * b);
int cmpfunc_2 (const void* a, const void* b);
double FloatMode(double* array, int N)
{
// printf("466: entering FloatMode function\n");
int i, j,z, maxCount;
double tmp, modeValue;

int* tally=(int*)malloc(sizeof(int)*N);
// printf("\n");
// for (i=0; i<N; i++)  printf("array[%i]=%f\n",i,array[i]);
//printf("\n");
qsort(array, N, sizeof(double), cmpfunc_2);


    for (i = 0; i <N; i++) 
    {
        for(z=i+1;z<N;z++)
        {
            
            if(array[i]==array[z])
            {
                tally[i]++;
            }
        }
    }
    

    
    maxCount = 0;
    modeValue = 0;
    for (i = 0; i <N; i++) 
    {
        if (tally[i] > maxCount) 
        {
            maxCount = tally[i];
            modeValue = array[i];
        }
    }
    return modeValue;
}    
double Mode(tLocus* Loci, int NUMBER_OF_LOCI, char* parameter)    
{



int i, j,z, maxCount;
double tmp, modeValue;

int* tally=(int*)malloc(sizeof(int)*NUMBER_OF_LOCI);

qsort(Loci, NUMBER_OF_LOCI, sizeof(tLocus), cmpfunc);

    for (i = 0; i <NUMBER_OF_LOCI; i++) 
    {
        for(z=i+1;z<NUMBER_OF_LOCI;z++)
        {
            
            if(Loci[i].value==Loci[z].value)
            {
                tally[i]++;
            }
        }
    }
    
    
    maxCount = 0;
    modeValue = 0;
    for (i = 0; i <NUMBER_OF_LOCI; i++) 
    {
        if (tally[i] > maxCount) 
        {
            maxCount = tally[i];
            modeValue = Loci[i].value;
        }
    }
//    printf("\nMode value is : %f", modeValue);
    return modeValue;
}
int cmpfunc (const void * a, const void * b)
{
   return ( *(int*)a - *(int*)b );
}
int cmpfunc_2 (const void* a, const void* b)
{
   return (* (double*) a - *(double*)b);
}

void PrintGenes(tLocus* Loci, int NUMBER_OF_LOCI)
{
	int i;
	for (i=0; i<NUMBER_OF_LOCI;i++)
	{
		printf("[%i]:gene=%5s, min=%1.5f, max=%1.5f\n",i,Loci[i].name, Loci[i].min, Loci[i].max);
	}
}

double* init_weights(int NUMBER_OF_LOCI, tLocus* ReadLoci, double w1, double w2)
{
	int k;
	double* weights = (double*) malloc(sizeof(double)*NUMBER_OF_LOCI);
	
	for (k=0; k<NUMBER_OF_LOCI; k++) weights[k]=0.0;
	
	for (k=0; k<NUMBER_OF_LOCI; k++)
	if (ReadLoci[k].omega==1) 
	{
		weights[k]=w1;
	}
	else
	{
		weights[k]=w2;
	}
return weights;
}
double WindowDispersionArray(double* signs, int start, int end)
{
	int k;
	
	double quadr_mean;
	double arithm_mean;
	
	quadr_mean=arithm_mean=0.0;
	double length = (double) (end-start);
	
//	printf("WindowDispersionArray: length = %f\n",length);
	for (k=start; k<end; k++)
	{
			arithm_mean	+= signs[k];
	}
	arithm_mean = arithm_mean / length;
	for (k=start; k<end; k++)
	{
		quadr_mean 	+= (signs[k]-arithm_mean)*(signs[k]-arithm_mean);
	}
	quadr_mean = quadr_mean / length;
	return quadr_mean;
}
double WindowDispersion(tLocus* Loci, int start, int end)
{
	int k;
	
	double quadr_mean;
	double arithm_mean;
	
	quadr_mean=arithm_mean=0.0;
	double length = (double) (end-start);
	for (k=start; k<end; k++)
	{
			arithm_mean	+= Loci[k].value;
	}
	arithm_mean = arithm_mean / length;
	for (k=start; k<end; k++)
	{
		quadr_mean 	+= (Loci[k].value-arithm_mean)*(Loci[k].value-arithm_mean);
	}
	quadr_mean = quadr_mean / length;
	return quadr_mean;
}
double SignsDispersion(int* values, int N)
{
	int k;
	double quadr_mean;
	double arithm_mean;
	
	quadr_mean=arithm_mean=0.0;
	for (k=0; k<N; k++)
	{	
		// 	printf("values[%i]=%i\n",k,values[k]);
// 			printf("%i:  SignsDispersion: arithm_mean = %f\n",k,arithm_mean);
			arithm_mean	+= (1.0)*values[k];
	}
// 	printf("SignsDispersion: arithm_mean = %f\n",arithm_mean);
	arithm_mean = arithm_mean / (1.0*N);
	for (k=0; k<N; k++)
	{
		quadr_mean 	+= pow((1.0)*values[k]-arithm_mean,2);
	}
	quadr_mean = quadr_mean / (1.0*N);
	return quadr_mean;
}
double FloatDispersion(double* values, int N)
{
	int k;
	double quadr_mean;
	double arithm_mean;
	
	quadr_mean=arithm_mean=0.0;
	for (k=0; k<N; k++)
	{
			arithm_mean	+= values[k];
	}
	arithm_mean = arithm_mean / N;
	for (k=0; k<N; k++)
	{
		quadr_mean 	+= pow(values[k]-arithm_mean,2);
	}
	quadr_mean = quadr_mean / N;
	return quadr_mean;

}
typedef struct
{
long long start;
long long end;
double dispersion;
} tRegion;

tRegion FindMinWindowDispersionArray(double* signs, long long start, long long end, int sensitivity)
{
	long long length, middle;
	double left, right;
	tRegion value;
	// printf("FindMinWindowDispersionArray: start = %lld,\t end = %lld\n",start,end);
	length = end-start+1;
// 	printf("Entered FindMinWindowDispersion\n");
// 	printf("length = %lld\n",length);	
//  	printf("FindMinWindowDispersionArray: length = %lld\n",length);
	if (length >  sensitivity)
	{
		middle = (start+end)/2;

		left  = WindowDispersionArray(signs, start, middle);
		right = WindowDispersionArray(signs, middle+1, end);
//		printf("left = %f,\t right = %f\n",left,right);
		if (left <= right)
		{
//			printf("left<right\n");
		//	printf("left<right\n");
//			printf("start = %lld, middle  = %lld\n",start, middle);
			value = FindMinWindowDispersionArray(signs, start, middle, sensitivity);
			
		}
		else if (left > right)
		{
	//		printf("left>right\n");
			value = FindMinWindowDispersionArray(signs, middle+1, end,sensitivity);
		
	//		printf("start = %lld,\t end = %lld\n",start, end);
		}
		// else
// 		{
// 			printf(" I do not know what to do\n");
// 		}
		return value;
	}// if (length >  sensitivity)
	value.start = start;
	value.end = end;
	value.dispersion=WindowDispersionArray(signs, start, end);
//	printf("value.dispersion = %f\n",value.dispersion);
	//printf("value.start = %lld,\t value.end = %lld\n",value.start,value.end);
	//printf("Exited FindMinWindowDispersion\n");
	return value;
}
tRegion FindMinWindowDispersion(tLocus* Loci, long long start, long long end, int sensitvity)
{
	long long length, middle;
	double left, right;
	tRegion value;

	length = end-start+1;
// 	printf("Entered FindMinWindowDispersion\n");
// 	printf("length = %lld\n",length);	
// 	printf("start = %lld,\t end = %lld\n",start, end);
	if (length >  sensitvity)
	{
		middle = (start+end)/2;

		left  = WindowDispersion(Loci, start, middle);
		right = WindowDispersion(Loci, middle+1, end);
//		printf("left = %f,\t right = %f\n",left,right);
		if (left <= right)
		{

//			printf("left<right\n");
			value = FindMinWindowDispersion(Loci, start, middle,sensitvity);
			value.dispersion = left;
	//		printf("value.start = %lld,\t value.end = %lld\n",value.start, value.end);
		}
		else if (left > right)
		{
	//		printf("left>right\n");
			value = FindMinWindowDispersion(Loci, middle+1, end,sensitvity);
			value.dispersion = right;
	//		printf("start = %lld,\t end = %lld\n",start, end);
		}
// 		else
// 		{
// 			printf(" I do not know what to do\n");
// 		}
		return value;
	}
	value.start = start;
	value.end = end;
	value.dispersion=WindowDispersion(Loci, start, middle);
	//printf("value.start = %lld,\t value.end = %lld\n",value.start,value.end);
	//printf("Exited FindMinWindowDispersion\n");
	return value;

}
void FindGenes(tLocus* Loci,  int start, int end)
{
	int k;

	int N = end/2; // nearest odd number, we will miss only one gene
	// range 0 ... N
	double left;
	double right;
	int left_range;
	int right_range;
	left = WindowDispersion(Loci, start,N);
	left_range = N-start+1;
	right_range = end-N;
	right = WindowDispersion(Loci,N+1, end);
	// we minimize dispersion
	
	if ( (left < right) && ( left_range>3 ))
	{
		FindGenes(Loci, start, N);
	}
	else if ( (left > right) && ((end-N)>3))
	{
		FindGenes(Loci, N+1, end);
	}
	else
	{
		
	}
	return;	
}
void oldFindGenes(tLocus* Loci)
{
	int k;
	int N;
	int start, end;
	double min = WindowDispersion(Loci, 5*k,5*(k+1));
	start = 0;
	end = 5;
	for (k=0; k<N; k++)
	{	
		
		if (min < WindowDispersion(Loci, 5*k,5*(k+1)))
		{
			start = 5*k;
			end = 5*(k+1);
		}
	}
	printf("start = %i\n", start);
	printf("end   = %i\n", end);
	return;
	
}
void Compute_SignsDispersion_per_Genes(tLocus* Genes, int N)
{
	int i;
	int j;
//	printf("In Compute_Dispersion_per_Genes\n");
	for (i=0; i<N; i++)
	{
		//printf("Genes[%i].signs==NULL=%i\n",i,Genes[i].signs==NULL);
// 		for (j=0; j<N; j++)
// 		{
// 			printf("Genes[%i].signs[%i]=%i\n",i, j,Genes[i].signs[j]);
// 		}
		Genes[i].dispersion = SignsDispersion(Genes[i].signs, Genes[i].matched);
// 		printf("Genes[%i].dispersion = %f\n",i, Genes[i].dispersion);
	}
	return;
}
void Compute_Mean_per_Genes(tLocus* Genes, int N)
{
	int k,i;
	double mean=0.0;
	for (k=0; k<N; k++)
	{
		mean = 0.0;
		for (i=0; i<N; i++)
		{
			mean  = mean+ Genes[k].values[i];
			printf("k= %i, i= %i: mean = %f\n",k,i,mean);
		}
		Genes[k].mean = mean/(1.0* Genes[k].matched);
		printf("Genes[%i].matched = %lld\n",k, Genes[k].matched);

		printf("Genes[%i].mean = %f\n",k, Genes[k].mean);
	}
	return;
}
void Compute_Mode_per_Genes(tLocus* Genes, int N)
{
	
	int k,i;
	//printf("In Compute_Mode_per_Genes\n");
//	printf("N = %i\n",N);
	for (k=0; k<N; k++)
	{
		
// 		if (Genes[k].values != NULL)
// 		{
			//printf("Genes[%i].matched = %lld\n",k,Genes[k].matched);
		// 	for (i=0; i<Genes[k].matched; i++)
// 			{
// 				printf("Genes[%i].values = %f\t",i, Genes[k].values[i]);
// 			}
// 			printf("\n");
			if (Genes[k].matched == 1)
			{
				//printf("Genes[%i].values[0] = %f\n",k, Genes[k].values[0]);
				Genes[k].mode = Genes[k].values[0];
			}
			else
			{
				double tmp;
				tmp = FloatMode(Genes[k].values,Genes[k].matched);
				if (tmp == 0.0) Genes[k].mode=Genes[k].values[0];
				else Genes[k].mode = tmp;
			}

			
// 		
// 		
// 			
// 		}
// 		else
// 		{
// 			printf("Genes[%i].values is NULL!\n",k);
// 		}
	}
	return;
}

void Compute_Min_Max_per_Genes(tLocus* Genes, int N)
{
	int k,i;
	//printf("Entering Compute_Min_Max_per_Genes\n");
	
	double swap;
	double min;
	double max;
	for (k=0; k<N; k++)
	{
		min = +200.0;
		max = -200.0;
		for (i=0; i<Genes[k].matched; i++)
		{
			if (Genes[k].values[i]<min) 
			{
				min=Genes[k].values[i];
			}
			if (Genes[k].values[i]>max)
			{
				max=Genes[k].values[i];
			} 
		}//now min and max are minimum and maximum value in values[i]
		Genes[k].min = min;
		Genes[k].max = max;
		// printf("Genes[%i].min = %f\n",k,Genes[k].min);
// 		printf("Genes[%i].max = %f\n",k,Genes[k].max);
	}//k
	//printf("Exiting Compute_Min_Max_per_Genes\n");
	return;
}

int FloatCompLocus(const void* elem1, const void* elem2)
{
	tLocus *Locus1 = (tLocus*) elem1;
	tLocus *Locus2 = (tLocus*) elem2;
	printf("%f\n",Locus1->value);
	return Locus1->value-Locus2->value;
}
 
//  int chr_num; 			//chromosome number in the genome
// 	char chr[15];			//chromosome number in string to include the case of X and Y chromosome
// 	long long start;		//start position in chromosome
// 	long long end;			// end position in chromosome
// 	char name[30]; //gene name
// 	double value;
// 	double min;
// 	double max;
// 	double mean;
// 	double mode;
// 	int omega;
// 	long long matched; // the number of values
// 	double* values;
void sort(tLocus* Loci, int NUMBER_OF_LOCI)
{

  
//    qsort(Loci, NUMBER_OF_LOCI, sizeof(tLocus), floatcomp);
	int i;
 
    // for(i = 0; i < 10; i++)
//         printf("%f\n", array[i]);
//     printf("\n");
 
    qsort(Loci, NUMBER_OF_LOCI, sizeof(tLocus), FloatCompLocus);
 
//     for(i = 0; i < 10; i++)
//         printf("%f\n", array[i]);
 return;
}
 void quick_sort(tLocus* Loci, int low, int high)
{
	int pivot, i,j;
	tLocus temp;
	
	if (low < high)
	{
		pivot = low;
		i = low;
		j = high;
		while (i < j)
		{
			while ((Loci[i].value <= Loci[pivot].value) && (i<high))
			{
				i++;
			}
			while (Loci[j].value>Loci[pivot].value)
			{
				j--;
			}
			if (i<j)
			{
				temp = Loci[i];
				Loci[i] = Loci[j];
				Loci[j] = temp;
			}
		}// while i<j
		
		temp = Loci[pivot];
		Loci[pivot] = Loci[j];
		Loci[j] = temp;
		quick_sort(Loci, low, j-1);
		quick_sort(Loci,j+1, high);
	}// if
	return;
}
int Use_Ref_Gene(char* file_name, tLocus* Loci, int N)
{
	
	FILE *file;

	if ( !(file=fopen(file_name,"r")))
	{
		printf("Cannot open file %s\n",file_name);
		exit(-1);
	}
	
	tLocus Locus;
	
	char chr[20];
	char name[20];

	long long start;
	long long end;
	int k;
	k=0;
	int j;
	int counted=0;
	int result;
		counted = 0;
	for (k=0;k<N; k++)
	{

//		printf("Loci[%i].chr = %s\n",k,Loci[k].chr);

    	while ( (result=fscanf(file,"%s %lld %lld %s",chr,&start,&end,name))!=EOF)
    	{
  //  		printf("result = %i\n",result);

//     		printf("chr = %s,\t start= %lld,\t end = %lld,\t, name = %s\n",chr,start,end,name);
    		strcpy(Locus.chr, chr);
//    		printf("counted=%i: Locus.chr = %s\n",counted,Locus.chr);
    		Locus.start = start;
    		Locus.end	= end;
    		strcpy(Locus.name,name);
  //  		printf("strcmp(Locus.chr,Loci[%i].chr)=%i\n",k,strcmp(Locus.chr,Loci[k].chr));
    	 	if (!strcmp(Locus.chr,Loci[k].chr))
 			{

	
 				if (LociOverLap(Locus, Loci[k]))
 				{
 				    		counted++;
 					//printf("Overlap:chr = %s,\t start= %lld,\t end = %lld,\t, name = %s\n",chr,start,end,name);
 					//printf("Locus.name  = %s\n",Locus.name);
 					strcpy(Loci[k].name, Locus.name);
 					 //printf("Loci[%i].name  = %s\n",k,Loci[k].name);
// 						values_buffer[matched] = ReadLoci[i].value;
// 						matched++;
// 						
// 						ReadLoci[i].omega=1;			
// 						strcpy(ReadLoci[i].name, SelectedLoci[k].name);
					break;
 				}

			}// if
    	}//while
    	rewind(file);
    	//printf("k=%i: counted = %i\n",k,counted);   
	}



	
//    while ( (fscanf(file,"%s %lf %lf %s",chr,&start,&end,name))!=EOF)
//    {
// 	   strcpy(SelectedGenes[k].chr,chr);
// 	   
// 	   strcpy(SelectedGenes[k].name,name);
// 
// 	   SelectedGenes[k].start=(long long) start;
// 
// 	   SelectedGenes[k].end=(long long ) end;
// 
// 	   SelectedGenes[k].min=0.0;
// 	   SelectedGenes[k].min=0.0;
// 	   k++;
// 	   if (k==NUM)
// 	   {
// 		   printf("Allocating more memory ... \n");
// 		   SelectedGenes=AllocateMoreMemToLoci(SelectedGenes, NUM,2*NUM);
// 		   NUM=2*NUM;
// 	   }
// 	   //printf("k=%i\n",k);
//    }
//    
//    *NUMBER_OF_SELECTED_GENES=k;
//    return SelectedGenes;
//    double* values_buffer = (double*) malloc (sizeof(double)*NUMBER_OF_LOCI);
//    int k,i;
//    for (k=0; k<NUMBER_OF_LOCI; k++) values_buffer[k]=0;
//    int num_matched, matched=0;
//    for (k=0; k<NUMBER_OF_SELECTED_LOCI;k++)
//    {
// 	   ReadLoci[k].omega=0;
//    }
//    double* weights = (double*) malloc(sizeof(double)*NUMBER_OF_LOCI);
//    for (k=0; k<NUMBER_OF_LOCI; k++) weights[k]=0;
// 
// 
//    for (k=0;k<NUMBER_OF_SELECTED_LOCI; k++)
//    {
// 	   for (i=0; i<NUMBER_OF_LOCI; i++)
// 	   {
// 
// 		   if (!strcmp(SelectedLoci[k].chr,ReadLoci[i].chr))
// 		   {
// 		   
// 			   if (LociOverLap(ReadLoci[i], SelectedLoci[k]))
// 			   {
// 					   
// 					   values_buffer[matched] = ReadLoci[i].value;
// 					   matched++;
// 					   
// 					   ReadLoci[i].omega=1;			
// 					   strcpy(ReadLoci[i].name, SelectedLoci[k].name);
// 			   }
// 
// 		   }
// 	   }//for i
//    
// 	   num_matched=num_matched+matched;
// 	   //printf("[%i]: mathched = %i\n",k,matched);
// 	   SelectedLoci[k].values = (double*) malloc(sizeof(double)*matched);
// //		printf("(SelectedLoci[k].values == NULL) = %i\n",SelectedLoci[k].values==NULL);
// 	   SelectedLoci[k].matched=matched;
// 	   for (i=0; i<matched; i++) 		
// 		   SelectedLoci[k].values[i]=values_buffer[i];
// 	   matched=0;
// 	   //free(values_buffer);
//    }//for k
// 
//    printf("I have compared %i genes\n",NUMBER_OF_SELECTED_LOCI);
//    printf("There are %i values identified in %i read  genes\n",num_matched, NUMBER_OF_SELECTED_LOCI);
//    gene_values_num = num_matched;
//    printf("There are %i values throughout the bedfile\n",NUMBER_OF_LOCI);
	return counted;
}

void get_monotone_regions(tLocus* Loci, int N)
{
	int i;
	Loci[0].sign=-1;//arbitrary choice
	for (i=0; i<(N-1); i++)
	{
		if (Loci[i].value>=Loci[i+1].value)
		{
			Loci[i].sign =-1;
		}
		else
		{
			Loci[i].sign = +1;
		}
	}
}
int get_size_of_monotone_regions(tLocus* Loci, int N, int* decrease, int* increase, char* parameter)
{
	int tmp_decrease = 0, tmp_increase = 0;
	int k;
	int decrease_length=0;
	int increase_length = 0;
	// for (k=0; k<N; k++)
// 	{
// 		printf("Loci[%i].name = %s,\t Loci[%i].omega = %i\n",k, Loci[k].name,k, Loci[k].omega);
// 	}
	//printf("GET_SIZE_OF_MONOTONE_REGIONS\n");
//	getchar();
	k=0;
	int counter=0;
	// if (strcmp(parameter, "Genes")==0)  printf("Genes\n");
// 	if (strcmp(parameter, "Loci")==0)  printf("Loci\n");
	int decrease_length_array_size = 100;
	int* decrease_length_array = (int*) malloc(sizeof(int)*decrease_length_array_size);
	
	int switch_register=1;
	if (strcmp(parameter,"Genes")==0)
	{
		while (k<N)
		{
			// printf("iteration k = %i\n",k);
			while (Loci[k].omega == 0)
			{
				if (k<N)
				{
					k++;
				}
				else
				{
					break;
				}
			}
			// printf("after k=%i\n",k);
			if ((switch_register) && (Loci[k].omega == 1))
			{
				while ( (Loci[k].sign == -1) && (Loci[k].omega == 1))
				{ 
					decrease_length_array[counter]++;
					// printf("Loci[%i].sign = %i\n",k, Loci[k].sign);
					k++;
				}
// 				printf("1253: k =%i\n",k);
				switch_register = 0;
			}
			else // either end of decrease region or end of Genes region
			{
				// printf("else branch\n");
// 				printf("switch_register = %i\n",switch_register);
// 				printf("counter = %i\n",counter);
// 				printf("Loci[%i].omega = %i\n",k, Loci[k].omega);
// 				printf("decrease_length_array[%i]=%i\n",counter, decrease_length_array[counter]);
				if (decrease_length_array[counter]>0) counter++;
				k++;
				if (counter==decrease_length_array_size) 
				{
					int * result;
					//printf("if (counter==decrease_length_array_size)\n");
				//	getchar();
					result = AllocateMoreMemToIntArray(decrease_length_array,decrease_length_array_size,2*decrease_length_array_size);
					if (result) decrease_length_array_size *= 2; else 
					{
						printf("Could not allocate enough memory\n");
						return -1;
					}
					decrease_length_array = result;
				}
				
				switch_register = 1;
			}
		}
	}
	else if (strcmp(parameter, "Loci")==0)
	{
		
		while (k<N)
		{
			
			// printf("Loci:iteration k = %i\n",k);
			while (Loci[k].omega == 1)
			{
				if (k<N)
				{
					k++;
				}
				else
				{
					break;
				}
			}
// 			printf("after k=%i\n",k);
			if ((switch_register) && (Loci[k].omega == 0))
			{
				while ( (Loci[k].sign == -1) && (Loci[k].omega == 0))
				{ 
					decrease_length_array[counter]++;
// 					printf("Loci:Loci[%i].sign = %i\n",k, Loci[k].sign);
					k++;
				}
// 				printf("1253: Loci:k =%i\n",k);
				switch_register = 0;
			}
			else // either end of decrease region or end of Genes region
			{
// 				printf("Loci:else branch\n");
// 				printf("Loci:switch_register = %i\n",switch_register);
// 				printf("Loci:counter = %i\n",counter);
// 				printf("Loci:Loci[%i].omega = %i\n",k, Loci[k].omega);
// 				printf("decrease_length_array[%i]=%i\n",counter, decrease_length_array[counter]);
				if (decrease_length_array[counter]>0) counter++;
				k++;
				if (counter==decrease_length_array_size) 
				{
					int * result;
					//printf("if (counter==decrease_length_array_size)\n");
					//getchar();
					result = AllocateMoreMemToIntArray(decrease_length_array,decrease_length_array_size,2*decrease_length_array_size);
					if (result) decrease_length_array_size *= 2; else 
					{
						printf("Could not allocate enough memory\n");
						return -1;
					}
					decrease_length_array = result;
				}
				
				switch_register = 1;
			}
		}//for
	}
	int i;
// 	if (strcmp(parameter, "Genes")==0)  printf("Genes\n");
// 	if (strcmp(parameter, "Loci")==0)  printf("Loci\n");
	// printf("counter = %i\n",counter);
// 	for (i=0; i<counter; i++)
// 	{
// 		printf("decrease_length_array[%i] = %i\n",i, decrease_length_array[i]);
// 	}
// 	
	double decrease_mean=0.0;
	int sensitive_areas=0;
	for (i=0; i<counter; i++)
	{
		if (decrease_length_array[i]>SENSITIVITY)
		{
			decrease_mean = +decrease_length_array[i];
			sensitive_areas++;
		}
	}
	if (sensitive_areas>0)
	decrease_mean = decrease_mean /(1.0*sensitive_areas);
	else decrease_mean = 0.0;
	// while (k<N)
// 	{
// 		//if (counter>20) break; else counter++;
// 		printf("k = %i\n",k);
// 		if (strcmp(parameter,"Genes")==0)
// 		{	
// 			
// 			printf("Loci[%i].omega = %i,\t Loci[%i].name = %s\n",k, Loci[k].omega, k, Loci[k].name);
// 			if (Loci[k].omega == 1)
// 			{	
// 				printf("if (Loci[k].omega == 1)\n");	
// 				decrease_length = 0;
// 				increase_length = 0;
// 				while (Loci[k].sign == -1) 
// 				{
// 					printf("while (Loci[k].sign == -1)\n");
// 					k++;
// 					decrease_length++;
// 				}
// 				printf("decrease_length = %i\n",decrease_length);
// 
// 				while (Loci[k].sign==+1) 
// 				{					
// 					k++;
// 					tmp_increase++;
// 				}
// 				printf("increase_length = %i\n",increase_length);
// 				
// 			}// if
// 			else k++;
// 		}	
// 		else if (strcmp(parameter, "Loci")==0)
// 		{
// 			if (Loci[k].omega == 0)
// 			{
// 				decrease_length = 0;
// 				increase_length = 0;
// 				while (Loci[k].sign == -1) 
// 				{
// 					k++;
// 					decrease_length++;
// 				}
// 				printf("decrease_length = %i\n",decrease_length);
// 
// 				while (Loci[k].sign==+1) 
// 				{
// 					k++;
// 					tmp_increase++;
// 				}
// 				printf("increase_length = %i\n",increase_length);
// 			}// if
// 			else k++;
// 		}
// 	}
	*decrease = tmp_decrease;
	*increase = tmp_increase;
	return decrease_mean;
}
void get_extremum_values(tLocus* Loci, int N, double* min, double* max, char* parameter)
{
	int i;
	
	
	*min = +200.0;
	*max = -200.0;
	
	for (i=0; i<N; i++)
	{
		if (strcmp(parameter,"Loci")==0)
		{
			if (Loci[i].omega==0)
			{
//				printf("@@@@^^^^\n");
				if (Loci[i].value < *min) 
				{
//					printf("@@@@^^^^\n");
					*min = Loci[i].value;
				}
				if (Loci[i].value > *max) 
				{
//					printf("@@@@^^^^\n");
					*max = Loci[i].value;
				}
			}
		}
		else if (strcmp(parameter,"Genes")==0)
		{
			if (Loci[i].omega==1)
			{
				if (Loci[i].value < *min) *min = Loci[i].value;
				if (Loci[i].value > *max) *max = Loci[i].value;
			}
		}
	}//	for (i=0; i<N; i++)
		
	return;
}
int main(int argc, char** argv)
{
	//  argv[1] loci  file
	//	argv[2]	genes file
	//  argv[3] reference genome
	
	char* 	genes_file			=	argv[1];
	char*	selected_genes_file	=	argv[2];

	int NUMBER_OF_LOCI;

	int NUMBER_OF_SELECTED_GENES;
	
	tLocus* Loci;
	tLocus* SelectedGenes; 
	
	
	Loci	=	ReadLocusValues(genes_file, &NUMBER_OF_LOCI);
 	get_monotone_regions(Loci,NUMBER_OF_LOCI);	
	SelectedGenes = ReadGenes(selected_genes_file,&NUMBER_OF_SELECTED_GENES);


	
	
	//Matching values from Loci with SelectedGenes
	MatchValues(Loci,NUMBER_OF_LOCI,SelectedGenes,NUMBER_OF_SELECTED_GENES);
	
	
//	printf("Printing loci ...\n");
	//PrintLocuss(Loci, NUMBER_OF_LOCI,"Loci");
	
	
	char* reference_genome = argv[3];
	//"hg19_RefGene.txt";
	
	
	printf("Looking missing genes in %s ...\n",reference_genome);	
	int more_genes_found;
	more_genes_found = Use_Ref_Gene("hg19_RefGene_test.txt", Loci, NUMBER_OF_LOCI);
	printf("%i new genes have been identified\n",more_genes_found);
	printf("Printing updated loci\n");
	 	int increase, decrease;
  printf("Calculating the size of decrease regions with length > %i in Genes ...\n",SENSITIVITY);
  	double decrease_mean;
	decrease_mean = get_size_of_monotone_regions(Loci, NUMBER_OF_LOCI, &decrease, &increase, "Genes");
	printf("Mean length of decrease regions with length > %i in Genes = %f\n",SENSITIVITY,decrease_mean);
//	printf("Genes: Decrease region size = %i,\t Increase region size = %i\n",decrease,increase);
	printf("Calculating the size of decrease regions with length >%i in Loci without Genes ...\n",SENSITIVITY);
	decrease_mean = get_size_of_monotone_regions(Loci, NUMBER_OF_LOCI, &decrease, &increase, "Loci");
	printf("Mean length of decrease regions with length >%i in Loci without Genes = %f\n",SENSITIVITY,decrease_mean);
//	get_size_of_monotone_regions(Loci, NUMBER_OF_LOCI, &decrease, &increase, "Loci");
//	printf("Loci: Decrease region size = %i,\t Increase region size = %i\n",decrease,increase);
  // 	PrintLocuss(Loci, NUMBER_OF_LOCI,"Loci");

	
// Minimum and maximum over all loci
//	MinMax(Loci, NUMBER_OF_LOCI, SelectedGenes, NUMBER_OF_SELECTED_GENES); 

	// int i, length; 
// 
 	Compute_Mode_per_Genes(SelectedGenes, NUMBER_OF_SELECTED_GENES);
 	Compute_Min_Max_per_Genes(SelectedGenes, NUMBER_OF_SELECTED_GENES);
	//Compute_Mean_per_Genes(SelectedGenes, NUMBER_OF_SELECTED_GENES);
// 	Compute_SignsDispersion_per_Genes(SelectedGenes, NUMBER_OF_SELECTED_GENES);
  	int i;
  	double min, max;
  	for (i=0; i<NUMBER_OF_SELECTED_GENES; i++)
 	{
		printf("Genes[%i]:\t name = %8s,\t min = %f,\t",i, SelectedGenes[i].name, SelectedGenes[i].min); 
		printf("mean = %f,\t mode = %f,\t", SelectedGenes[i].mean,SelectedGenes[i].mode);
		printf("max = %f,\t dispersion = %f\n", SelectedGenes[i].max, SelectedGenes[i].dispersion);
 	}	
// // 	tRegion Region = FindMinWindowDispersion(Loci, 0, NUMBER_OF_LOCI, 5);
// // 	printf("Region.start  = %lld,\t Region.end = %lld, Region.Dispersion = %f\n",Region.start, Region.end, Region.dispersion);
// 	
// 	
// 
// 	double* signs = (double*) malloc(sizeof(double)*NUMBER_OF_LOCI);
// 	for (i=0; i<NUMBER_OF_LOCI; i++)
// 	{
// 	// 	printf("[%3i]:gene=%10s, chr=%5s, start=%12lld, end=%12lld, value = %1.8lf\n",i,Loci[i].name, Loci[i].chr, Loci[i].start, Loci[i].end,Loci[i].value);
// 	//	printf("sign[%i]=%i\n",i, Loci[i].sign);
// 		signs[i] = Loci[i].sign;
// 	}
// 	tRegion Region = FindMinWindowDispersionArray(signs, 0, NUMBER_OF_LOCI, 6);
// 	
// 	printf("Printing genes ...\n");
// 	
// 	//PrintLocuss(SelectedGenes, NUMBER_OF_SELECTED_GENES,"Genes");
// 	
// 	printf("Printing loci ...\n");
// 	PrintLocuss(Loci, NUMBER_OF_LOCI,"Loci");
// 	
//  	
  	double LociMeanValue;
//  	double LociModeValue;
// 	double min, max;
 	double GenesMeanValue; 
//  	double GenesModeValue; 
// 
double* weights;
  	w2 = w1 = 1.0;
  	weights = init_weights(NUMBER_OF_LOCI,Loci,w1,w2);
   	GenesMeanValue = MeanValue(Loci,  NUMBER_OF_LOCI, weights,"Genes");
   	LociMeanValue = MeanValue(Loci,  NUMBER_OF_LOCI, weights,"Loci");
//   	GenesModeValue = Mode(Loci,  NUMBER_OF_LOCI,"Genes");
//  	printf("\n");	
   	printf("Genes MeanValue=%f\n", GenesMeanValue);
   	printf("Loci MeanValue=%f\n", LociMeanValue);
//  	printf("\n");	
//  	printf("\n");	
//   //	printf("Genes ModeValue=%f\n",	GenesModeValue);
//   	printf("\n");	
//   	
   	printf("Finding minimum and maximun throughout gene values ...\n");
   	get_extremum_values(Loci, NUMBER_OF_LOCI, &min, &max, "Genes");
   	printf("Genes region: min = %f,\t max= %f\n",min,max);
   	printf("Finding minimum and maximun throughout loci values excluding gene values ...\n");
   	get_extremum_values(Loci, NUMBER_OF_LOCI, &min, &max, "Loci");
   	printf("Loci without Genes regions: min = %f,\t max= %f\n",min,max);
// 	
// 	for (i=0; i<NUMBER_OF_LOCI; i++)
// 	{
// 		printf("Loci[%i].sign = %i\n",i, Loci[i].sign);
// 	}
// 	
	return 0;
}
	
