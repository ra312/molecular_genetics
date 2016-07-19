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
	char sex;
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
	int* decrease_regions;
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

typedef struct
{
	int N; // number of Loci
	
	tLocus* Loci;
	
	char sex; // X,Y, or 0.

} tChrom;
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
	
	free(array);
	//printf("Inside  int* AllocateMoreMemToIntArray(int* array, int old_size, int new_size)\n");
	// for (i=0; i<old_size; i++)
// 	{
// 		printf("new_array[%i]=%i\n",i, new_array[i]);
// 	}
	return new_array;
}
// tChrom* MoreChrom(tChrom* array, int old_size, int new_size)
// {
// 	tChrom* new_array = (tChrom*)malloc(sizeof(tChrom)*new_size);
// 	int i;
// 	int j;
// 	int k;
// 	
// 	if (new_size<old_size)
// 	{
// 		printf("Please enter bigger array size\n");
// 		return NULL;
// 	}
// 	for (i=0; i<old_size; i++)
// 	{
// 		new_array[i].sex=array[i].sex;
// 		for (k=0; k<array[i].N; k++)
// 		{
// 		new_array[i].Loci[k].values = (int*) malloc(sizeof(int)*array[i].Loci[k].matched);
// 		new_array[i].Loci[k].signs = (int*) malloc(sizeof(int)*array[i].Loci[k].matched);
// 		new_array[i].Loci[k].decrease_regions = (int*) malloc(sizeof(int)*array[i].Loci[k].matched);
// 		
// 		for (j=0; j<array[i].matched; i++)
// 		{
// 			new_array[i].Loci[k]values[j]=array[i].values[j];
// 			new_array[i].signs[j]=array[i].signs[j];
// 			new_array[i].decrease_regions[j]=array[i].decrease_regions[j];
// 		}
// 			free(array[i].Loci[k].values);
// 			free(array[i].Loci[k].signs);
// 			free(array[i].Loci[k].decrease_regions);
// 
// 		}
// 	}
// 	free(array);
// 	//printf("Inside  int* AllocateMoreMemToIntArray(int* array, int old_size, int new_size)\n");
// 	// for (i=0; i<old_size; i++)
// // 	{
// // 		printf("new_array[%i]=%i\n",i, new_array[i]);
// // 	}
// 	return new_array;
// }
tChrom* ReadLocusValues(char* genes_file)
{
	int NUM = 30000;
	
	char* chr;
	char* chr_num;
	char sex;
	chr = (char*) malloc(sizeof(char)*5);
	tChrom* Chrom = (tChrom*) malloc(sizeof(tChrom)*23);
	int k,i;
	for (k=0; k<23; k++) 
	{
		Chrom[k].Loci = (tLocus*) malloc(sizeof(tLocus)*NUM);
		Chrom[k].N = NUM;
	}

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
	
	k=0;
	int Loci_Number_at_Chr[23];
	for (k=0; k<23; k++) 
	{
		Loci_Number_at_Chr[k]=0;
		if ( (Chrom[k].Loci = (tLocus*) malloc(sizeof(tLocus)*NUM))==NULL)
		{
			printf("Memory could not be allocated\n");
		}
		printf("testing\n");
		for (i=0; i<NUM; i++) Chrom[k].Loci[i].value =0;
		printf("testing\n");
	}
	int LociNumber;
	//reading Loci and saving them to the correct chromosome
	while ( (fscanf(file,"%s %lld %lld %lf",chr,&start,&end,&value))!=EOF)
	{
		//chrom_num=atoi(chr+3);
printf("192££££\n");
		if (strstr(chr,"Y")) 
		{
			sex = 'Y';
			chrom_num=23;
		}
		else if (strstr(chr,"X"))
		{
			sex = 'X';
			chrom_num=23;
		}
		else 
		{
			sex = 0; // no sex;
			chrom_num=atoi(chr+3);
		
		}
		//printf("Catching seg fault\n");
		printf("210:Loci_Number_at_Chr[%i]=%i\n",chrom_num, Loci_Number_at_Chr[chrom_num]);
//		strcpy(Loci[k].chr,chr);
//		Loci[k].chr_num=chrom_num;
		LociNumber = Loci_Number_at_Chr[chrom_num];
		printf("214L LociNumber = %i\n",LociNumber);
		printf("chrom_num = %i\n",chrom_num);
		printf("start = %i\n",start);
		Chrom[chrom_num].Loci[LociNumber].start=start;
				printf("216:Catching seg fault\n");
//		Chrom[chrom_num].Loci[LociNumber].end = end;

//		Chrom[chrom_num].sex=sex;
//		Chrom[chrom_num].Loci[LociNumber].value = value;		
		Loci_Number_at_Chr[chrom_num]++;
		LociNumber = Loci_Number_at_Chr[chrom_num];
		printf("220:LociNumber = %i\n",LociNumber);
				printf("Catching seg fault\n");
		if (LociNumber == NUM)
		{
			printf("Allocating more memory ... \n");
//			int NUMBER = LociNumber;
			Chrom[i].Loci=AllocateMoreMemToLoci(Chrom[chrom_num].Loci ,LociNumber, 2*LociNumber);
			if (Chrom[i].Loci)
			{
				Loci_Number_at_Chr[chrom_num] = LociNumber*2;
			}
			else
			{
				printf("230: Memory could not be allocated\n");
			}
		}
		
	}
//	*NUMBER_OF_LOCI=k;
	
	fclose(file);
	return Chrom;
}
void PrintChromosome(tChrom* Chrom, int chromosome_number)
{
	int i,j;
	int k = chromosome_number;
	printf("Chromosome %i\n",k);
	//Chrom[k].N = number of Loci in N-th chromosome
	//Chrom[k].Loci[i] = i-th Locus with Chrom[k].Loci[i].matched values
	for (i=0; i<Chrom[k].N; i++)
	{

		for (j=0; j<Chrom[k].Loci[i].matched; j++)
		{
			printf("gene = %5s", Chrom[k].Loci[i].name);
			printf("start = %10lld", Chrom[k].Loci[i].start);
			printf("end = %10lld", Chrom[k].Loci[i].end);
			printf(" value = 1.7%f\n",Chrom[k].Loci[i].value);			
		}
	}
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
	tChrom* Genes = (tChrom*) malloc(sizeof(tChrom)*23);
	int chr_index[23];
	for (k=0; k<23; k++) chr_index[k]=0;
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
// void MinMax(tLocus* ReadLoci, int NUMBER_OF_LOCI, tLocus* SelectedLoci, int NUMBER_OF_SELECTED_LOCI)
// {
// 	//MinMax computes MIN and MAX per each gene in SelectedLoci
// 	int k,j;
// 	double swap;
// 	double min;
// 	double max;
// 	for (k=0; k<NUMBER_OF_SELECTED_LOCI; k++)
// 	{
// 		min=+200.0;
// 		max=-200.0;
// 		for (j=0; j<NUMBER_OF_LOCI; j++)
// 		{
// 			if (!strcmp(SelectedLoci[k].name,ReadLoci[j].name))
// 			{
// 				
// 				if (ReadLoci[j].value<min) 
// 				{
// //					printf("Min\n");
// 					min=ReadLoci[j].value;
// 				}
// 				if (ReadLoci[j].value>max)
// 				{
// 							max=ReadLoci[j].value;
// 				} 
// 			}
// 		}//for j
// 		for (j=0; j<NUMBER_OF_LOCI; j++)
// 			if (!strcmp(SelectedLoci[k].name,ReadLoci[j].name))
// 			{
// 				ReadLoci[j].min=min;
// 				ReadLoci[j].max=max;
// 			}
// 	
// 	}//for k
// 			
// 	
// 
// 	return;
// }
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
						printf("LociOverlap: i = %i\n",i);

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
    	while ( (result=fscanf(file,"%s %lld %lld %s",chr,&start,&end,name))!=EOF)
    	{
    		strcpy(Locus.chr, chr);
    		Locus.start = start;
    		Locus.end	= end;
    		strcpy(Locus.name,name);
    	 	if (!strcmp(Locus.chr,Loci[k].chr))
 			{
 				if (LociOverLap(Locus, Loci[k]))
 				{
 				    		counted++;
 					strcpy(Loci[k].name, Locus.name);
 					break;
 				}

			}// if
    	}//while
    	rewind(file);
    	//printf("k=%i: counted = %i\n",k,counted);   
	}
	return counted;
}

void get_monotone_regions(tLocus* Loci, int N, char* parameter)
{
	int i, error;

	for (i=0; i<N; i++)
	{
		if (Loci[i].value<= 0)
		{
			Loci[i].sign =-1;
		}
		else
		{
			Loci[i].sign = +1;
		}
	}
	i=0;
	int j;
	int start, length;
	if (strcmp(parameter, "Genes")==0) 
	{
		error = 15;
		for (j=0; j<N; j++)
		{
			while (i<Loci[j].matched)
			{
				while (Loci[j].signs[i++] == -1) ;
	
				start = i;
				length = 0;
				while (Loci[j].signs[i++] == +1) length++;
				if (length <error)
				{
					for (i=start; i<length; i++)
					{
						Loci[j].signs[i] = -1;
					}
				}
			}
			
		}
	}
	if (strcmp(parameter,"Loci")==0) 
	{
		error = 2;
		while (i<N)
		{
			while (Loci[i++].sign == -1) ;
			start = i;
			length = 0;
			while (Loci[i++].sign == +1) length++;
			if (length <error)
			{
				for (i=start; i<length; i++)
				{
					Loci[i].sign = -1;
				}
			}
		}
	}
	return;
}
double get_size_of_monotone_regions(tLocus* Loci, int N, int* decrease, int* increase, char* parameter)
{
	int tmp_decrease = 0, tmp_increase = 0;
	int k,i;
	int decrease_length=0;
	int increase_length = 0;
	k=0;
	int counter=0;
	int decrease_length_array_size = 100;
	int* decrease_length_array = (int*) malloc(sizeof(int)*decrease_length_array_size);
	
 	for (k=0; k<decrease_length_array_size; k++)
 	{
 		decrease_length_array[k]=0;
 	}
	int switch_register=1;

	double decrease_mean=0.0;
 	int sensitive_areas=0;

  	if (strcmp(parameter, "Loci")==0)	
  	{
		while (k<N)
		{				
			int tmp_start = k, tmp_length;
			while (Loci[k].omega == 1) if (k<N) k++; else break;
			if ((switch_register) && (Loci[k].omega == 0))
			{
				while ( (Loci[k].sign == -1) && (Loci[k].omega == 0))
				{ 
					decrease_length_array[counter]++;
					k++;
				}
				switch_register = 0;
			}
			else // either end of decrease region or end of Genes region
			{
				if (decrease_length_array[counter]>0) counter++;
				k++;
				if (counter==decrease_length_array_size) 
				{
					int * result;
					result = AllocateMoreMemToIntArray(decrease_length_array,decrease_length_array_size,2*decrease_length_array_size);
					if (result) decrease_length_array_size *= 2; else 
					{
						printf("Could not allocate enough memory\n");
						return -1;
					}
					decrease_length_array = result;
				}
				switch_register = 1;
			}//else
		}//while
	}
	else if (strcmp(parameter,"Genes")==0)
	{	
		while (k<N)
 		{
 			int tmp_start = k, tmp_length;
			while (Loci[k].omega == 1) if (k<N) k++; else break;		
			if ((switch_register) && (Loci[k].omega == 0))
			{
				while ( (Loci[k].sign == -1) && (Loci[k].omega == 0))
				{ 
					decrease_length_array[counter]++;
					k++;
				}
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
		}//while
	}
	
	printf("!!!!!End get_monotone_regions: decrease_mean = %f\n",decrease_mean);
// 	*decrease = tmp_decrease;
// 	*increase = tmp_increase;
	if (strcmp(parameter,"Genes")==0)
	{
	for (i=0; i<counter; i++)
	{
		if (decrease_length_array[i]>=1)
		{
			printf("Decrease_lengt_array[%i]=%i\n",i, decrease_length_array[i]);
			decrease_mean += ((double)decrease_length_array[i]);
		//	printf("Loop i = %i, decrease_mean = %f\n",i, decrease_mean);
			sensitive_areas++;
		}
	}
	}
	printf("decrease_mean = %f\n",decrease_mean);
	if (sensitive_areas>0)
	decrease_mean = decrease_mean /(1.0*sensitive_areas);
	else decrease_mean = 0.0;
	printf("decrease_mean = %f\n",decrease_mean);
	free(decrease_length_array);
	return decrease_mean;
}

// void get_extremum_values(tLocus* Loci, int N, double* min, double* max, char* parameter)
// {
// 	int i;
// 	
// 	
// 	*min = +200.0;
// 	*max = -200.0;
// 	
// 	for (i=0; i<N; i++)
// 	{
// 		if (strcmp(parameter,"Loci")==0)
// 		{
// 			if (Loci[i].omega==0)
// 			{
// //				printf("@@@@^^^^\n");
// 				if (Loci[i].value < *min) 
// 				{
// //					printf("@@@@^^^^\n");
// 					*min = Loci[i].value;
// 				}
// 				if (Loci[i].value > *max) 
// 				{
// //					printf("@@@@^^^^\n");
// 					*max = Loci[i].value;
// 				}
// 			}
// 		}
// 		else if (strcmp(parameter,"Genes")==0)
// 		{
// 			if (Loci[i].omega==1)
// 			{
// 				if (Loci[i].value < *min) *min = Loci[i].value;
// 				if (Loci[i].value > *max) *max = Loci[i].value;
// 			}
// 		}
// 	}//	for (i=0; i<N; i++)
// 		
// 	return;
// }
void DeleteLoci(tLocus* Loci, int N)
{
	
	int i;
	
	for (i=0; i<N; i++)
	{
		if (Loci[i].values) free(Loci[i].values);
		if (Loci[i].signs) free(Loci[i].signs);
	}
	free(Loci);
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
	
	tChrom* Chrom=	ReadLocusValues(genes_file);
	//PrintChromosome(Chrom, 5);
	
 	// SelectedGenes = ReadGenes(selected_genes_file,&NUMBER_OF_SELECTED_GENES);
// 
//   	get_monotone_regions(Loci,NUMBER_OF_LOCI,"Loci");
//   	get_monotone_regions(SelectedGenes,NUMBER_OF_SELECTED_GENES,"Genes");

	
	
	//Matching values from Loci with SelectedGenes
// 	 MatchValues(Loci,NUMBER_OF_LOCI,SelectedGenes,NUMBER_OF_SELECTED_GENES);

// 	 char* reference_genome = argv[3];
	//"hg19_RefGene.txt";
	
	// 
//  	printf("Looking missing genes in %s ...\n",reference_genome);	
//  	int more_genes_found;
//  	more_genes_found = Use_Ref_Gene("hg19_RefGene_test.txt", Loci, NUMBER_OF_LOCI);
//  	printf("%i new genes have been identified\n",more_genes_found);
// 	printf("Printing updated loci\n");
// 
//   printf("Calculating the size of decrease regions with length > %i in Genes ...\n",SENSITIVITY);
// 
//  	decrease_mean = get_size_of_monotone_regions(Loci, NUMBER_OF_LOCI, &decrease, &increase, "Genes");
//  	printf("Mean length of decrease regions with length > %i in Genes = %f\n",SENSITIVITY,decrease_mean);
// 	printf("Genes: Decrease region size = %i,\t Increase region size = %i\n",decrease,increase);
//  	printf("Calculating the size of decrease regions with length >%i in Loci without Genes ...\n",SENSITIVITY);
//  	decrease_mean = get_size_of_monotone_regions(Loci, NUMBER_OF_LOCI, &decrease, &increase, "Loci");
//  	printf("Mean length of decrease regions with length >%i in Loci without Genes = %f\n",SENSITIVITY,decrease_mean);
// //	get_size_of_monotone_regions(Loci, NUMBER_OF_LOCI, &decrease, &increase, "Loci");
// //	printf("Loci: Decrease region size = %i,\t Increase region size = %i\n",decrease,increase);
//   // 	PrintLocuss(Loci, NUMBER_OF_LOCI,"Loci");
// 
// 	
// // Minimum and maximum over all loci
// //	MinMax(Loci, NUMBER_OF_LOCI, SelectedGenes, NUMBER_OF_SELECTED_GENES); 
// 
// 	// int i, length; 
// // 
//   	Compute_Mode_per_Genes(SelectedGenes, NUMBER_OF_SELECTED_GENES);
//   	Compute_Min_Max_per_Genes(SelectedGenes, NUMBER_OF_SELECTED_GENES);
//  	Compute_Mean_per_Genes(SelectedGenes, NUMBER_OF_SELECTED_GENES);
// // 	Compute_SignsDispersion_per_Genes(SelectedGenes, NUMBER_OF_SELECTED_GENES);
//   	int i;
//   	double min, max;
//   	for (i=0; i<NUMBER_OF_SELECTED_GENES; i++)
//  	{
// 		printf("Genes[%i]:\t name = %8s,\t min = %f,\t",i, SelectedGenes[i].name, SelectedGenes[i].min); 
// 		printf("mean = %f,\t mode = %f,\t", SelectedGenes[i].mean,SelectedGenes[i].mode);
// 		printf("max = %f,\t dispersion = %f\n", SelectedGenes[i].max, SelectedGenes[i].dispersion);
//  	}	
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
  	// double LociMeanValue;
// //  	double LociModeValue;
// // 	double min, max;
//  	double GenesMeanValue; 
// //  	double GenesModeValue; 
// // 
// 	double* weights;
//   	w2 = w1 = 1.0;
//   	weights = init_weights(NUMBER_OF_LOCI,Loci,w1,w2);
//    	GenesMeanValue = MeanValue(Loci,  NUMBER_OF_LOCI, weights,"Genes");
//    	LociMeanValue = MeanValue(Loci,  NUMBER_OF_LOCI, weights,"Loci");
// //   	GenesModeValue = Mode(Loci,  NUMBER_OF_LOCI,"Genes");
// //  	printf("\n");	
//    	printf("Genes MeanValue=%f\n", GenesMeanValue);
//    	printf("Loci MeanValue=%f\n", LociMeanValue);
// //  	printf("\n");	
// //  	printf("\n");	
// //   //	printf("Genes ModeValue=%f\n",	GenesModeValue);
// //   	printf("\n");	
// //   	
//    	printf("Finding minimum and maximun throughout gene values ...\n");
//    	get_extremum_values(Loci, NUMBER_OF_LOCI, &min, &max, "Genes");
//    	printf("Genes region: min = %f,\t max= %f\n",min,max);
//    	printf("Finding minimum and maximun throughout loci values excluding gene values ...\n");
//    	get_extremum_values(Loci, NUMBER_OF_LOCI, &min, &max, "Loci");
//    	printf("Loci without Genes regions: min = %f,\t max= %f\n",min,max);
// // 	
// // 	for (i=0; i<NUMBER_OF_LOCI; i++)
// // 	{
// // 		printf("Loci[%i].sign = %i\n",i, Loci[i].sign);
// // 	}
// 	
	DeleteLoci(Loci, NUMBER_OF_LOCI);
	DeleteLoci(SelectedGenes, NUMBER_OF_SELECTED_GENES);
	return 0;
}
	
