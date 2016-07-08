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
} tLocus;


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
			arithm_mean	+= weights[k]*fabs(Loci[k].value)*Loci[k].omega;
			volume += weights[k]*Loci[k].omega;
		}
		arithm_mean = arithm_mean / volume; // conditional expectation.
// 		printf("Genes: arithm_mean = %f\n",arithm_mean);
		for (k=0; k<NUMBER_OF_LOCI; k++)
		{
			// printf("Loci[%i].value=%f\n", k, Loci[k].value);
// 			printf("Loci[%i].omega=%i\n",k,Loci[k].omega);
// 			printf("weights[%i]=%f\n",k, weights[k]);
			quadr_mean 	+= pow(0.0+Loci[k].value-arithm_mean, 2.0)*Loci[k].omega*weights[k];
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
			arithm_mean	+= Loci[k].value*weights[k]*(1-0*Loci[k].omega);
			volume		+= weights[k]*(1-0*Loci[k].omega);
		}
		arithm_mean = arithm_mean / volume;
		for (k=0; k<NUMBER_OF_LOCI; k++)
		{
			quadr_mean +=pow(0.0+Loci[k].value-arithm_mean, 2.0)*weights[k]*(1-0*Loci[k].omega);
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
	int k,i;
	for (k=0; k<NUMBER_OF_LOCI; k++) values_buffer[k]=0;
	int num_matched, matched=0;
	for (k=0; k<NUMBER_OF_SELECTED_LOCI;k++)
	{
		ReadLoci[k].omega=0;
	}
	double* weights = (double*) malloc(sizeof(double)*NUMBER_OF_LOCI);
	for (k=0; k<NUMBER_OF_LOCI; k++) weights[k]=0;


	for (k=0;k<NUMBER_OF_SELECTED_LOCI; k++)
	{
		for (i=0; i<NUMBER_OF_LOCI; i++)
		{

			if (!strcmp(SelectedLoci[k].chr,ReadLoci[i].chr))
			{
			
				if (LociOverLap(ReadLoci[i], SelectedLoci[k]))
				{
						
						values_buffer[matched] = ReadLoci[i].value;
						matched++;
						
						ReadLoci[i].omega=1;			
						strcpy(ReadLoci[i].name, SelectedLoci[k].name);
				}

			}
		}//for i
	
		num_matched=num_matched+matched;
		//printf("[%i]: mathched = %i\n",k,matched);
		SelectedLoci[k].values = (double*) malloc(sizeof(double)*matched);
//		printf("(SelectedLoci[k].values == NULL) = %i\n",SelectedLoci[k].values==NULL);
		SelectedLoci[k].matched=matched;
		for (i=0; i<matched; i++) 		
			SelectedLoci[k].values[i]=values_buffer[i];
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
int i, j,z, maxCount;
double tmp, modeValue;

int* tally=(int*)malloc(sizeof(int)*N);
  
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
void FindMinWindowDispersion(tLocus* Loci, int start, int end)
{
	return;
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
void Compute_Mode_per_Genes(tLocus* Genes, int N)
{
	
	int k,i;

	for (k=0; k<N; k++)
	{
		
		if (Genes[k].values != NULL)
		{
			Genes[k].mode=FloatMode(Genes[k].values,Genes[k].matched);
		
		
			
		}
		else
		{
			printf("Genes[%i].values is NULL!\n",k);
		}
	}
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
int main(int argc, char** argv)
{

	char* 	genes_file			=	argv[1];
	char*	selected_genes_file	=	argv[2];

	int NUMBER_OF_LOCI;

	int NUMBER_OF_SELECTED_GENES;
	
	tLocus* Loci;
	tLocus* SelectedGenes; 
	
	
	Loci	=	ReadLocusValues(genes_file, &NUMBER_OF_LOCI);
	
	SelectedGenes = ReadGenes(selected_genes_file,&NUMBER_OF_SELECTED_GENES);
	Compute_Mode_per_Genes(SelectedGenes, NUMBER_OF_SELECTED_GENES);

	//quick_sort(Loci, 0, NUMBER_OF_LOCI);
	MatchValues(Loci,NUMBER_OF_LOCI,SelectedGenes,NUMBER_OF_SELECTED_GENES);
	//printf("NUMBER_OF_LOCI = %i\n",NUMBER_OF_LOCI);
	printf("Printing loci ...\n");
	//PrintLocuss(Loci, NUMBER_OF_LOCI,"Loci");
	char* reference_genome = argv[3];
	//"hg19_RefGene.txt";
	printf("Looking missing genes in %s ...\n",reference_genome);	
	int more_genes_found;
	more_genes_found = Use_Ref_Gene("hg19_RefGene_test.txt", Loci, NUMBER_OF_LOCI);
	printf("% i genes have been identified using %s file \n",more_genes_found, reference_genome);
//	PrintLocuss(Loci, NUMBER_OF_LOCI,"Loci");
	MinMax(Loci, NUMBER_OF_LOCI, SelectedGenes, NUMBER_OF_SELECTED_GENES);
	
 	//Compute_Mode_per_Genes(SelectedGenes, NUMBER_OF_SELECTED_GENES);
	printf("Printing genes ...\n");
//	PrintLocuss(SelectedGenes, NUMBER_OF_SELECTED_GENES,"Genes");
	//printf("Printing loci ...\n");
	//PrintLocuss(Loci, NUMBER_OF_LOCI,"Loci");
	
	double* weights;
	double LociMeanValue;
	double LociModeValue;
	double LociDispersion;
	double GenesMeanValue; 
	double GenesModeValue; 
 	double GenesDispersion;
 	
 	printf("\nThe statistical parameters\n");
 	w2 = w1 = 1.0;
	weights = init_weights(NUMBER_OF_LOCI,Loci,w1,w2);
	printf("w1 = %f,\t w2 = %f\n",w1,w2);

	LociMeanValue = MeanValue(Loci, NUMBER_OF_LOCI,weights, "Loci");
	LociDispersion = Dispersion(Loci,  NUMBER_OF_LOCI,weights,"Loci");
	LociModeValue = Mode(Loci,  NUMBER_OF_LOCI,"Loci");
 	GenesMeanValue = MeanValue(Loci,  NUMBER_OF_LOCI, weights,"Genes");
 	GenesModeValue = Mode(Loci,  NUMBER_OF_LOCI,"Genes");
 	GenesDispersion = Dispersion(Loci, NUMBER_OF_LOCI, weights,"Genes");
	
	printf("Loci MeanValue=%f\n", LociMeanValue);
	printf("Loci ModeValue=%f\n",	LociModeValue);
	printf("Loci Dispersion =%f\n",	LociDispersion);
	 	
	printf("Genes MeanValue=%f\n",	 GenesMeanValue);
	printf("Genes ModeValue=%f\n",	GenesModeValue);
 	printf("Genes Dispersion =%f\n",GenesDispersion);
//  	printf("Before sorting ...\n");
//  	PrintLocuss(Loci, NUMBER_OF_LOCI,"Loci");
	//sort(Loci, 4);

	printf("After sorting ...\n");
	//PrintLocuss(SelectedGenes, NUMBER_OF_LOCI,"Genes");
	

// 
//  	printf("\nThe statistical parameters\n");
//  	w2 = 1.0;
//  	w1 = 2.0;
//  	weights = init_weights(NUMBER_OF_LOCI,Loci,w1,w2);
// 	printf("w1 = %f,\t w2 = %f\n",w1,w2);
// 
//   	
//  	LociMeanValue = MeanValue(Loci, NUMBER_OF_LOCI,weights, "Loci");
// 	LociDispersion = Dispersion(Loci,  NUMBER_OF_LOCI,weights,"Loci");
//  	GenesMeanValue = MeanValue(Loci,  NUMBER_OF_LOCI, weights,"Genes");
//  	GenesDispersion = Dispersion(Loci, NUMBER_OF_LOCI, weights,"Genes");
// 	
// 	printf("Loci MeanValue=%f\n", LociMeanValue);
//  	printf("Loci Dispersion =%f\n",	 LociDispersion);
// 	
//  	printf("Genes MeanValue=%f\n",	 GenesMeanValue);
//  	printf("Genes Dispersion =%f\n", GenesDispersion);
// 
//  	printf("\nThe statistical parameters\n"); 	
//  	w1	= 3.0;	
//  	w2  = 1.0;
//  	weights = init_weights(NUMBER_OF_LOCI,Loci,w1,w2);
//  	printf("w1 = %f,\t w2 = %f\n",w1,w2);
//  	
//  	
//  	LociMeanValue = MeanValue(Loci, NUMBER_OF_LOCI,weights, "Loci");
// 	LociDispersion = Dispersion(Loci,  NUMBER_OF_LOCI,weights,"Loci");
// 	GenesMeanValue = MeanValue(Loci,  NUMBER_OF_LOCI, weights,"Genes");
//  	GenesDispersion = Dispersion(Loci, NUMBER_OF_LOCI, weights,"Genes");
// 
// 	printf("Loci MeanValue=%f\n", LociMeanValue);
//  	printf("Loci Dispersion=%f\n",	 LociDispersion);
//  	 	
//  	printf("Genes MeanValue=%f\n",	 GenesMeanValue);
//  	printf("Genes Dispersion =%f\n", GenesDispersion);
// 
//  	printf("\nThe statistical parameters \n"); 	
//  	w1	= 500000.0;	
//  	w2  = 1.0;
//  	weights = init_weights(NUMBER_OF_LOCI,Loci,w1,w2);
//  	printf("w1 = %f,\t w2 = %f\n",w1,w2);
//  	
//  	
//  	LociMeanValue = MeanValue(Loci, NUMBER_OF_LOCI,weights, "Loci");
// 	LociDispersion = Dispersion(Loci,  NUMBER_OF_LOCI,weights,"Loci");
// 	GenesMeanValue = MeanValue(Loci,  NUMBER_OF_LOCI, weights,"Genes");
//  	GenesDispersion = Dispersion(Loci, NUMBER_OF_LOCI, weights,"Genes");
// 
// 	printf("Loci MeanValue=%f\n", LociMeanValue);
//  	printf("Loci Dispersion=%f\n",	 LociDispersion);
//  	 	
//  	printf("Genes MeanValue=%f\n",	 GenesMeanValue);
//  	printf("Genes Dispersion =%f\n", GenesDispersion);
// 

	//free(Loci);
	//free(SelectedGenes);
	//free(weights);
	return 0;
}
	
