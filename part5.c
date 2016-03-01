//
//  parser.c
//  
//
//  Created by Spiros Chatzikotoulas 
//  and Dimitris Kalfountzos 
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "hash_2.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "csparse.h"


#define SIZE 20
#define HASH_SIZE 65536


int team2_element=0,nzASparse=0;//plithos omadas 2
double *Apinakas;
double *bpinakas;
double *xpinakas;
double itol=1e-3;
struct entry_s *groundNode;
gsl_matrix_view Apin;
gsl_vector *x;
gsl_permutation * p;
int flagLU=1, flagCholesky=0, dcSweep=0, dcPlot=0, position=0,flag_geiwsis=0, flagCG=0, flagBiCG=0, flagItol=0, flagSparseLU=0, flagSparseChol=0, flagSparseCG=0, flagSparseBiCG=0;
char stepStr[SIZE], beginStr[SIZE], endStr[SIZE], elStr[SIZE], namStr[SIZE];
struct entry_s *node_sweep;
cs *Apinakas_sparse_compress=NULL;




//List for I, V, R, L, C, D elements of a circuit

struct listElements{
    
    char element,name[SIZE]; //name of the element, id, positive and negative nodes of the element
    float value, length, width; //value of the element
    struct entry_s  *pNode, *nNode,*mosD, *mosG, *mosS, *mosB, *bjtC, *bjtB, *bjtE; 
    struct listElements *next;
    struct listElements *previous;
    
};

hashtable_t *hashtable;




int foundComment(char *str); //Function to detect comments after ---> '*'
struct listElements * listInitialize(struct listElements *rootNode); //Initialize List of elements I, V, R, C, L and D
struct listElements * addElementsList(struct listElements *rootNode, char element, char *name, struct entry_s *pNode, struct entry_s *nNode, struct entry_s *mosD, struct entry_s *mosG, struct entry_s *mosS, struct entry_s *mosB, struct entry_s *bjtC, struct entry_s *bjtB, struct entry_s *bjtE, float value, float length, float width);
struct listElements * trimInput(char *line2, struct listElements *rootNode);//Function to detect whitespace characters and call the addElementsList to add the detected values
void printList(struct listElements *rootNode);//Prints the elements of the circuit
int analyse_MNA(struct listElements *rootNode);
void print_Pinakes(int counter_team2,int dimension);
void LU(int dimension);
void Cholesky(int dimension);
void solve(int dimension,int flagCholesky);
int detectName(char type, struct listElements *rootNode, char *nameCheck);
int DC(struct listElements *rootNode,int flagCholesky);
void sweepDC(int dimension,struct listElements *rootNode);
void solve_diagwnio(double *dianisma1, double *dianisma2, double *apotelesma, int length);
void BiCG(int dimension);
void CG(int dimension);
void create_Inverse_Pinaka_Preconditioner(double *inverse_pinakas_preconditioner, int length);
int analyse_MNA_sparse(struct listElements *rootNode);
void LU_sparse(int dimension);
void Cholesky_sparse(int dimension);




int main(int argc, char* argv[])
{
    
    FILE *f;
    char line;
    char *line2=NULL;
    int lengthC=0, i=0, comment=0, dimension=0;
    
      if(argc!=2){
	
	printf("wrong number of arguments\n");
	return(1);
	
	}
    
    
    f=fopen(argv[1],"r+");
    
    
    if(f==NULL){
    
        printf("ERROR: Couldn't open the FILE\n"); //Error FILE
        return(0);
    }
    
    hashtable = ht_create( HASH_SIZE );
    
    groundNode=(struct entry_s *)malloc(sizeof(struct entry_s));
    groundNode->key = "0\0";
    groundNode->id_node = 0;
    
    last_node_num=0;
    ground = 0;
    
    struct listElements *rootNode=NULL;
    rootNode = listInitialize(rootNode);

    
    while(1){
        
        i=0;
        
        while(1){
            
            line=fgetc(f);
            
            if(line==EOF){

                break;
            }
            else if(line=='\n'){

                break;
            }
            else{
                line2=(char*)realloc(line2,(lengthC+1)*sizeof(char));
                line2[i]=line;
                i++;
                lengthC++;
                }
            
            }
        line2=(char*)realloc(line2,(lengthC+1)*sizeof(char));
        line2[i]='\0';
        
            comment=foundComment(line2);
            if(comment==1){
               //printf("there is a comment\n");
            }
            else if(comment==2){
                //printf("Empty Row\n");
            }
            else{
                rootNode=trimInput(line2, rootNode);
            }
            
            
            lengthC=0;
            
            if(feof(f)){
               //printf("End of file reached\n");
                break;
            }
            
            
    }
        //free(line2);
        fclose(f);
        //detectGround(rootNode);
	
	if (ground){
	  //printList(rootNode);
	}
	  else{
	    
	    printf("there is no ground\n");
	  }

   
	dimension=DC(rootNode,flagCholesky);
	
	if(dcSweep==1){
		sweepDC(dimension,rootNode);
	}
	

    
    return(0);

    
}
    

    

//Detect comments in the file

    
int foundComment(char *line2){
    
        char comment='*';
        
        while(1){
            
            while(isspace(*line2)){
                line2+=1;
                continue;
            }
            if(*line2=='\0'){
                
                    return(2); //found the end of the string
            }
            else if(*line2==comment){
                    
                    return(1); //found comment
                }
            else{
                return(0);
                line2+=1;
                }
               
        }
}
               
//initialize the list
               
struct listElements * listInitialize(struct listElements *rootNode){
                   
                   
    rootNode=(struct listElements *)malloc(sizeof(struct listElements));   //allocate memory
    rootNode->next = rootNode;
    rootNode->previous = rootNode;
    
       
    return(rootNode);
    
}



//add elements

struct listElements * addElementsList(struct listElements *rootNode, char element, char *name, struct entry_s *pNode, struct entry_s *nNode, struct entry_s *mosD, struct entry_s *mosG, struct entry_s *mosS, struct entry_s *mosB, struct entry_s *bjtC, struct entry_s *bjtB, struct entry_s *bjtE, float value, float length, float width) {
        
    struct listElements *elements;
    
    elements=(struct listElements *)malloc(sizeof(struct listElements));
    
    elements->element = element;
    strcpy(elements->name , name);
    elements->pNode = pNode;
    elements->nNode = nNode;
    elements->value=value;
    elements->mosD = mosD;
    elements->mosG = mosG;
    elements->mosS = mosS;
    elements->mosB = mosB;
    elements->bjtC = bjtC;
    elements->bjtB = bjtB;
    elements->bjtE = bjtE;
    elements->length=length;
    elements->width=width;
    
    elements->next=rootNode;
    elements->previous=rootNode->previous;
    rootNode->previous->next=elements;
    rootNode->previous=elements;
    

    return(rootNode);
}


            

//Trim whitespaces from the input detect numbers and strings and add them to the list
    
struct listElements * trimInput(char *line2, struct listElements *rootNode){
        
    int i=0,j=0;   //counters
    char elementStr[SIZE], tmp_str[SIZE];
    char *optionsStr;
    char nameStr[SIZE];  //numbers to string format
    char valueStr[SIZE];
    int pNode=0, nNode=0;
    //char stepStr[SIZE];
    //char beginStr[SIZE];
    //char endStr[SIZE];
    /*char mosGStr[SIZE];
    char mosSStr[SIZE];
    char mosBStr[SIZE];
    char bjtCStr[SIZE];
    char bjtBStr[SIZE];
    char bjtEStr[SIZE];*/
    struct entry_s *pNodeStr, *nNodeStr, *mosDStr, *mosGStr, *mosSStr, *mosBStr, *bjtCStr, *bjtBStr, *bjtEStr;
    char lengthStr[SIZE];
    char widthStr[SIZE];
    //flagCholesky=0;
    
    

    
    //detect the input as it is expected
    
        while(isspace(line2[i])){
            i++;
            continue;
        }
    if((line2[i]=='I')||(line2[i]=='V')||(line2[i]=='R')||(line2[i]=='L')||(line2[i]=='C')||(line2[i]=='i')||(line2[i]=='v')||(line2[i]=='r')||(line2[i]=='l')||(line2[i]=='c')){
	    
	    if((line2[i]=='V')||(line2[i]=='L')||(line2[i]=='v')||(line2[i]=='l')){
	      team2_element++;
	      //printf("\n%d\n",team2_element);
	    }
      
            elementStr[j]=line2[i];
            i++;
            j++;

        elementStr[j]='\0'; //found the end of the element
        j=0;

    
        while(isspace(line2[i])){ //ignore multiple spaces or tabs
            i++;
        }
    
    
        while(isalnum(line2[i])){
            nameStr[j]=line2[i];
            //printf("%s\n",idStr);
            j++;
            i++;
        }
        
        nameStr[j]='\0'; // found the end of the id	    
        j=0;
    
        while(isspace(line2[i])){ //ignore multiple spaces or tabs
            i++;
        }
    
        while(isdigit(line2[i])||isalpha(line2[i])){
            tmp_str[j]=line2[i];
            j++;
            i++;
        }
        tmp_str[j]='\0'; //found the positive node
        pNode=atoi(tmp_str);
        pNodeStr = ht_set( hashtable, tmp_str, groundNode);
        j=0;

        while(isspace(line2[i])){ //ignore multiple spaces or tabs
            i++;
        }
    
        while(isdigit(line2[i])||isalpha(line2[i])){
            tmp_str[j]=line2[i];
            j++;
            i++;
        }
        tmp_str[j]='\0'; //found the negative node
        nNode=atoi(tmp_str);
        nNodeStr = ht_set( hashtable, tmp_str, groundNode);
        
        j=0;

        while(isspace(line2[i])){ //ignore multiple spaces or tabs
            i++;
        }
    
        while((line2[i]!='\0')&&!(isspace(line2[i]))){
        
        valueStr[j]=line2[i];
        //printf("%s\n", valueStr);
        j++;
        i++;
        
        }
    
        valueStr[j]='\0';  //found the value of the element
        j=0;
    
    
        
       
        
        
    //detection of input line ended, add it to the list
	
	if(!strcasecmp("R",elementStr))
	{
	   if(nNode >0 && pNode>0)	
	   {
		 nzASparse+=4; 
	   }
	   else if(nNode >0) {
		 nzASparse += 1;   
	   }
	   else if(pNode>0) {
		nzASparse += 1;   
	   }
	}
	
	if(!strcasecmp("V",elementStr))
	{
	   if(nNode >0 && pNode>0)	
	   {
		 nzASparse+=4; 
	   }
	   else if(nNode >0) {
		 nzASparse += 2;   
	   }
	   else if(pNode>0) {
		nzASparse += 2;   
	   }
	}
	
	if(!strcasecmp("L",elementStr))
	{
	   if(nNode >0 && pNode>0)	
	   {
		 nzASparse+=4; 
	   }
	   else if(nNode >0) {
		 nzASparse += 2;   
	   }
	   else if(pNode>0) {
		nzASparse += 2;   
	   }
	}
	
	printf("%d\n", nzASparse);
         
        
       // printf("the I element %d\n", nzASparse);

        
        rootNode=addElementsList(rootNode, elementStr[0], nameStr, pNodeStr, nNodeStr, NULL, NULL, NULL, NULL, NULL, NULL, NULL, atof(valueStr), atof(lengthStr),atof(widthStr));
    
    }
    else if(line2[i]=='D'){
        
        while(isalpha(line2[i])){
            elementStr[j]=line2[i];
            i++;
            j++;
        }
        elementStr[j]='\0'; //found the end of the element
        j=0;
        //printf("%s\n",elementStr);
        
        
        while(isspace(line2[i])){ //ignore multiple spaces or tabs
            i++;
        }
        
        
        while(isdigit(line2[i])||isalpha(line2[i])){
            nameStr[j]=line2[i];
            //printf("%s\n",idStr);
            j++;
            i++;
        }
        
        nameStr[j]='\0'; // found the end of the id
        j=0;
        
        while(isspace(line2[i])){ //ignore multiple spaces or tabs
            i++;
        }
        
        while(isdigit(line2[i])||isalpha(line2[i])){
            tmp_str[j]=line2[i];
            j++;
            i++;
        }
        tmp_str[j]='\0'; //found the positive node
        pNodeStr = ht_set( hashtable, tmp_str, groundNode);
        j=0;

        while(isspace(line2[i])){ //ignore multiple spaces or tabs
            i++;
        }
    
        while(isdigit(line2[i])||isalpha(line2[i])){
            tmp_str[j]=line2[i];
            j++;
            i++;
        }
        tmp_str[j]='\0'; //found the negative node
        nNodeStr = ht_set( hashtable, tmp_str, groundNode);
        j=0;
        
        
        
        rootNode=addElementsList(rootNode, elementStr[0], nameStr, pNodeStr, nNodeStr, NULL, NULL, NULL, NULL, NULL, NULL, NULL, atof(valueStr), atof(lengthStr),atof(widthStr));
        
    }
    else if(line2[i]=='M'){
        
        while(isalpha(line2[i])){
            elementStr[j]=line2[i];
            i++;
            j++;
        }
        elementStr[j]='\0'; //found the end of the element
        j=0;
        
        
        while(isspace(line2[i])){ //ignore multiple spaces or tabs
            i++;
        }
        
        while(isdigit(line2[i])||isalpha(line2[i])){
            nameStr[j]=line2[i];
            j++;
            i++;
        }
        
        nameStr[j]='\0'; // found the end of the id
        //printf("%s\n",idStr);
        j=0;
        
        while(isspace(line2[i])){
            i++;
        }
        
        while(isdigit(line2[i])||isalpha(line2[i])){
            tmp_str[j]=line2[i];
            j++;
            i++;
        }
        tmp_str[j]='\0'; // found the end of the mosD
        mosDStr = ht_set( hashtable, tmp_str, groundNode);
        //printf("%s\n",mosDStr);
        j=0;
        
        while(isspace(line2[i])){
            i++;
        }
        
        while(isdigit(line2[i])||isalpha(line2[i])){
            tmp_str[j]=line2[i];
            j++;
            i++;
        }
        
        tmp_str[j]='\0'; // found the end of the mosG
        mosGStr = ht_set( hashtable, tmp_str, groundNode);
        j=0;
        
        while(isspace(line2[i])){
            i++;
        }
        
        while(isdigit(line2[i])||isalpha(line2[i])){
            tmp_str[j]=line2[i];
            j++;
            i++;
        }
        
        tmp_str[j]='\0'; // found the end of the mosS
        mosSStr = ht_set( hashtable, tmp_str, groundNode);
        j=0;
        
        while(isspace(line2[i])){
            i++;
        }
        
        while(isdigit(line2[i])||isalpha(line2[i])){
            tmp_str[j]=line2[i];
            
            j++;
            i++;
        }
        
        tmp_str[j]='\0'; // found the end of the mosB
        mosBStr = ht_set( hashtable, tmp_str, groundNode);
        j=0;
        
        
        while(isspace(line2[i])){
            i++;
        }
        
        while(isdigit(line2[i])||isalpha(line2[i])){
            lengthStr[j]=line2[i];
            j++;
            i++;
        }
        
        lengthStr[j]='\0'; // found the end of the width
        j=0;
        
        while(isspace(line2[i])){
            i++;
        }
        
        while((line2[i]!='\0')&&!(isspace(line2[i]))){
            widthStr[j]=line2[i];
            j++;
            i++;
        }
        
        widthStr[j]='\0';
        j=0;
        
        
       rootNode=addElementsList(rootNode, elementStr[0], nameStr, NULL, NULL, mosDStr, mosGStr, mosSStr, mosBStr, NULL, NULL, NULL, atof(valueStr), atof(lengthStr),atof(widthStr));
        
    }
    else if(line2[i]=='Q'){
      
        
        while(isalpha(line2[i])){
            elementStr[j]=line2[i];
            i++;
            j++;
        }
        elementStr[j]='\0'; //found the end of the element
        j=0;
        
        
        while(isspace(line2[i])){
            i++;
        }
        
        
        while(isdigit(line2[i])||isalpha(line2[i])){
            nameStr[j]=line2[i];
            j++;
            i++;
        }
        
        nameStr[j]='\0'; // found the end of the id
        j=0;
        
        while(isspace(line2[i])){
            i++;
        }
        
        while(isdigit(line2[i])||isalpha(line2[i])){
            tmp_str[j]=line2[i];
            j++;
            i++;
        }
        tmp_str[j]='\0'; //found the Collector
        bjtCStr = ht_set( hashtable, tmp_str, groundNode);
        j=0;
        
        while(isspace(line2[i])){
            i++;
        }
        
        while(isdigit(line2[i])||isalpha(line2[i])){
            tmp_str[j]=line2[i];
            j++;
            i++;
        }
        tmp_str[j]='\0'; //found the Base
        bjtBStr = ht_set( hashtable, tmp_str, groundNode);
        j=0;
        
        while(isspace(line2[i])){
            i++;
        }
        
        while((line2[i]!='\0')&&!(isspace(line2[i]))){
            
            tmp_str[j]=line2[i];
            j++;
            i++;
            
        }
        
        tmp_str[j]='\0';  //found the value of the ekpompos
        bjtEStr = ht_set( hashtable, tmp_str, groundNode);
        j=0;
        
        
        rootNode=addElementsList(rootNode, elementStr[0], nameStr, NULL, NULL, NULL, NULL, NULL, NULL, bjtCStr, bjtBStr, bjtEStr, atof(valueStr), atof(lengthStr),atof(widthStr));


        
    }
    
    //flag for SPD, DC, PLOT put extra recognition for tokens different V and I
    else if(line2[i]=='.'){
        
        optionsStr=strtok(line2,". \t\n");
        
        if(!strcmp(optionsStr, "OPTIONS")){
            
            optionsStr = strtok(NULL," ");
            
            if(optionsStr!=NULL){
                
                
                if(optionsStr[0]=='S'){
                    if(optionsStr[2]=='D'){
                        optionsStr[3]='\0';
                    }
                    else if(optionsStr[2]=='A'){
                        optionsStr[6]='\0';
                    }
                }
                else if(optionsStr[0]=='I'){
                    optionsStr[4]='\0';
                }
                if(!strcmp(optionsStr, "SPD")){
                    
                    flagCholesky=1;  //[flag to enabe cholesky option]
                    flagLU=0;
                    flagBiCG=0;
                    flagCG=0;
                    
                    printf("flagcholesky %d\n",flagCholesky);
                    optionsStr=strtok(NULL," ");
                    if(!strcmp(optionsStr=="ITER")){
                        
                        flagSparseLU=0;
                        flagSparseChol=0;
                        flagSparseCG=0;
                        flagSparseBiCG=0;
                        flagCG=1;
                        flagCholesky=0;
                        flagLU=0;
                        flagBiCG=0;
                        printf("I'm here for CG %d\n",flagCG);

                        
                    }
                    else if((optionsStr!=NULL)&&(optionsStr[2]=='A')){
                        optionsStr[6]='\0';
                        if(!strcmp(optionsStr, "SPARSE")){
                            
                            flagSparseLU=0;
                            flagSparseChol=1;
                            flagSparseCG=0;
                            flagSparseBiCG=0;
                            flagCG=0;
                            flagCholesky=0;
                            flagLU=0;
                            flagBiCG=0;
                            
                            printf("I'm here for SPD SPARSE %d\n",flagSparseChol);
                            
                        }
                    }
                    
                }
                else if(!strcmp(optionsStr,"ITER")){
                    
                    flagBiCG=1;
                    flagCholesky=0;
                    flagLU=0;
                    flagCG=0;
                    //printf("%s\n",optionsStr);
                    printf("the flag BI-CG is %d\n",flagBiCG);
                    optionsStr=strtok(NULL," ");
                    if((optionsStr!=NULL)&&(optionsStr[2]=='A')){
                        optionsStr[6]='\0';
                        if(!strcmp(optionsStr, "SPARSE")){//we have to make changes for the flags!!! global SPARSE flag only
                            
                            flagSparseLU=0;
                            flagSparseChol=0;
                            flagSparseCG=0;
                            flagSparseBiCG=1;
                            flagCG=0;
                            flagCholesky=0;
                            flagLU=0;
                            flagBiCG=0;
                            printf("I'm here for ITER SPARSE %d\n",flagSparseBiCG);
                            
                        }
                    }
                    
                }
                else if(!strcmp(optionsStr,"SPARSE")){
                    
                    
                    flagSparseLU=1;
                    flagSparseChol=0;
                    flagSparseCG=0;
                    flagSparseBiCG=0;
                    flagCG=0;
                    flagCholesky=0;
                    flagLU=0;
                    flagBiCG=0;
                    printf("I'm here for SPRASE  %d\n",flagSparseLU);
                    
                }
                else if(!strcmp(optionsStr,"ITOL")){
                    
                    optionsStr=strtok(NULL," ");
                    
                    if(optionsStr!=NULL){
                        
                        while(optionsStr[i]!='\0'){
                            
                            valueStr[j]=optionsStr[i];
                            j++;
                            i++;
                        }
                        valueStr[j]='\0';
                        itol=atof(valueStr);
                        flagItol=1;
                        printf("the itol flag is %d\n",flagItol);
                        printf("the itol is %e\n",itol);
                        j=0;
                        i=0;
                    }
                    else{
                        
                        flagItol=1;
                        //printf("the itol flag is %d\n",flagItol);
                        //printf("the itol is %e\n",itol);
                    }
                    
                }
                
                
                
                if((optionsStr!=NULL)&&(optionsStr[0]=='I')){
                    
                    if(!strcmp(optionsStr, "ITER")){
                        
                        flagCG=1;
                        flagCholesky=0;
                        flagLU=0;
                        flagBiCG=0;
                        printf("the flag CG is %d\n",flagCG);
                        //printf("%d\n",flagCG);
                        //printf("%d\n",flagCholesky);
                        optionsStr=strtok(NULL," ");
                        optionsStr[6]='\0';
                        if((optionsStr!=NULL)&&(optionsStr[2]=='A')){
                            
                            if(!strcmp(optionsStr, "SPARSE")){
                                
                                flagSparseLU=0;
                                flagSparseChol=0;
                                flagSparseCG=1;
                                flagSparseBiCG=0;
                                flagCG=0;
                                flagCholesky=0;
                                flagLU=0;
                                flagBiCG=0;
                                
                                printf("I'm here for SPD ITER SPARSE %d\n",flagSparseCG);
                                
                            }
                        }
                        
                        
                        
                    }
                }
                
            }
        }
        
        else if(!strcmp(optionsStr,"DC")){
         
            
                optionsStr = strtok(NULL," ");
                elStr[j]=optionsStr[i];
                j++;
                elStr[j]='\0';
                j=0;
            
                if((elStr[i]=='V')||(elStr[i]=='v')||(elStr[i]=='I')||(elStr[i]=='i')){
            
			char type = elStr[i];
                        i++;
                        while(optionsStr[i]!='\0'){
                    
                                namStr[j]=optionsStr[i];
                                j++;
                                i++;
                        }
                    
                        namStr[j]='\0';
                        //printf("%s\n",nameStr);
                        j=0;
                        optionsStr = strtok(NULL," ");
                        i=0;
                    
      
                    
                        while(optionsStr[i]!='\0'){
                    
                                beginStr[j]=optionsStr[i];
                                j++;
                                i++;
                    
                        }
    
                        beginStr[j]='\0';
                        //printf("%s\n",beginStr);
                        j=0;
                
                        optionsStr = strtok(NULL," ");
                        i=0;
                
                
                        while(optionsStr[i]!='\0'){
                    
                            endStr[j]=optionsStr[i];
                            j++;
                            i++;
                    
                        }
                        endStr[j]='\0';
                        //printf("%s\n",endStr);
                        j=0;
                
                
                        optionsStr = strtok(NULL," ");
                        i=0;
                
                        while(optionsStr[i]!='\0'){
                    
                        stepStr[j]=optionsStr[i];
                        j++;
                        i++;
                    
                        }
                
                        stepStr[j]='\0';
                       // printf("%s\n",stepStr);

                        j=0;
                        j=detectName(type, rootNode,namStr);
                    
                    
                    if(j==0){
                        
                        printf("There isn't an element with that id in the Netlist\n");
                        exit(-1);
                        }
			dcSweep++;
                }

            
        }
    
        else if(!strcmp(optionsStr,"PLOT")){
            
                optionsStr = strtok(NULL," ");
                i=1;
                if(optionsStr[i]=='('){
                
                    i++;
                    while(isdigit(optionsStr[i])||isalpha(optionsStr[i])){
                        nameStr[j]=optionsStr[i];
                        i++;
                        j++;

                    }
                
                }
                nameStr[j]='\0';
                //printf("%s\n",nameStr);
                j=0;
                //j=detectName(rootNode,nameStr);

                if(j==0){
                
                    printf("there isn't an element with that id in the Netlist\n");
                
                }
		dcPlot++;
            }
        }
    else{
        
        printf("unknown element\n");
        
        }
    
    
    return(rootNode);
}

// Name Detection
int detectName(char type, struct listElements *rootNode, char *nameCheck){
    
    struct listElements *elementP;
    int flag=0;
    
    for(elementP=rootNode->next;elementP->next!=rootNode->next;elementP=elementP->next){
        
        if((elementP->element== type)&&(!(strcmp(elementP->name,nameCheck)))){
            flag=flag+1;
        }

    }

    return(flag);
    
}


//Print the elements of te circuit

void printList(struct listElements *rootNode){
    
    struct listElements *elementP;
    
    for(elementP=rootNode->next;elementP->next!=rootNode->next;elementP=elementP->next){
        
        if((elementP->element=='I')||(elementP->element=='V')||(elementP->element=='R')||(elementP->element=='L')||(elementP->element=='C')||(elementP->element=='D')||(elementP->element=='i')||(elementP->element=='v')||(elementP->element=='r')||(elementP->element=='l')||(elementP->element=='c')||(elementP->element=='d'))
        
            printf("%c%s - %s - %s - %e \n",elementP->element, elementP->name, (elementP->pNode)->key, (elementP->nNode)->key, elementP->value);
        else if(elementP->element=='M'){
            
            printf("%c%s - %s - %s - %s - %s - %e - %e\n", elementP->element, elementP->name, (elementP->mosD)->key, (elementP->mosG)->key, (elementP->mosS)->key, (elementP->mosB)->key, elementP->length, elementP->width);
        }
        else if(elementP->element=='Q'){
            
            printf("%c%s - %s - %s - %s \n",elementP->element, elementP->name, (elementP->bjtC)->key, (elementP->bjtB)->key, (elementP->bjtE)->key);

            
        }

    }
}


int analyse_MNA(struct listElements *rootNode){
    
    int counter_team2=0;
    int dimension=0;
    struct listElements *tmp_deiktis;
    
    //printf("\nPlihtos komvwn xwris tin pigi einai %d\n",last_node_num);
    //printf("\nplithos stoixeiwn omadas 2 einai %d\n",team2_element);
    
    dimension=(last_node_num)+team2_element;
    printf("\n to dimension einai %d\n\n",dimension);
    
    Apinakas = (double*)calloc(dimension*dimension,sizeof(double));
    bpinakas = (double*)calloc(dimension,sizeof(double));
    xpinakas= (double *) calloc(dimension, sizeof(double));
    
    if(Apinakas==NULL || bpinakas==NULL){
      
	printf("NO ALLOCATE MEMORY\n");
	exit(-1);
    }
    
    
    tmp_deiktis=rootNode;
    
    
    tmp_deiktis=tmp_deiktis->next;
    while(tmp_deiktis!=rootNode){
      
	switch(tmp_deiktis->element){
	  case'r':
	  case'R':
		  if(tmp_deiktis->pNode->id_node!=0 && tmp_deiktis->nNode->id_node!=0){
		  
		    
		  Apinakas[(tmp_deiktis->pNode->id_node-1)*dimension + (tmp_deiktis->pNode->id_node-1)] += 1/(tmp_deiktis->value);
		  Apinakas[(tmp_deiktis->nNode->id_node-1)*dimension + (tmp_deiktis->nNode->id_node-1)] += 1/(tmp_deiktis->value);
		  Apinakas[(tmp_deiktis->pNode->id_node-1)*dimension + (tmp_deiktis->nNode->id_node-1)] -= 1/(tmp_deiktis->value);
		  Apinakas[(tmp_deiktis->nNode->id_node-1)*dimension + (tmp_deiktis->pNode->id_node-1)] -= 1/(tmp_deiktis->value);
		 // printf("+ - stoiveio %d, %d\n", (tmp_deiktis->pNode->id_node-1), (tmp_deiktis->pNode->id_node-1));
		  }
		  
		  if(tmp_deiktis->pNode->id_node!=0 && tmp_deiktis->nNode->id_node==0){
		    
			Apinakas[(tmp_deiktis->pNode->id_node-1)*dimension + (tmp_deiktis->pNode->id_node-1)] += 1/(tmp_deiktis->value);
			//printf("+ stoiveio %d +%lf\n", (tmp_deiktis->pNode->id_node-1), 1/(tmp_deiktis->value));
		    
		  }
		  
		  if(tmp_deiktis->pNode->id_node==0 && tmp_deiktis->nNode->id_node!=0){
		    
			Apinakas[(tmp_deiktis->nNode->id_node-1)*dimension + (tmp_deiktis->nNode->id_node-1)] += 1/(tmp_deiktis->value);
			//printf("- stoiveio %d\n", (tmp_deiktis->nNode->id_node-1)*dimension + (tmp_deiktis->nNode->id_node-1));
		    
		  }
		  //printf("\nkala paw RRR\n");
		  
		  
		  break;
		  
	  case'v':	  
	  case'V':
		  
		  
		 if(tmp_deiktis->pNode->id_node!=0 && tmp_deiktis->nNode->id_node!=0){
			Apinakas[(counter_team2 + last_node_num)*dimension + (tmp_deiktis->pNode->id_node-1)] += 1;
			Apinakas[(counter_team2 + last_node_num)*dimension + (tmp_deiktis->nNode->id_node-1)] -= 1;
			Apinakas[(counter_team2 + last_node_num) + dimension*(tmp_deiktis->pNode->id_node-1)] += 1;
			Apinakas[(counter_team2 + last_node_num) + dimension*(tmp_deiktis->nNode->id_node-1)] -= 1;
		   
		   
		 }
		  
		 if(tmp_deiktis->pNode->id_node!=0 && tmp_deiktis->nNode->id_node==0){
			Apinakas[(counter_team2 + last_node_num)*dimension + (tmp_deiktis->pNode->id_node-1)] += 1;
			Apinakas[(counter_team2 + last_node_num) + dimension*(tmp_deiktis->pNode->id_node-1)] += 1;
		   
		   
		} 
		
		if(tmp_deiktis->pNode->id_node==0 && tmp_deiktis->nNode->id_node!=0){
			Apinakas[(counter_team2 + last_node_num)*dimension + (tmp_deiktis->nNode->id_node-1)] -= 1;
			Apinakas[(counter_team2 + last_node_num) + dimension*(tmp_deiktis->nNode->id_node-1)] -= 1;
		  
		  
		}
		
		bpinakas[(counter_team2 + last_node_num)] += (tmp_deiktis->value);
		
		if (strcmp(namStr,tmp_deiktis->name)==0){
			
			position=counter_team2;
		}
		
		counter_team2++;
		
		  
		  //printf("\nkala paw VVV\n");
		  break;
		  
	  case'i':	  
	  case'I':
		  if(tmp_deiktis->pNode->id_node!=0 && tmp_deiktis->nNode->id_node!=0){
		  
			bpinakas[(tmp_deiktis->pNode->id_node-1)] -= (tmp_deiktis->value);
			bpinakas[(tmp_deiktis->nNode->id_node-1)] += (tmp_deiktis->value);
		    
		  }
		  
		  if(tmp_deiktis->pNode->id_node!=0 && tmp_deiktis->nNode->id_node==0){
		   
		    
			bpinakas[(tmp_deiktis->pNode->id_node-1)] -= (tmp_deiktis->value);
		  }
		  
		  if(tmp_deiktis->pNode->id_node==0 && tmp_deiktis->nNode->id_node!=0){
			
			bpinakas[(tmp_deiktis->nNode->id_node-1)] += (tmp_deiktis->value);
		  }
		  //printf("\nkala paw III\n");
		  
	
		  break;
		  
		  
	  case'l':	  
	  case'L':  
		  
		  
		if(tmp_deiktis->pNode->id_node!=0 && tmp_deiktis->nNode->id_node!=0){
			Apinakas[(counter_team2 + last_node_num)*dimension + (tmp_deiktis->pNode->id_node-1)] += 1;
			Apinakas[(counter_team2 + last_node_num)*dimension + (tmp_deiktis->nNode->id_node-1)] -= 1;
			Apinakas[(counter_team2 + last_node_num) + dimension*(tmp_deiktis->pNode->id_node-1)] += 1;
			Apinakas[(counter_team2 + last_node_num) + dimension*(tmp_deiktis->nNode->id_node-1)] -= 1;
		   
		   
		 }
		  
		 if(tmp_deiktis->pNode->id_node!=0 && tmp_deiktis->nNode->id_node==0){
			Apinakas[(counter_team2 + last_node_num)*dimension + (tmp_deiktis->pNode->id_node-1)] += 1;
			Apinakas[(counter_team2 + last_node_num) + dimension*(tmp_deiktis->pNode->id_node-1)] += 1;
		   
		   
		} 
		
		if(tmp_deiktis->pNode->id_node==0 && tmp_deiktis->nNode->id_node!=0){
			Apinakas[(counter_team2 + last_node_num)*dimension + (tmp_deiktis->nNode->id_node-1)] -= 1;
			Apinakas[(counter_team2 + last_node_num) + dimension*(tmp_deiktis->nNode->id_node-1)] -= 1;
		  
		  
		}
		
		bpinakas[(counter_team2 + last_node_num)] += (tmp_deiktis->value);
		  
		  counter_team2++;
		  //printf("\nkala paw LLL\n");
		  break;
	}
      
      
      tmp_deiktis=tmp_deiktis->next;
      
    }
    
	//print_Pinakes(counter_team2,dimension);
	
	return(dimension);
}

int analyse_MNA_sparse(struct listElements *rootNode){
    
    int counter_team2=0,i=0;
    int dimension=0;
    struct listElements *tmp_deiktis;
    cs *Apinakas_sparse=NULL;
    //printf("\nPlihtos komvwn xwris tin pigi einai %d\n",last_node_num);
    //printf("\nplithos stoixeiwn omadas 2 einai %d\n",team2_element);
    
    dimension=(last_node_num)+team2_element;
    printf("\n to dimension einai %d\n\n",dimension);
    
    //nzASparse = 39;
    
    printf("%d\n", nzASparse);
    
    Apinakas_sparse = cs_spalloc(dimension,dimension,nzASparse,1,1);
    bpinakas = (double*)calloc(dimension,sizeof(double));
    xpinakas= (double *)calloc(dimension, sizeof(double));
    
    if(Apinakas_sparse==NULL){
      
	printf("NO ALLOCATE MATRIX\n");
	exit(-1);
    }
 
    tmp_deiktis=rootNode;
    
    
    tmp_deiktis=tmp_deiktis->next;
    while(tmp_deiktis!=rootNode){
      
	switch(tmp_deiktis->element){

	  case'r':
	  case'R':
		  if(tmp_deiktis->pNode->id_node!=0 && tmp_deiktis->nNode->id_node!=0){
		  
			cs_entry(Apinakas_sparse, tmp_deiktis->pNode->id_node-1,tmp_deiktis->nNode->id_node-1, -(1/tmp_deiktis->value));
			//Apinakas_sparse-> i[i]=tmp_deiktis->pNode->id_node-1;
			//Apinakas_sparse-> p[i]=tmp_deiktis->nNode->id_node-1;
			//Apinakas_sparse-> x[i]=-(1/tmp_deiktis->value);
			//i++;
			cs_entry(Apinakas_sparse, tmp_deiktis->pNode->id_node-1,tmp_deiktis->pNode->id_node-1, (1/tmp_deiktis->value));
			cs_entry(Apinakas_sparse, tmp_deiktis->nNode->id_node-1,tmp_deiktis->pNode->id_node-1, -(1/tmp_deiktis->value));
			cs_entry(Apinakas_sparse, tmp_deiktis->nNode->id_node-1,tmp_deiktis->nNode->id_node-1, (1/tmp_deiktis->value));
			//Apinakas_sparse-> i[i]=tmp_deiktis->nNode->id_node-1;
			//Apinakas_sparse-> p[i]=tmp_deiktis->pNode->id_node-1;
			//Apinakas_sparse-> x[i]=-(1/tmp_deiktis->value);
			//i++;
			
		  }
		  
		  if(tmp_deiktis->pNode->id_node!=0 && tmp_deiktis->nNode->id_node==0){
		    
			//Apinakas_sparse-> i[i]=tmp_deiktis->pNode->id_node-1;
			//Apinakas_sparse-> p[i]=tmp_deiktis->pNode->id_node-1;
			//Apinakas_sparse-> x[i]=-(1/tmp_deiktis->value);
			//i++;
			cs_entry(Apinakas_sparse, tmp_deiktis->pNode->id_node-1,tmp_deiktis->pNode->id_node-1, (1/tmp_deiktis->value));
		    
		  }
		  
		  if(tmp_deiktis->pNode->id_node==0 && tmp_deiktis->nNode->id_node!=0){
		    
			//Apinakas_sparse-> i[i]=tmp_deiktis->nNode->id_node-1;
			//Apinakas_sparse-> p[i]=tmp_deiktis->nNode->id_node-1;
			//Apinakas_sparse-> x[i]=-(1/tmp_deiktis->value);
			//i++;
			cs_entry(Apinakas_sparse, tmp_deiktis->nNode->id_node-1,tmp_deiktis->nNode->id_node-1, (1/tmp_deiktis->value));
		    
		  }
		  //printf("\nkala paw RRR\n");
		  
		  
		  break;
	  
	  case'v':	  
	  case'V':
		  
		  
		 if(tmp_deiktis->pNode->id_node!=0 && tmp_deiktis->nNode->id_node!=0){
			
			 cs_entry(Apinakas_sparse, tmp_deiktis->pNode->id_node-1, counter_team2+last_node_num, 1);
			 //Apinakas_sparse->i[i]=tmp_deiktis->pNode->id_node-1;
			 //Apinakas_sparse->p[i]=(counter_team2+last_node_num);
			 //Apinakas_sparse->x[i]=1;
			 //i++;
			 cs_entry(Apinakas_sparse, counter_team2+last_node_num, tmp_deiktis->pNode->id_node-1, 1);
			 cs_entry(Apinakas_sparse, tmp_deiktis->nNode->id_node-1, counter_team2+last_node_num, -1);
			 cs_entry(Apinakas_sparse, counter_team2+last_node_num, tmp_deiktis->nNode->id_node-1, -1);
			 
			 //Apinakas_sparse->p[i]=tmp_deiktis->pNode->id_node-1;
			 //Apinakas_sparse->i[i]=(counter_team2+last_node_num);
			 //Apinakas_sparse->x[i]=1;
			 //i++; 
			 
			 //Apinakas_sparse->i[i]=tmp_deiktis->nNode->id_node-1;
			 //Apinakas_sparse->p[i]=(counter_team2+last_node_num);
			 //Apinakas_sparse->x[i]=-1;
			 //i++;
			 
			 //Apinakas_sparse->p[i]=tmp_deiktis->nNode->id_node-1;
			 //Apinakas_sparse->i[i]=(counter_team2+last_node_num);
			 //Apinakas_sparse->x[i]=-1;
			 //i++;
		   
		 }
		  
		 if(tmp_deiktis->pNode->id_node!=0 && tmp_deiktis->nNode->id_node==0){

			 //Apinakas_sparse->i[i]=tmp_deiktis->pNode->id_node-1;
			 //Apinakas_sparse->p[i]=(counter_team2+last_node_num);
			 //Apinakas_sparse->x[i]=1;
			 //i++;
			 
			 
			 //Apinakas_sparse->p[i]=tmp_deiktis->pNode->id_node-1;
			 //Apinakas_sparse->i[i]=(counter_team2+last_node_num);
			 //Apinakas_sparse->x[i]=1;
			 //i++;
 			 cs_entry(Apinakas_sparse, tmp_deiktis->pNode->id_node-1, counter_team2+last_node_num, 1);
			 //Apinakas_sparse->i[i]=tmp_deiktis->pNode->id_node-1;
			 //Apinakas_sparse->p[i]=(counter_team2+last_node_num);
			 //Apinakas_sparse->x[i]=1;
			 //i++;
			 cs_entry(Apinakas_sparse, counter_team2+last_node_num, tmp_deiktis->pNode->id_node-1, 1);

		   
		} 
		
		if(tmp_deiktis->pNode->id_node==0 && tmp_deiktis->nNode->id_node!=0){

			//Apinakas_sparse->i[i]=tmp_deiktis->nNode->id_node-1;
			//Apinakas_sparse->p[i]=(counter_team2+last_node_num);
			//Apinakas_sparse->x[i]=-1;
			//i++;
			
			//Apinakas_sparse->p[i]=tmp_deiktis->nNode->id_node-1;
			//Apinakas_sparse->i[i]=(counter_team2+last_node_num);
			//Apinakas_sparse->x[i]=-1;
			//i++;
 			cs_entry(Apinakas_sparse, tmp_deiktis->nNode->id_node-1, counter_team2+last_node_num, -1);
			//Apinakas_sparse->i[i]=tmp_deiktis->pNode->id_node-1;
			//Apinakas_sparse->p[i]=(counter_team2+last_node_num);
			//Apinakas_sparse->x[i]=1;
			//i++;
			cs_entry(Apinakas_sparse, counter_team2+last_node_num, tmp_deiktis->nNode->id_node-1, -1);
			
		  
		}
		
		bpinakas[(counter_team2 + last_node_num)] += (tmp_deiktis->value);
		
		if (strcmp(namStr,tmp_deiktis->name)==0){
			
			position=counter_team2;
		}
		
		counter_team2++;
		
		  
		  //printf("\nkala paw VVV\n");
		  break;
		  
	  case'i':	  
	  case'I':
	  if(tmp_deiktis->pNode->id_node!=0 && tmp_deiktis->nNode->id_node!=0){
		  
			bpinakas[(tmp_deiktis->pNode->id_node-1)] -= (tmp_deiktis->value);
			bpinakas[(tmp_deiktis->nNode->id_node-1)] += (tmp_deiktis->value);
		    
		  }
		  
		  if(tmp_deiktis->pNode->id_node!=0 && tmp_deiktis->nNode->id_node==0){
		   
		    
			bpinakas[(tmp_deiktis->pNode->id_node-1)] -= (tmp_deiktis->value);
		  }
		  
		  if(tmp_deiktis->pNode->id_node==0 && tmp_deiktis->nNode->id_node!=0){
			
			bpinakas[(tmp_deiktis->nNode->id_node-1)] += (tmp_deiktis->value);
		  }
		  //printf("\nkala paw III\n")
		  
	
		  break;
		  
		  
	  case'l':	  
	  case'L':  
		  
		  
		 if(tmp_deiktis->pNode->id_node!=0 && tmp_deiktis->nNode->id_node!=0){
			
			 cs_entry(Apinakas_sparse, tmp_deiktis->pNode->id_node-1, counter_team2+last_node_num, 1);
			 //Apinakas_sparse->i[i]=tmp_deiktis->pNode->id_node-1;
			 //Apinakas_sparse->p[i]=(counter_team2+last_node_num);
			 //Apinakas_sparse->x[i]=1;
			 //i++;
			 cs_entry(Apinakas_sparse, counter_team2+last_node_num, tmp_deiktis->pNode->id_node-1, 1);
			 cs_entry(Apinakas_sparse, tmp_deiktis->nNode->id_node-1, counter_team2+last_node_num, -1);
			 cs_entry(Apinakas_sparse, counter_team2+last_node_num, tmp_deiktis->nNode->id_node-1, -1);
			 
			 //Apinakas_sparse->p[i]=tmp_deiktis->pNode->id_node-1;
			 //Apinakas_sparse->i[i]=(counter_team2+last_node_num);
			 //Apinakas_sparse->x[i]=1;
			 //i++; 
			 
			 //Apinakas_sparse->i[i]=tmp_deiktis->nNode->id_node-1;
			 //Apinakas_sparse->p[i]=(counter_team2+last_node_num);
			 //Apinakas_sparse->x[i]=-1;
			 //i++;
			 
			 //Apinakas_sparse->p[i]=tmp_deiktis->nNode->id_node-1;
			 //Apinakas_sparse->i[i]=(counter_team2+last_node_num);
			 //Apinakas_sparse->x[i]=-1;
			 //i++;
		   
		 }
		  
		 if(tmp_deiktis->pNode->id_node!=0 && tmp_deiktis->nNode->id_node==0){

			 //Apinakas_sparse->i[i]=tmp_deiktis->pNode->id_node-1;
			 //Apinakas_sparse->p[i]=(counter_team2+last_node_num);
			 //Apinakas_sparse->x[i]=1;
			 //i++;
			 
			 
			 //Apinakas_sparse->p[i]=tmp_deiktis->pNode->id_node-1;
			 //Apinakas_sparse->i[i]=(counter_team2+last_node_num);
			 //Apinakas_sparse->x[i]=1;
			 //i++;
 			 cs_entry(Apinakas_sparse, tmp_deiktis->pNode->id_node-1, counter_team2+last_node_num, 1);
			 //Apinakas_sparse->i[i]=tmp_deiktis->pNode->id_node-1;
			 //Apinakas_sparse->p[i]=(counter_team2+last_node_num);
			 //Apinakas_sparse->x[i]=1;
			 //i++;
			 cs_entry(Apinakas_sparse, counter_team2+last_node_num, tmp_deiktis->pNode->id_node-1, 1);

		   
		} 
		
		if(tmp_deiktis->pNode->id_node==0 && tmp_deiktis->nNode->id_node!=0){

			//Apinakas_sparse->i[i]=tmp_deiktis->nNode->id_node-1;
			//Apinakas_sparse->p[i]=(counter_team2+last_node_num);
			//Apinakas_sparse->x[i]=-1;
			//i++;
			
			//Apinakas_sparse->p[i]=tmp_deiktis->nNode->id_node-1;
			//Apinakas_sparse->i[i]=(counter_team2+last_node_num);
			//Apinakas_sparse->x[i]=-1;
			//i++;
 			cs_entry(Apinakas_sparse, tmp_deiktis->nNode->id_node-1, counter_team2+last_node_num, -1);
			//Apinakas_sparse->i[i]=tmp_deiktis->pNode->id_node-1;
			//Apinakas_sparse->p[i]=(counter_team2+last_node_num);
			//Apinakas_sparse->x[i]=1;
			//i++;
			cs_entry(Apinakas_sparse, counter_team2+last_node_num, tmp_deiktis->nNode->id_node-1, -1);
			
		  
		}
		
		bpinakas[(counter_team2 + last_node_num)] += 0;//(tmp_deiktis->value);
		  
		  counter_team2++;
		  //printf("\nkala paw LLL\n");
		  break;

	}
	
	//printf("%d\n", i);

      
      tmp_deiktis=tmp_deiktis->next;
      
    }

    
 
    Apinakas_sparse->nz=nzASparse;
    
    
    for(i=0;i<=nzASparse;i++){
    printf("the %d -st element from sparse matrix: %d\n",i, Apinakas_sparse->p[i]);
    }
    Apinakas_sparse_compress=cs_compress(Apinakas_sparse);
   
    cs_spfree(Apinakas_sparse);  
    
    cs_dupl(Apinakas_sparse_compress);
    
       
	//print_Pinakes(counter_team2,dimension);
	
	return(dimension);
}

void print_Pinakes(int counter_team2,int dimension){
  
    int i,j;
  //test prints
 printf("\nMNA\n\n");

    printf("\npinakas A:\n\n");
    for(counter_team2 = 0; counter_team2< dimension; counter_team2++) {
        for(i = 0; i < dimension; i++) {
   printf("%12lf ", Apinakas[counter_team2 * dimension + i]);
   if(i == last_node_num-1) {
    printf("   |");
   }
  }
  printf("\n");
 
  if(counter_team2 == last_node_num-1) {
   for(j = 0; j < dimension * 14; j++) {
    printf("_");
   }
   printf("\n\n");
  }
 }
 
 printf("\n\n\npinakas b:\n\n");
 for(counter_team2 = 0; counter_team2 < dimension; counter_team2++) {
  printf("%12lf ", bpinakas[counter_team2]);
  printf("\n");
  if(counter_team2 == last_node_num-1) {
   printf("_____________\n\n");
  }
    }
    
    printf("\n\n");
}


void LU(int dimension){
	
	
	//dimension=analyse_MNA(rootNode);

	Apin = gsl_matrix_view_array (Apinakas, dimension, dimension);
	
	int s;

	p = gsl_permutation_alloc (dimension);

	gsl_linalg_LU_decomp (&Apin.matrix, p, &s);



	//gsl_permutation_free (p);
	//gsl_vector_free (x);
		
  
}


void Cholesky(int dimension){
	
	
	//dimension=analyse_MNA(rootNode);
	
	Apin = gsl_matrix_view_array (Apinakas, dimension, dimension);
	
	gsl_linalg_cholesky_decomp (&Apin.matrix);
	
	
}


void solve(int dimension,int flagCholesky){
	
	//int i;
	
	gsl_vector_view bpin
	= gsl_vector_view_array(bpinakas, dimension);
	
	//printf("ela eimai kala\n");
	
	if(flagCholesky==1){
		//printf("pro chol\n");
		gsl_linalg_cholesky_solve (&Apin.matrix , &bpin.vector, x);
		//printf("meta chol\n");
	}
	if(flagLU==1){
		//printf("pro lu\n");
		gsl_linalg_LU_solve (&Apin.matrix, p, &bpin.vector, x);
		//printf("meta lu\n");
	}
	//printf("meta apo lu kai chhol\n");
	
	/*for(i=0;i<dimension;i++){
		
		printf("%12lf\n",gsl_vector_get(x,i));
	}*/
	printf ("x :  \n");
	gsl_vector_fprintf (stdout, x, "%g");
	
}


int DC(struct listElements *rootNode,int flagCholesky){
	
	int dimension;
	
	
	
	if(flagCholesky==1){
		
		dimension=analyse_MNA(rootNode);
		x=gsl_vector_alloc (dimension);
		
		Cholesky(dimension);
		solve(dimension,flagCholesky);
		printf("Kalesa Chol\n");
	}
	if(flagLU==1){
		dimension=analyse_MNA(rootNode);
		x=gsl_vector_alloc (dimension);
		
		
		LU(dimension);
		solve(dimension,flagCholesky);
		printf("Kalesa LU\n");
	}
	
	
	if(flagBiCG==1){
		
		dimension=analyse_MNA(rootNode);
		x=gsl_vector_alloc (dimension);
		
		
		printf("Mpika BiCG\n");
		BiCG (dimension);
		printf("Kalesa BiCG\n");
	}
	
	if(flagCG==1){
		
		dimension=analyse_MNA(rootNode);
		x=gsl_vector_alloc (dimension);
		
		CG(dimension);
		printf("Kalesa CG\n");
	}
	
	if(flagSparseLU==1){
		printf("prin kalesw MNA gia sparse \n");
		dimension=analyse_MNA_sparse(rootNode);
		printf("prin kalesw LU sparse \n");
		LU_sparse(dimension);
		printf("kalesa LU SPARSE\n");
	}
	
	if(flagSparseChol==1){
		printf("prin kalesw MNA gia sparse \n");
		dimension=analyse_MNA_sparse(rootNode);
		printf("prin kalesw Chol sparse \n");
		Cholesky_sparse(dimension);
		printf("kalesa Chol SPARSE\n");
	}
	
	
	
	return(dimension);
	
}


void sweepDC(int dimension,struct listElements *rootNode){
	
	
	double i, begin=atof(beginStr), end=atof(endStr), step=atof(stepStr);
	struct listElements *tmp_deiktis;
	
	tmp_deiktis=rootNode;
    
	tmp_deiktis=tmp_deiktis->next;
	
	while(tmp_deiktis!=rootNode){
		
		if((tmp_deiktis->element==elStr[0]) && strcmp(tmp_deiktis->name,namStr)==0 ){
			
			break;
		}
		
		tmp_deiktis=tmp_deiktis->next;
		
		
	}
	
	
	
	for(i=begin;i<=end;i+=step){
		
		printf("\n\nstep = %lf\n\n",i);
		
		if(elStr[0]=='V'||elStr[0]=='v'){
			bpinakas[last_node_num+position]=i;
			
		}
		
			
		else if(elStr[0]=='I'||elStr[0]=='i'){
			
			if(tmp_deiktis->pNode!=groundNode){
				
				bpinakas[tmp_deiktis->pNode->id_node-1]=-i;
			}
			
			if(tmp_deiktis->nNode!=groundNode){
				
				bpinakas[tmp_deiktis->nNode->id_node-1]=i;
			}
			
		}	
		else{
			
			printf("oute V oute I\n");
			exit(-1);
		}
		
		if(flagCholesky==1 || flagLU==1){
		solve(dimension,flagCholesky);
		}
		
		if(flagBiCG==1){
			
			BiCG(dimension);
		}
		
		if(flagCG==1){
			
			CG(dimension);
		}
	}
	
	
	
}




void create_Inverse_Pinaka_Preconditioner(double *inverse_pinakas_preconditioner, int length) {
  int i;
  double tmp;
  
  for(i = 0; i < length; i++) {
  	tmp = Apinakas[i * length + i]; //pairnei ta diagwnia stoixeia tou pinaka A
  	inverse_pinakas_preconditioner[i] = 1 / (tmp == 0 ? 1 : tmp); // kanei tin antisrofi
  }      
}


void solve_diagwnio(double *dianisma1, double *dianisma2, double *apotelesma, int length) {
  int i;
  
    for(i = 0; i < length; i++) {
        apotelesma[i] = dianisma1[i] * dianisma2[i];
    }
    
}

void CG(int dimension) {
  int k;
  int iteration = 0;
  
  double tmp, r, r1, alfa, vita, tmp3;
  double *dianisma_r = NULL, *dianisma_z = NULL, *inverse_pinakas_preconditioner = NULL, *dianisma_p = NULL, *dianisma_q = NULL;
  
    
 
    dianisma_r = (double *) calloc(dimension, sizeof(double));
    dianisma_z = (double *) malloc(dimension * sizeof(double));
    inverse_pinakas_preconditioner = (double *) malloc(dimension * sizeof(double));
    dianisma_p = (double *) malloc(dimension * sizeof(double));
    dianisma_q = (double *) malloc(dimension * sizeof(double));
    
  
    if(dianisma_r == NULL || dianisma_z == NULL || inverse_pinakas_preconditioner == NULL || dianisma_p == NULL || dianisma_q == NULL) {
        printf("Den exoun arxikopoihthei ta dianismata.\n");
        exit(-1);
    }
    
    create_Inverse_Pinaka_Preconditioner(inverse_pinakas_preconditioner, dimension);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, dimension, dimension, 1.0, Apinakas, dimension, xpinakas, 1, 0.0, dianisma_r, 1);  
    cblas_dscal(dimension, -1.0, dianisma_r, 1);
    cblas_daxpy(dimension, 1.0, bpinakas, 1, dianisma_r, 1);
    tmp = cblas_ddot(dimension, bpinakas, 1, bpinakas, 1);
    tmp = itol * itol *  (tmp == 0 ? 1 : tmp);
    
    while((tmp3 = cblas_ddot(dimension, dianisma_r, 1, dianisma_r, 1)) > tmp && iteration < dimension) {
        iteration++; 
        
        solve_diagwnio(inverse_pinakas_preconditioner, dianisma_r, dianisma_z, dimension); // z = inv(M) * r to z einai antistrofos tou M epi r
      
        r = cblas_ddot(dimension, dianisma_r, 1, dianisma_z, 1);
        
        if(iteration == 1) {
            memcpy(dianisma_p, dianisma_z, dimension * sizeof(double)); // kanoume to dianusma p iso me z
        }
        else {
             vita = r / r1;
            
            cblas_dscal(dimension, vita, dianisma_p, 1); // p=b*p
            cblas_daxpy(dimension, 1.0, dianisma_z, 1, dianisma_p, 1);
        }
        
        r1 = r;

        cblas_dgemv(CblasRowMajor, CblasNoTrans, dimension, dimension, 1.0, Apinakas, dimension, dianisma_p, 1, 0.0, dianisma_q, 1);
        alfa = r / cblas_ddot(dimension, dianisma_p, 1, dianisma_q, 1);
        cblas_daxpy(dimension, alfa, dianisma_p, 1, xpinakas, 1);
        cblas_daxpy(dimension, -alfa, dianisma_q, 1, dianisma_r, 1);

    }
   
   
    printf("\ndianisma x:\n\n");
    for(k = 0; k < dimension; k++) {
        printf("%12lf\n", xpinakas[k]);
    }
    
    printf("\n\n");

    
    free(dianisma_r);
    free(dianisma_z);
    free(inverse_pinakas_preconditioner);
    free(dianisma_p);
    free(dianisma_q);
  
}

void BiCG(int dimension) {
  int iteration = 0;
  
  double tmp, r, r1, alfa, vita, tmp3, omega;
  double EPS = 1e-14;
  double *dianisma_r = NULL, *dianisma_r2 = NULL, *dianisma_z = NULL, *dianisma_z2 = NULL, *inverse_pinakas_preconditioner = NULL, *dianisma_p = NULL, *dianisma_p2 = NULL, *dianisma_q = NULL, *dianisma_q2 = NULL;
  int k;
  
 
    dianisma_r = (double *) calloc(dimension, sizeof(double));
    dianisma_r2 = (double *) malloc(dimension * sizeof(double));
    dianisma_z = (double *) malloc(dimension * sizeof(double));
    dianisma_z2 = (double *) malloc(dimension * sizeof(double));
    inverse_pinakas_preconditioner = (double *) malloc(dimension * sizeof(double));
    dianisma_p = (double *) malloc(dimension * sizeof(double));
    dianisma_p2 = (double *) malloc(dimension * sizeof(double));
    dianisma_q = (double *) malloc(dimension * sizeof(double));
    dianisma_q2 = (double *) malloc(dimension * sizeof(double));
    
      
    if(dianisma_r == NULL || dianisma_r2 == NULL || dianisma_z == NULL || dianisma_z2 == NULL || inverse_pinakas_preconditioner == NULL || dianisma_p == NULL || dianisma_p2 == NULL || dianisma_q == NULL || dianisma_q2 == NULL) {
        printf("Den exoun arxikopoihthei ta dianismata.\n");
        exit(-1);
    }
    
  
    
    create_Inverse_Pinaka_Preconditioner(inverse_pinakas_preconditioner, dimension);  // antistrofos M = 1/diagwnios(G)
    
    cblas_dgemv(CblasRowMajor, CblasNoTrans, dimension, dimension, 1.0, Apinakas, dimension, xpinakas, 1, 0.0, dianisma_r, 1);
    cblas_dscal(dimension, -1.0, dianisma_r, 1);
    cblas_daxpy(dimension, 1.0, bpinakas, 1, dianisma_r, 1); 
    memcpy(dianisma_r2, dianisma_r, dimension * sizeof(double));
    tmp = cblas_ddot(dimension, bpinakas, 1, bpinakas, 1);
    tmp = itol * itol *  (tmp == 0 ? 1 : tmp);
    

    while((tmp3 = cblas_ddot(dimension, dianisma_r, 1, dianisma_r, 1)) > tmp && iteration < dimension) {
        iteration++; 
        solve_diagwnio(inverse_pinakas_preconditioner, dianisma_r, dianisma_z, dimension);  //epilush diagwniou sustymatos
        solve_diagwnio(inverse_pinakas_preconditioner, dianisma_r2, dianisma_z2, dimension);
        r = cblas_ddot(dimension, dianisma_z, 1, dianisma_r2, 1);
        
        if(fabs(r) < EPS) { 
            printf("Entopistike diairesi me 0,opote lathos.\n");
            exit(-1);
        }
        
        if(iteration == 1) {
            memcpy(dianisma_p, dianisma_z, dimension * sizeof(double));
            memcpy(dianisma_p2, dianisma_z2, dimension * sizeof(double));
        }
        else {
            vita = r / r1;
            cblas_dscal(dimension, vita, dianisma_p, 1);
            cblas_daxpy(dimension, 1.0, dianisma_z, 1, dianisma_p, 1);
            cblas_dscal(dimension, vita, dianisma_p2, 1);
            cblas_daxpy(dimension, 1.0, dianisma_z2, 1, dianisma_p2, 1);
        }
        
        r1 = r;
        cblas_dgemv(CblasRowMajor, CblasNoTrans, dimension, dimension, 1.0, Apinakas, dimension, dianisma_p, 1, 0.0, dianisma_q, 1);
        cblas_dgemv(CblasRowMajor, CblasTrans, dimension, dimension, 1.0, Apinakas, dimension, dianisma_p2, 1, 0.0, dianisma_q2, 1);
     
        omega = cblas_ddot(dimension, dianisma_p2, 1, dianisma_q, 1);
        
        if(fabs(omega) < EPS) { // ean to wmega einai 0 exit
           printf("Entopistike diairesi me 0,opote lathos.\n");
            exit(-1);
        }
        
        alfa = r / omega;
        cblas_daxpy(dimension, alfa, dianisma_p, 1, xpinakas, 1);
        cblas_daxpy(dimension, -alfa, dianisma_q, 1, dianisma_r, 1);
        cblas_daxpy(dimension, -alfa, dianisma_q2, 1, dianisma_r2, 1);

    }
    
    
    printf("dianisma x:\n");
    for(k = 0; k < dimension; k++) {
        printf("%12lf\n", xpinakas[k]);
    }
    
    printf("\n\n");
    
    
    free(dianisma_r);
    free(dianisma_r2);
    free(dianisma_z);
    free(dianisma_z2);
    free(inverse_pinakas_preconditioner);
    free(dianisma_p);
    free(dianisma_p2);
    free(dianisma_q);
    free(dianisma_q2);
  
}


void LU_sparse(int dimension){
	
	
	int i;
	css  *S=NULL;
	csn  *N=NULL;
	double *temp_dianisma=NULL;//proswrino dianisma
	
	
	temp_dianisma=(double *)calloc(dimension, sizeof(double));

	S=cs_sqr(2,Apinakas_sparse_compress,0);

	N=cs_lu(Apinakas_sparse_compress,S,1);
	
	cs_spfree(Apinakas_sparse_compress);
	

	
	memcpy(xpinakas, bpinakas, dimension * sizeof(double));
	
	cs_ipvec(N -> pinv, xpinakas, temp_dianisma, dimension);
        cs_lsolve(N -> L, temp_dianisma);
        cs_usolve(N -> U, temp_dianisma);
        cs_ipvec(S -> q, temp_dianisma, xpinakas, dimension);

	for (i = 0; i < dimension; i++)
	{
		printf("x[%d]: %f\n", i, xpinakas[i]);
	}
	
	free(temp_dianisma);
		
}

void Cholesky_sparse(int dimension){
	
	int i;
	css  *S=NULL;
	csn  *N=NULL;
	double *temp_dianisma=NULL;//proswrino dianisma
	
	
	temp_dianisma=(double *)calloc(dimension, sizeof(double));
	

	S=cs_schol(1,Apinakas_sparse_compress);
	N=cs_chol(Apinakas_sparse_compress,S);
	//cs_spfree(Apinakas_sparse_compress);
	if(N==NULL || S==NULL)
	{
		printf("Mallon trexeis arxeio pou dn einai gia cholesky\n");
		printf("Telos\n");
					fflush(stdout);
					exit(0);
	}
	memcpy(xpinakas, bpinakas, dimension * sizeof(double));
	cs_ipvec(S->pinv,xpinakas,temp_dianisma,dimension);
	cs_lsolve(N->L,temp_dianisma);
	cs_ltsolve(N->L,temp_dianisma);
	cs_pvec(S->pinv,temp_dianisma,xpinakas,dimension);
	
        for (i = 0; i < dimension; i++)
	{
		printf("x[%d]: %f\n", i, xpinakas[i]);
	}
	
	free(temp_dianisma);
}



void BiCG_sparse(int dimension){
	
	int i;
	
	
	
	
	
	
	
	
	
	
	
	for (i = 0; i < dimension; i++){
		printf("x[%d]: %f\n", i, xpinakas[i]);
	}
}
