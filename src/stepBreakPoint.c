/****************************************************************************
 ****************************************************************************
 *                                                                          *
 *  Copyright (C) 2020  Genome Research Ltd.                                *
 *                                                                          *
 *  Author: Zemin Ning (zn1@sanger.ac.uk)                                   *
 *                                                                          *
 *  This file is part of covidPileup pipeline.                              *
 *                                                                          *
 *  covidPileup is a free software: you can redistribute it and/or modify it*
 *  under the terms of the GNU General Public License as published by the   *
 *  Free Software Foundation, either version 3 of the License, or (at your  *
 *  option) any later version.                                              *
 *                                                                          *
 *  This program is distributed in the hope that it will be useful, but     *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of              *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU        *
 *  General Public License for more details.                                *
 *                                                                          *
 *  You should have received a copy of the GNU General Public License along *
 *  with this program.  If not, see <http://www.gnu.org/licenses/>.         *
 *                                                                          *
 ****************************************************************************
 ****************************************************************************/
/****************************************************************************/

#include <values.h>
#include <stdio.h>
#include <netinet/in.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <ctype.h>
#include "fasta.h"

#define GT '>'
#define GT4 (((((GT<<8)+GT)<<8)+GT)<<8)+GT

#define ENDS_EXTRA 0
#define PADCHAR '-'
#define MAX_N_BRG 50000 
#define MAX_N_ROW 40000 
#define nfm 800000
#define nfm_sub 500000
#define Max_N_NameBase 60
#define Max_N_Pair 100
static long *h_dna;

/* SSAS default parameters   */
static int IMOD=0;
static int len_covid=200000000;
static int set_qual=15;
static int set_len=10;
static int a_len = 200;
static int long_flag =1;

int main(int argc, char **argv)
{
    FILE *fp,*namef,*namef2;
    long dataSize,totalBases;
    int i,j,k,m,n,nSeq,nRead,args,num_SNPs=0,num_GAPs=0,num_INDs;
    int proce_flag,n_Sbase,locX,locY,locD,locYY,flag_H,num_index;
    int nseq,num_align,*ctg_index,*ctg_hitst,*cig_base,*cig_code,*snp_offset,*snp_freq,*snp_refer,*ind_offset;
    fasta *seq,*seqp;
    long IntSeg,IntSegRC,IntBase,IntBaseRC;
    char *ptr,**ctgname,line[100000],cigarline[10000],*seqline,*read_base,*refe_base;
    char score1[60],score2[60],indel_base[500],indel_base2[500],refname[500];
    char **cmatrix(long nrl,long nrh,long ncl,long nch);
    fasta *segg;
    long Size_pdata,Size_q_pdata;
    unsigned int *pidata;
    int num_seqque;
    char line2[100000],line3[100000];
    char *pdata,*st,*ed;

    seq=NULL;
    fflush(stdout);
    if(system("ps aux | grep stepBreakPoint; date") == -1)
    {
//        printf("System command error:\n);
    }

    locX = 0;
    locY = 0;
    n_Sbase = 7;
    proce_flag = 1;
    if(argc < 2)
    {
      printf("Usage: %s <-length 40000> <input_reference_fasta> <alignment_file> <output_SNP_file>\n",argv[0]);
      exit(1);
    }

    nSeq=0;
    args=1;
    for(i=1;i<argc;i++)
    {
       if(!strcmp(argv[i],"-mod"))
       {
         sscanf(argv[++i],"%d",&IMOD); 
         args=args+2;
       }
       else if(!strcmp(argv[i],"-long"))
       {
         sscanf(argv[++i],"%d",&long_flag);
         args=args+2;
       }
    }

    /*
    if((fp=fopen(argv[args],"rb"))==NULL) printf("Cannot open file\n");
    fseek(fp, 0, SEEK_END);
    Size_q_pdata = ftell(fp) + 1;
    fclose(fp);
    if((pdata=(char*)calloc(Size_q_pdata,sizeof(char)))==NULL)
      printf("calloc pdata\n");
    num_seqque = extractFastq(argv[args],pdata,Size_q_pdata);
    if((segg=(fasta*)calloc((num_seqque+1),sizeof(fasta)))==NULL)
      printf("calloc segg\n");
    if((seq=decodeFastq(argv[args],&num_seqque,&totalBases,pdata,Size_q_pdata,segg))==NULL)
      printf("no query data found.\n");
    nseq=0;
    nSeq = num_seqque;
    printf("Number of shotgun reads  %d \n",nSeq);

    if(totalBases == 0)
      return EXIT_SUCCESS;
                                    */
/*  input read alignment info line   */
    if((namef = fopen(argv[args],"r")) == NULL)
    {
      printf("ERROR main:: reads group file: %s \n",argv[args+1]);
      exit(1);
    }

    nRead = 0;
    while(!feof(namef))
    {
      if(fgets(line,100000,namef) == NULL)
        printf("Data input file: %s\n",argv[args+1]);
      if(feof(namef)) break;
      nRead++;
    }
    fclose(namef); 

    printf("Number of genomes  %d %d\n",nSeq,nRead);

    ctgname=cmatrix(0,nRead,0,500);
    if((cig_base= (int *)calloc(10000,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: calloc - ctg_index\n");
      exit(1);
    }
    if((cig_code= (int *)calloc(10000,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: calloc - ctg_index\n");
      exit(1);
    }
    if((ind_offset= (int *)calloc(len_covid,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: calloc - ctg_index\n");
      exit(1);
    }
    if((snp_offset= (int *)calloc(len_covid,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: calloc - ctg_index\n");
      exit(1);
    }
    if((snp_freq= (int *)calloc(len_covid,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: calloc - ctg_index\n");
      exit(1);
    }
    if((snp_refer= (int *)calloc(totalBases,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: calloc - ctg_index\n");
      exit(1);
    }
    if((seqline= (char *)calloc(len_covid,sizeof(char))) == NULL)
    {
      printf("ERROR Memory_Allocate: calloc - seqline\n");
      exit(1);
    }
    if((read_base= (char *)calloc(len_covid,sizeof(char))) == NULL)
    {
      printf("ERROR Memory_Allocate: calloc - read_base\n");
      exit(1);
    }
    if((refe_base= (char *)calloc(len_covid,sizeof(char))) == NULL)
    {
      printf("ERROR Memory_Allocate: calloc - read_base\n");
      exit(1);
    }
    if((ctg_index= (int *)calloc(nRead,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: calloc - ctg_index\n");
      exit(1);
    }
    if((ctg_hitst= (int *)calloc(nRead,sizeof(int))) == NULL)
    {
      printf("ERROR Memory_Allocate: calloc - ctg_index\n");
      exit(1);
    }

    printf("Memory assigned!\n");

    dataSize=1000;
/*  process contigs one by one   */
    if((namef = fopen(argv[args],"r")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }
    if((namef2 = fopen(argv[args+1],"w")) == NULL)
    {
      printf("ERROR main:: reads group file \n");
      exit(1);
    }

/*  read the SNP output file         */
    num_align=0;
    while(!feof(namef))
    {
      int nPair=0,seq_len,hitst;
      char base[500],cbase[20];

      if(fgets(line,100000,namef) == NULL)
        printf("Data input file:\n");
      if(feof(namef)) break;
      strcpy(line2,line);
      strcpy(line3,line);
      proce_flag = 1;
      if(proce_flag)
      { 
        for(ptr=strtok(line," ");ptr!=NULL;ptr=strtok((char *)NULL," "),nPair++)
        {
        }
        i=0;
        for(ptr=strtok(line2," ");ptr!=NULL;ptr=strtok((char *)NULL," "),i++)
        {
           if(i==0)
           {
             memset(base,'\0',500);
             strcat(base,ptr);
             strcpy(ctgname[num_align],base);
           }
           else if(i==1)
           {
             memset(base,'\0',500);
             strcat(base,ptr);
             num_index = atoi(base);
           }
           else if(i==2)
           {
             memset(base,'\0',500);
             strcat(base,ptr);
             strcpy(refname,base);
           }
           else if(i==3)
           {
             memset(base,'\0',500);
             strcat(base,ptr);
             ctg_hitst[num_align] = atoi(base);
             hitst = atoi(base)-1;
           }
           else if(i==5)
           {
             memset(cigarline,'\0',10000);
             strcat(cigarline,ptr);
	   }
           else if(i==6)
           {
             int c_len,r_len,stopflag; 
             int num_hits,offset,s_len,h_len,r_offset,s_offset; 
             int gbase = 0;
             char *st,*ed;

/*
             memset(read_base,'\0',len_covid);
             memset(refe_base,'\0',len_covid);
             memset(seqline,'\0',len_covid);
             memset(snp_offset,'\0',4*len_covid);
             memset(ind_offset,'\0',4*len_covid);  */

             memset(seqline,'\0',10000);
             strcat(seqline,ptr);
             seq_len = strlen(seqline);
             st = cigarline;
             c_len = strlen(cigarline);
	     memset(cig_code,0,4*c_len);
	     memset(cig_base,0,4*c_len);
             c_len = strlen(cigarline)-1;
	     st = strchr(cigarline,'S');
	     flag_H = 0;
             for(k=0;k<c_len;k++)
             {
		if(cigarline[k] == 'M')
		  cig_code[k] = 1;
		else if(cigarline[k] == 'I')
		  cig_code[k] = 2;
		else if(cigarline[k] == 'D')
		  cig_code[k] = 3;
		else if(cigarline[k] == 'S')
		  cig_code[k] = 4;
		else if(cigarline[k] == 'H')
		  cig_code[k] = 5;
		else
		  cig_code[k] = 0;
             }
             num_hits = 0;
             offset = 0;
             r_offset = hitst;
             s_offset = 0;
	     s_len = 0;
	     h_len = 0;
	     for(k=0;k<c_len;k++)
             {
		stopflag = 0;
		j = k+1;
		while((j < c_len)&&(stopflag == 0))
                {
                   if(cig_code[j]==0)
                   {
                     j++;
                   }
                   else
                     stopflag=1;
                }
                memset(cbase,'\0',20);
                if(num_hits == 0)
                {
                  for(m=k;m<j;m++)
		    cbase[m-k] = cigarline[m];
                }
                else
                {
                  for(m=(k+1);m<j;m++)
		    cbase[m-k-1] = cigarline[m];
                } 
                r_len = atoi(cbase);
                if(cigarline[j] == 'S')
                {
	          s_len = r_len;
	          if(j<c_len)
                    s_offset = r_len;
//               printf("Cigar: %d %d %d %s %d %d %d %d\n",num_align,c_len,j,ctgname[num_align],hitst,offset,r_offset,s_offset);
                }
                else if(cigarline[j] == 'H')
                {
                  /*for(n=0;n<r_len;n++)
                  {
		     read_base[offset+n] = 'N';
		  }
                  offset = offset + r_len; */
		  flag_H = 1;
		  h_len = r_len;
//                  s_offset = r_len;
                }
                else if(cigarline[j] == 'M')
                {
/*                  for(n=0;n<r_len;n++)
                  {
                     refe_base[offset+n] = seq->data[r_offset+n]; 
                     read_base[offset+n] = seqline[s_offset+n];
                     snp_offset[offset+n] = r_offset+n+1; 
                  }   */
//                  printf("Base0: %d %d %d %c %c %c\n",num_align,n,offset,refe_base[n],read_base[n],seqline[n]);
                  offset = offset + r_len;
                  r_offset = r_offset + r_len;
                  s_offset = s_offset + r_len;
                }
                else if(cigarline[j] == 'I')
                {
/*                  for(n=0;n<r_len;n++)
                  {
                     refe_base[offset+n] = '-'; 
                     read_base[offset+n] = seqline[s_offset+n]; 
                     snp_offset[offset+n] = offset+1; 
                     ind_offset[offset+n] = r_offset+n+1; 
                  }   */
                  s_offset = s_offset+r_len;
                  offset = offset + r_len;
                }
                else if(cigarline[j] == 'D')
                {
/*                  for(n=0;n<r_len;n++)
                  {
                     refe_base[offset+n] = seq->data[r_offset+n]; 
                     read_base[offset+n] = '-'; 
                     snp_offset[offset+n] = r_offset+n+1; 
                     ind_offset[offset+n] = r_offset+n; 
                  }   */
                  offset = offset + r_len;
                  r_offset = r_offset+r_len;
                  
                }
                num_hits++;
                k = j-1;
             }
         
             num_SNPs = 0; 
             num_INDs = 0; 
             num_GAPs = 0;
	     if((strncmp(refname,"chr",3)) == 0)
	     {
	       int idt = 0;
	       if((strncmp(refname,"chrX",4)) == 0)
	       {
	         idt = 23;
	       }
	       else if((strncmp(refname,"chrY",4)) == 0)
	       {
	         idt = 24;
	       }
               else
	       {
	         ed = strrchr(refname,'r');
		 idt = atoi(ed+1);
	       }
//                 printf("Break1: %d %d %s %d %d %d %d %d\n",num_align,idt,ctgname[num_align],hitst,hitst,r_offset,s_offset,seq_len-1);
	       if((strlen(refname) < 10)&&(idt > 0))
	       {
		 if(long_flag == 1)
                   fprintf(namef2,"Break1: %d %d %s %d %d %d\n",num_align,idt,ctgname[num_align],hitst,hitst,r_offset);
		 else
	         {
  		     if((s_offset < (seq_len-1))||(seq_len < 120)||(s_len >= 10)||(h_len >= 10))
		     {
//		      printf("Break1: %d %d %s %d %d %d %d %d\n",num_align,idt,ctgname[num_align],hitst,hitst,r_offset,s_offset,seq_len-1);
                       fprintf(namef2,"Break1: %d %d %s %d %d %d\n",num_align,idt,ctgname[num_align],hitst,hitst,r_offset);
		     }
		 }
	       }
//               fprintf(namef2,"Break2: %d %d %s %d %d %d\n",num_align,idt,ctgname[num_align],r_offset,hitst,r_offset);
//               printf("Break1: %d %s %s %d %d %d\n",num_align,refname,ctgname[num_align],hitst,offset,seq_len-1);
//               printf("Break2: %d %s %s %d %d %d\n",num_align,refname,ctgname[num_align],r_offset,s_offset,seq_len-1);
//             fprintf(namef2,"Break1: %d %s %s %d %d %d\n",num_align,refname,ctgname[num_align],hitst,offset,seq_len-1);
//             fprintf(namef2,"Break2: %d %s %s %d %d %d\n",num_align,refname,ctgname[num_align],r_offset,s_offset,seq_len-1);
	     }
	     else if((strncmp(refname,"tar",3)) == 0)
	     {
	       int idt;
	       ed = strrchr(refname,'_');
               idt = atoi(ed+1)+1;
	       if(idt > 0)
	       {
	         if(long_flag == 1)
                   fprintf(namef2,"Break1: %d %d %s %d %d %d\n",num_align,idt,ctgname[num_align],hitst,hitst,r_offset);
	         else
	         {
  		   if((s_offset < (seq_len-1))||(seq_len < 120)||(s_len >= 10)||(h_len >= 10))
                     fprintf(namef2,"Break1: %d %d %s %d %d %d\n",num_align,idt,ctgname[num_align],hitst,hitst,r_offset);
	         }
	       }
//               fprintf(namef2,"Break2: %d %d %s %d %d %d\n",num_align,idt,ctgname[num_align],r_offset,hitst,r_offset);
	     }
           }
        }
      }
      num_align++;
    }
    fclose(namef);
    fclose(namef2);

/*
    if(seq){
        free(seq->name);
        free(seq);
        seq = NULL;
    }   

    nRead = 0;
    for(i=0;i<totalBases;i++)
    {
       if(snp_refer[i] > 0)
         nRead++;
    } 
    printf("SNP calls: %ld %d\n",totalBases,nRead); */

    return EXIT_SUCCESS;

}
/* end of the main */

#define SWAP(a,b) temp=(a);(a)=b;(b)=temp;

/*   Subroutine to sort an array arr[0,...,n-1] into ascending order while
     making the corresponding reaarangement of the array brr[0,...,n-1]
     by the use of Quicksort (Sedgwick, R. 1978, Communications o fthe ACM,
     vol. 21, pp. 847-857) also see Numerical Recipes in C                  */  

/* =============================== */
void ArraySort_Long(int n, long *arr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,NSTACK=50,istack[NSTACK];
     long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
             }
             arr[i+1]=a;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/* =============================== */
void ArraySort_Int(int n, int *arr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
             }
             arr[i+1]=a;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/* =============================== */
void ArraySort_Mix(int n, long *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/* =============================== */
void ArraySort_Int2(int n, int *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/*   function to sort an array into a decreasing order:  a>b>c>....    */  
/* =============================== */
void ArraySort2_Int2(int n, int *arr, int *brr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]>=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);

          if(arr[m]<arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
          }

          if(arr[m+1]<arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
          }

          if(arr[m]<arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          for(;;)
          {
             do i++; while (arr[i]>a);
             do j--; while (arr[j]<a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}

/* =============================== */
void ArraySort_Mix3(int n, long *arr, int *brr, int *crr)
/* =============================== */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,c,NSTACK=50,istack[NSTACK];
     long a,temp,MIN=7;

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             a=arr[j];
             b=brr[j];
             c=crr[j];
             for(i=j-1;i>=m;i--)
             {
                if(arr[i]<=a) break;
                arr[i+1]=arr[i];
                brr[i+1]=brr[i];
                crr[i+1]=crr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
             crr[i+1]=c;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          SWAP(arr[k],arr[m+1]);
          SWAP(brr[k],brr[m+1]);
          SWAP(crr[k],crr[m+1]);

          if(arr[m]>arr[ir])
          {
            SWAP(arr[m],arr[ir]);
            SWAP(brr[m],brr[ir]);
            SWAP(crr[m],crr[ir]);
          }

          if(arr[m+1]>arr[ir])
          {
            SWAP(arr[m+1],arr[ir]);
            SWAP(brr[m+1],brr[ir]);
            SWAP(crr[m+1],crr[ir]);
          }

          if(arr[m]>arr[m+1])
          {
            SWAP(arr[m],arr[m+1]);
            SWAP(brr[m],brr[m+1]);
            SWAP(crr[m],crr[m+1]);
          }

          i=m+1;
          j=ir;
          a=arr[m+1];
          b=brr[m+1];
          c=crr[m+1];
          for(;;)
          {
             do i++; while (arr[i]<a);
             do j--; while (arr[j]>a);
             if(j<i) break;
             SWAP(arr[i],arr[j]);
             SWAP(brr[i],brr[j]);
             SWAP(crr[i],crr[j]);
          }
          arr[m+1]=arr[j];
          arr[j]=a;
          brr[m+1]=brr[j];
          brr[j]=b;
          crr[m+1]=crr[j];
          crr[j]=c;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/*   to swap the string arrays           */
/* ============================================= */
void s_swap(char Pair_Name[][Max_N_NameBase], int i, int j)
/* ============================================= */
{
     char temp[Max_N_NameBase];

     strcpy(temp,Pair_Name[i]);
     strcpy(Pair_Name[i],Pair_Name[j]);
     strcpy(Pair_Name[j],temp);
}


/*   to sort the string array in order          */
/* ============================================= */
void ArraySort_String(int n, char Pair_Name[][Max_N_NameBase], int *brr)
/* ============================================= */
{
     int i,ir=n-1,j,k,m=0,jstack=0,b,NSTACK=50,istack[NSTACK];
     int temp,MIN=7;
     char p[Max_N_NameBase];

     for(;;)
     {
/*      Insertion sort when subarray is small enough    */
        if(ir-m<MIN)
        {
          for(j=m+1;j<=ir;j++)
          {
             strcpy(p,Pair_Name[j]);
             b=brr[j];
             for(i=j-1;i>=m;i--)
             {
                if(strcmp(Pair_Name[i],p)<=0) break;
                strcpy(Pair_Name[i+1],Pair_Name[i]);
                brr[i+1]=brr[i];
             }
             strcpy(Pair_Name[i+1],p);
             brr[i+1]=b;
          }
          if(!jstack) return;
          ir=istack[jstack--];
          m=istack[jstack--];
        }
        else
        {
          k=(m+ir)>>1;
          s_swap(Pair_Name,k,m+1);
          SWAP(brr[k],brr[m+1]);

          if(strcmp(Pair_Name[m],Pair_Name[ir])>0)
          {
            s_swap(Pair_Name,m,ir);
            SWAP(brr[m],brr[ir]);
          }

          if(strcmp(Pair_Name[m+1],Pair_Name[ir])>0)
          {
            s_swap(Pair_Name,m+1,ir);
            SWAP(brr[m+1],brr[ir]);
          }

          if(strcmp(Pair_Name[m],Pair_Name[m+1])>0)
          {
            s_swap(Pair_Name,m,m+1);
            SWAP(brr[m],brr[m+1]);
          }

          i=m+1;
          j=ir;
          strcpy(p,Pair_Name[m+1]);
          b=brr[m+1];
          for(;;)
          {
             do i++; while (strcmp(Pair_Name[i],p)<0);
             do j--; while (strcmp(Pair_Name[j],p)>0);
             if(j<i) break;
             s_swap(Pair_Name,i,j);
             SWAP(brr[i],brr[j]);
          }
          strcpy(Pair_Name[m+1],Pair_Name[j]);
          strcpy(Pair_Name[j],p);
          brr[m+1]=brr[j];
          brr[j]=b;
          jstack+=2;

/*        Push pointers to larger subarray on stack      */
/*        process smaller subarray immediately           */
          if(jstack>NSTACK)
          {
             printf("Stack error: NSTACK too small\n");
             exit(0);
          }
          if(ir-i+1>=j-m)
          {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
          }
          else
          {
            istack[jstack]=j-1;
            istack[jstack-1]=m;
            m=i;
          }
        }
     }
}


/* creat an int matrix with subscript ange m[nrl...nrh][ncl...nch]  */
int     **imatrix(long nrl,long nrh,long ncl,long nch)
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        int  **m;

        /* allocate pointers to rows        */
        if((m=(int **)calloc(nrow,sizeof(int*)))==NULL)
        {
           printf("error imatrix: calloc error No. 1 \n");
           return(NULL);
        }
        m+=0;
        m-=nrl;

        /* allocate rows and set pointers to them        */
        if((m[nrl]=(int *)calloc(nrow*ncol,sizeof(int)))==NULL)
        {
           printf("error imatrix: calloc error No. 2 \n");
           return(NULL);
        }
        m[nrl]+=0;
        m[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
           m[i]=m[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return m;
}

/* creat char matrix with subscript ange cm[nrl...nrh][ncl...nch]  */
char    **cmatrix(long nrl,long nrh,long ncl,long nch)
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        char **cm;

        /* allocate pointers to rows        */
        if((cm=(char **)calloc(nrow,sizeof(char*)))==NULL)
        {
           printf("error cmatrix: calloc error No. 1 \n");
           return(NULL);
        }
        cm+=0;
        cm-=nrl;

        /* allocate rows and set pointers to them        */
        if((cm[nrl]=(char *)calloc(nrow*ncol,sizeof(char)))==NULL)
        {
           printf("error cmatrix: calloc error No. 2 \n");
           return(NULL);
        }
        cm[nrl]+=0;
        cm[nrl]-=nrl;

        for(i=nrl+1;i<=nrh;i++)
           cm[i]=cm[i-1]+ncol;
        /* return pointer to array of pointers to rows   */
        return cm;
}

