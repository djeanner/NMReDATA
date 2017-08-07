#include <stdio.h>
#include <string.h>

FILE *file;
FILE *file2;
int coh,coc,con,coo;
long int i,j,k,i1,i2,i3,last;
int numb;
char *tmpi;
char label[10];
double ppmi, ppma,v1[6];
double occ[2000],vir[2000],miin,maax;
int no,nv;
int yes1,yes2,yes3,yes4,conv,scf;
char report[50000];
char buff[10000];
double dou,doudou;
double coup[200][200];
int atom[500];
char tmpc[1000];
char tmpd[1000];
double xyz[500][13];
int mincount,count=0;
char s1[100];
char s2[100];
char s3[100];
char s4[100];
char s5[100];
char s6[100];
char s7[100];
char s8[100];
char s9[100];
char s10[100];
char s11[100];
char s12[100];
char s13[100];
char s14[100];
char s15[100];
char s16[100];
char s17[100];
double len,en,minen,nu1,nu2,dipole;
int stnum,maxnum,sanp,maxpoin,lsanp;
int okall,ok1,ok2,ok3,ok4;



int main(int argc, char *argv[])
{
  while (--argc > 0){
  ok1=0;
  lsanp=1;
  minen=1e6;
  conv=0;
  numb=0;no=0;nv=0;
  dou=0;doudou=0;for(i=0;i>200;i++){occ[i]=0;vir[i]=0;}miin=1e50;maax=-1e50;
      strcpy(tmpc,argv[argc]);
      file=fopen(tmpc,"r");
      if(file!=0){
          while( fgets( buff, sizeof buff, file ) != NULL ) {
  	  if(strstr(buff,"X           Y           Z")){
 		count++;
//		fprintf(stderr,"geom %d\n",count);
	        if ( fgets( buff, sizeof buff, file ) == NULL ) break; //this like is a separator
	           	for(i=0;i<500;i++){
	      	  		if ( fgets( buff, sizeof buff, file ) == NULL ) break;
				xyz[i][10]=xyz[i][7]; xyz[i][11]=xyz[i][8]; xyz[i][12]=xyz[i][9];// prevous matrix
				xyz[i][ 7]=xyz[i][4]; xyz[i][ 8]=xyz[i][5]; xyz[i][9]=xyz[i][6];// prevous matrix
     			    	if ( sscanf( buff, "%ld %d %ld %lf %lf %lf",&i1,&atom[i],&i3,&xyz[i][4], &xyz[i][5], &xyz[i][6]) != 6) break;
   				last=i+1;
			}
	 }
           if(strstr(buff,"on scan point")){
                  sscanf(buff,"%s %s %d %s %s %s %s %s %d %s %s %s %d %s %s %d",s1,s2,&stnum,s3,s4,s5,s6,s7,&maxnum,s8,s9,s10,&sanp,s11,s12,&maxpoin);
		  if(sanp!=lsanp){
//		fprintf(stderr,"save geom %d\n",lsanp);
	  	strcpy(tmpd,tmpc);
		sprintf(tmpd,"%s_s%d.XYZt",tmpc,lsanp-1);
      		file2=fopen(tmpd,"w");
  		for(i=0;i<last;i++){
			fprintf(file2,"%d    %10.6lf %10.6lf %10.6lf\n",atom[i],xyz[i][10],xyz[i][11],xyz[i][12]);
		}
			fprintf(file2,"\n");
		fclose(file2);
		sprintf(tmpd,"%s_s%d.hartree",tmpc,lsanp-1);
      		file2=fopen(tmpd,"w");
		fprintf(file2,"%.8f\n",len);
		fprintf(stderr,"%20.8f\n",len);
		fclose(file2);
		sprintf(tmpd,"%s_s%d.kcal",tmpc,lsanp-1);
      		file2=fopen(tmpd,"w");
		fprintf(file2,"%.8f\n",len*627.509);
		fclose(file2);
		  	lsanp=sanp;
		  }
             }
	     // do the same as above when finishing with "cpu time"
	     //
           if(strstr(buff,"Dipole moment")){
	  	strcpy(tmpd,tmpc);
		strcat(tmpd,".dipole");
      		file2=fopen(tmpd,"w");
          	if ( fgets( buff, sizeof buff, file ) == NULL ) break;
		 		sscanf(buff+0, "%s %lf %s %lf %s %lf %s %lf",s1,&len,s1,&len,s1,&len,s1,&dipole);
		fprintf(file2,"%.8f\n",dipole);
		fclose(file2);
	}

           if(strstr(buff,"Mulliken atomic charges:")){
	  	strcpy(tmpd,tmpc);
		strcat(tmpd,".Mull");
      		file2=fopen(tmpd,"w");
         	if ( fgets( buff, sizeof buff, file ) == NULL ) break;
  		for(i=0;i<last;i++){
          		if ( fgets( buff, sizeof buff, file ) == NULL ) break;
		    		sscanf(buff+10, " %lf ",&len);
		fprintf(file2,"%.8f\n",len);
		}
		fclose(file2);
	}


           if(strstr(buff,"cpu time")){
//		fprintf(stderr,"save geom %d\n",lsanp);
	  	strcpy(tmpd,tmpc);
		sprintf(tmpd,"%s_s%d.XYZt",tmpc,lsanp-1);
      		file2=fopen(tmpd,"w");
  		for(i=0;i<last;i++){
			fprintf(file2,"%d    %10.6lf %10.6lf %10.6lf\n",atom[i],xyz[i][10],xyz[i][11],xyz[i][12]);
		}
			fprintf(file2,"\n");
		fclose(file2);
		sprintf(tmpd,"%s_s%d.hartree",tmpc,lsanp-1);
      		file2=fopen(tmpd,"w");
		fprintf(file2,"%.8f\n",len);
		fprintf(stderr,"%20.8f\n",len);
		fclose(file2);
		sprintf(tmpd,"%s_s%d.kcal",tmpc,lsanp-1);
      		file2=fopen(tmpd,"w");
		fprintf(file2,"%.8f\n",len*627.509);
		fclose(file2);
		  	lsanp=sanp;
		}
	     // do the same as above when finishing with "cpu time"
	     //
	     //


	  if(strstr(buff,"Fermi Contact (FC) contribution to K (Hz):"))
;
	  if(strstr(buff,"Fermi Contact (FC) contribution to J (Hz):"))
;
	  if(strstr(buff,"Spin-dipolar (SD) contribution to K (Hz):"))
;
	  if(strstr(buff,"Spin-dipolar (SD) contribution to J (Hz):"))
;
	  if(strstr(buff,"Paramagnetic spin-orbit (PSO) contribution to K (Hz):"))
;
	  if(strstr(buff,"Paramagnetic spin-orbit (PSO) contribution to J (Hz):"))
;
	  if(strstr(buff,"Diamagnetic spin-orbit (DSO) contribution to K (Hz):"))
;
	  if(strstr(buff,"Diamagnetic spin-orbit (DSO) contribution to J (Hz):"))
;
	  if(strstr(buff,"Total nuclear spin-spin coupling K (Hz):"))
;
	  if(strstr(buff,"Total nuclear spin-spin coupling J (Hz):"))
;
	  if(strstr(buff,"occ")){
		    		no+=sscanf(buff+29, "%lf %lf %lf %lf %lf",&occ[no],&occ[no+1],&occ[no+2],&occ[no+3],&occ[no+4]);
	  }
	  if(strstr(buff,"virt")){
		    		nv+=sscanf(buff+29, "%lf %lf %lf %lf %lf",&vir[nv],&vir[nv+1],&vir[nv+2],&vir[nv+3],&vir[nv+4]);
	  }
	  if(strstr(buff,"Total nuclear spin-spin coupling J (Hz):")){
	  	strcpy(tmpd,tmpc);
		strcat(tmpd,".txtj");
      		file2=fopen(tmpd,"w");
		  for(i=numb;i>0;i-=5){
          		if ( fgets( buff, sizeof buff, file ) == NULL ) break;
			for(j=numb-i;j<numb;j++){
          			if ( fgets( buff, sizeof buff, file ) == NULL ) break;
				tmpi=buff;
	//			fprintf(stdout,"%s\n",buff);
				for(k=0;*tmpi!=0;k++) {
	//				fprintf(stdout,"%d ",*tmpi);
					if(*tmpi==68) {// correct D for E for floating point format
						*tmpi=69; 
					}
					tmpi++;
				}
	//			fprintf(stdout,"%s\n",buff);
		    		sscanf(buff, "%*d %lf %lf %lf %lf %lf",&v1[0],&v1[1],&v1[2],&v1[3],&v1[4]);
				for(k=0;k<(j+1-(numb-i));k++) {
					if(k==5) break;
					coup[j][k+numb-i]=v1[k];
					coup[k+numb-i][j]=v1[k];
//					fprintf(stderr, "%2d-%2d ",j+1,k+1+numb-i);
				}
			}
		  }
		  for(i=0;i<numb;i++){
		  	for(j=0;j<numb;j++) fprintf(file2,"%10.6le ",coup[i][j]);
			fprintf(file2,"\n");
		}
		fclose(file2);
	 }


	  
	  if(strstr(buff,"SCF GIAO Magnetic shielding tensor")){
	  	strcpy(tmpd,tmpc);
		strcat(tmpd,".txt");
      		file2=fopen(tmpd,"w");
		for(i=0;i<10000000;i++){
          		if ( fgets( buff, sizeof buff, file ) == NULL ) break;
		    	if(sscanf(buff, "%d %s %*s %*s %lf %*s %*s %lf",&numb,label,&ppmi,&ppma) < 4) break;
			if(strstr(label,"H")) fprintf(file2,"%d %s %.6lf\n",numb,label,-ppmi+31.90376667);
				else
				{if(strstr(label,"C")) fprintf(file2,"%d %s %.6lf\n",numb,label,-ppmi+183.4183);
					else
					fprintf(file2,"%d %s %.6lf\n",numb,label,ppmi-183.4183);}
          		if ( fgets( buff, sizeof buff, file ) == NULL ) break;
          		if ( fgets( buff, sizeof buff, file ) == NULL ) break;
          		if ( fgets( buff, sizeof buff, file ) == NULL ) break;
          		if ( fgets( buff, sizeof buff, file ) == NULL ) break;
		}
		fclose(file2);
	  	strcpy(tmpd,tmpc);
		strcat(tmpd,".nb");
      		file2=fopen(tmpd,"w");
			fprintf(file2,"%ld\n",i);
		fclose(file2);
	  }
	  if(strstr(buff,"Low frequenc")){
		sprintf(tmpd,"%40s : %s",tmpc,buff);
	  	strcat(report,tmpd);
	  }
	  if(strstr(buff,"maginary frequ")){
		sprintf(tmpd,"%40s : %s",tmpc,buff);
	  	strcat(report,tmpd);
	  }
	  if(strstr(buff,"Maximum Force        ")){conv=1; sscanf( buff+50, "%s",tmpd);if(strstr(tmpd,"YES")) yes1=1; else yes1=0; }
	  if(strstr(buff,"Maximum Displacement ")){conv=1; sscanf( buff+50, "%s",tmpd);if(strstr(tmpd,"YES")) yes2=1; else yes2=0; }
	  if(strstr(buff,"RMS     Force        ")){conv=1; sscanf( buff+50, "%s",tmpd);if(strstr(tmpd,"YES")) yes3=1; else yes3=0; }
	  if(strstr(buff,"RMS     Displacement ")){conv=1; sscanf( buff+50, "%s",tmpd);if(strstr(tmpd,"YES")) yes4=1; else yes4=0; }
	  if(strstr(buff,"SCF Done")){
    		  // sscanf( buff+26 ,"%lf",&doudou);scf=1;
		  len=en;
		   sscanf(buff,"%s %s %s %s %lf %s %s ",s1,s2,s3,s4,&en,s6,s7);
		   doudou=en;scf=1;
		   if(en<minen){
		//fprintf(stderr,"en %lf\n",en);
  			minen=en;mincount=count;
	  		for(i=0;i<last;i++){
				okall=ok1*ok2*ok3*ok4;
				xyz[i][1]=xyz[i][4];
				xyz[i][2]=xyz[i][5];
				xyz[i][3]=xyz[i][6];
			}
		}
	   }
	  if(strstr(buff,"Sum of electronic and zero-point Energies")){
		   sscanf( buff+45, "%lf",&dou); sprintf(tmpd,"%40s : Energy after zero point correction : %12.7lf  or %12.7lf kcal\n",tmpc,dou,dou*627.509);strcat(report,tmpd);
	   }
	  /*	for(i=0;i<500;i++){
          		if ( fgets( buff, sizeof buff, file ) == NULL ) break;
		    	if ( sscanf( buff, "%d %d %d %lf %lf %lf",&i1,&atom[i],&i3,&xyz[i][1], &xyz[i][2], &xyz[i][3]) != 6) break;
//			fprintf(stderr,"%d    %10.6lf %10.6lf %10.6lf\n",atom[i],xyz[i][1],xyz[i][2],xyz[i][3]);
			last=i+1;
		}*/
	 }
	if(conv){
		if(yes1+yes2+yes3+yes4 == 4) 
			sprintf(tmpd,"%40s : Convergence OK!\n",tmpc);
		else
			sprintf(tmpd,"%40s : Convergence **** NOT **** Not OK!\n",tmpc);
		strcat(report,tmpd);
	}
	if(scf)   sprintf(tmpd,"%40s : SCF : %12.7lf  or %12.7lf kcal\n",tmpc,minen,minen*627.509);strcat(report,tmpd);
	sprintf(tmpd,"\n",tmpc); strcat(report,tmpd);
	if(scf) if(minen<1e16) {
	  	strcpy(tmpd,tmpc);
		strcat(tmpd,".hartree");
      		file2=fopen(tmpd,"w");
		fprintf(file2,"%20.8lf\n",minen);
		fclose(file2);
	  	strcpy(tmpd,tmpc);
		strcat(tmpd,".kcal");
      		file2=fopen(tmpd,"w");
		fprintf(file2,"%.8lf\n",minen*627.509);
		fclose(file2);
	  	strcpy(tmpd,tmpc);
		strcat(tmpd,".min.XYZt");
      		file2=fopen(tmpd,"w");
		coc=0,coh=0,coo=0,con=0;
  		for(i=0;i<last;i++){
			fprintf(file2,"%d    %10.6lf %10.6lf %10.6lf\n",atom[i],xyz[i][1],xyz[i][2],xyz[i][3]);
			if(atom[i]==6)coc++;
			if(atom[i]==1)coh++;
			if(atom[i]==8)coo++;
			if(atom[i]==7)con++;
		}
		fprintf(file2,"\n");
		fclose(file2);
		fprintf(stderr,"formula : C%d H%d O%d N%d \n",coc,coh,coo,con);
	}

		/* energies of occ and virt orbitals */
	if(no+nv>0){
	  	strcpy(tmpd,tmpc);
		strcat(tmpd,".txte");
      		file2=fopen(tmpd,"w");
	  	fprintf(file2,"%d ",no); 
		for(i=0;i<no;i++){ fprintf(file2,"%lf ",occ[i]);if(miin>occ[i])miin=occ[i];if(maax<occ[i])maax=occ[i];}
	  	fprintf(file2,"\n");
	  	fprintf(file2,"%d ",nv); 
		for(i=0;i<nv;i++){ fprintf(file2,"%lf ",vir[i]);if(miin>vir[i])miin=vir[i];if(maax<vir[i])maax=vir[i];}
	  	fprintf(file2,"\n");
		fclose(file2);
	  	strcpy(tmpd,tmpc);
		strcat(tmpd,".txteminmax");
      		file2=fopen(tmpd,"w");
		fprintf(file2,"%lf\n",miin);
		fprintf(file2,"%lf\n",maax);
		fclose(file2);
	}
	}
	fclose(file);
   }
   printf("%s",report);
}
/*    for(k=1;k<13;k++)
    for(l=1;l<32;l++){
		printf("%2d %2d ",k,l); for(i=0;i<24;i++){tm=0;for(j=0;j<60;j++) tm+=table[k][l][i*60+j];if(tm>9999) printf(" xxxx");else printf("%5d",tm);} printf("\n");
		}*/
