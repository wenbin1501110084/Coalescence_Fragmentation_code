#include<iostream>
#include<cmath>
#include<fstream>
#include<string>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<time.h>
#include <math.h>
#define random(x) (rand()%x)
#define ptcutm 3.0

using namespace std;

int findarray(int*p,int len,int val);
int main( int argc, char* argv[])
{
	int nrun = atoi(argv[1]);
	//string path1=string(argv[2]);
	//string path2=string(argv[3]);
	/*
	char infilefrag1[128];
        char infilefrag2[128];
        char infilefrag3[128];
        char infilefrag4[128];
        char infilefrag5[128];
	*/
                srand((double)clock());
                srand(unsigned(time(0)));

	char infilethermal[128];
	char infilejet[128];
	char outfiles[128];
	int id, mid,thermalid[50000]={0},thermalcharge[50000]={0};

	int Netherml,Nejet,mark,Nfrag1,Nfrag2,Nfrag3,Nfrag4,Nfrag5,Nfrag6,Nfrag7;
	double thermalpx[50000]={0.0},thermalpy[50000]={0.0},thermalpz[50000]={0.0},thermalEng[50000]={0.0},thermalmass[50000]={0.0},thermalx[50000]={0.0},thermaly[50000]={0.0},thermalz[50000]={0.0},thermalt[50000]={0.0};
       // double px2[50000]={0.0},py2[50000]={0.0},pz2[50000]={0.0},energy2[50000]={0.0},mass2[50000]={0.0},xx2[50000]={0.0},yy2[50000]={0.0},zz2[50000]={0.0},tt2[50000]={0.0};
        double px,py,pz,energy,mass,xx,yy,zz,tt;
	srand((double)clock());
	/*SSSSSSSS
	double* thermalid;
	double* thermalEng;
	double* thermalpx;
	double* thermalpy;
	double* thermalpz;
	double* thermalmass;
	double* thermalx;
	double* thermaly;
	double* thermalz;
	double* thermalt;
	*/
        double masspi=0.13957,masspi0=0.13498;
        double massK=0.49368,massK0=0.49765;
        double massp=0.93827,massn=0.93957;
        double massLambda=1.11568,massphi=1.01946,massOmega=1.67243;

        sprintf(outfiles,"oscar.dat");
	ofstream outfile(outfiles);
	//pgd particle
        char infilepdg[128];
	int pdgid[600]={0};
        sprintf(infilepdg,"chosen_particles.dat");
        ifstream inhypdg(infilepdg);
	int lenght=251;
	for(int ll=0;ll<lenght;ll++){
		inhypdg>>pdgid[ll];
	}

        char ipput_filename1[128];
        sprintf(ipput_filename1,"./hadrons_frag1.dat");
        ifstream infilef1(ipput_filename1);

        char ipput_filename2[128];
        sprintf(ipput_filename2,"./hadrons_frag2.dat");
        ifstream infilef2(ipput_filename2);

        char ipput_filename3[128];
        sprintf(ipput_filename3,"./hadrons_frag3.dat");
        ifstream infilef3(ipput_filename3);

        char ipput_filename4[128];
        sprintf(ipput_filename4,"./hadrons_frag4.dat");
        ifstream infilef4(ipput_filename4);

        char ipput_filename5[128];
        sprintf(ipput_filename5,"./hadrons_frag5.dat");
        ifstream infilef5(ipput_filename5);

        char ipput_filename6[128];
        sprintf(ipput_filename6,"./hadrons_frag6.dat");
        ifstream infilef6(ipput_filename6);

        char ipput_filename7[128];
        sprintf(ipput_filename7,"./hadrons_frag7.dat");
        ifstream infilef7(ipput_filename7);


/*
	string ipput_filename1;
	ipput_filename1 = "../hadrons_frag1.dat";// oooooooooooooutput files of final hadrons
        sprintf(infilefrag1,  ipput_filename1.c_str());

        string ipput_filename2;
        ipput_filename2 = "../hadrons_frag2.dat";// oooooooooooooutput files of final hadrons
        sprintf(infilefrag2,  ipput_filename2.c_str());

        string ipput_filename3;
        ipput_filename3 = "../hadrons_frag3.dat";// oooooooooooooutput files of final hadrons
        sprintf(infilefrag3,  ipput_filename3.c_str());

        string ipput_filename4;
        ipput_filename4 = "../hadrons_frag4.dat";// oooooooooooooutput files of final hadrons
        sprintf(infilefrag4,  ipput_filename4.c_str());

        string ipput_filename5;
        ipput_filename5 = "../hadrons_frag5.dat";// oooooooooooooutput files of final hadrons
        sprintf(infilefrag5,  ipput_filename5.c_str());
i/
        FILE* infilef1;
        infilef1 = fopen(infilefrag1,"r");

        FILE* infilef2;
        infilef2 = fopen(infilefrag2,"r");

        FILE* infilef3;
        infilef3 = fopen(infilefrag3,"r");

        FILE* infilef4;
        infilef4 = fopen(infilefrag4,"r");

        FILE* infilef5;
        infilef5 = fopen(infilefrag5,"r");
*/

        string ipput_filenamec;
        ipput_filenamec = "./coaleced_hadron.dat";// oooooooooooooutput files of final hadrons

        sprintf(infilejet, ipput_filenamec.c_str());
        FILE* infilec;
        infilec = fopen(infilejet,"r");

//	outfile<<"OSC1997A"<<endl;
//	outfile<<"final_id_p_x"<<endl;
//	outfile<<" 3DHydro       1.1  (197,    79)+(197,    79)  eqsp  0.1000E+03         1"<<endl;
	int totalnum=0;
	for(int ie=0;ie!=nrun;ie++){
                char infilehypersf[128];
                sprintf(infilehypersf,"./mc_particle_list");
                ifstream inhy(infilehypersf);
                fscanf(infilec,"%d %d\n",&mid,&Nejet);

		if(!infilef1.eof()){infilef1>>mid>>Nfrag1;}else{Nfrag1=0;}
                if(!infilef2.eof()){infilef2>>mid>>Nfrag2;}else{Nfrag2=0;}
                if(!infilef3.eof()){infilef3>>mid>>Nfrag3;}else{Nfrag3=0;}
                if(!infilef4.eof()){infilef4>>mid>>Nfrag4;}else{Nfrag4=0;}
                if(!infilef5.eof()){infilef5>>mid>>Nfrag5;}else{Nfrag5=0;}
                if(!infilef6.eof()){infilef6>>mid>>Nfrag6;}else{Nfrag6=0;}
                if(!infilef7.eof()){infilef7>>mid>>Nfrag7;}else{Nfrag7=0;}

		/*
                fscanf(infilef1,"%d %d\n",&mid,&Nfrag1);
                fscanf(infilef2,"%d %d\n",&mid,&Nfrag2);
                fscanf(infilef3,"%d %d\n",&mid,&Nfrag3);
                fscanf(infilef4,"%d %d\n",&mid,&Nfrag4);
                fscanf(infilef5,"%d %d\n",&mid,&Nfrag5);
		*/
                //getline(inhy,stempa);
                int count=0;
                /*
		while(!inhy.eof()){
                	inhy>>thermalt[count]>>thermalx[count]>>thermaly[count]>>thermalz[count]>>thermalmass[count]>>thermalEng[count]>>thermalpx[count]>>thermalpy[count]>>thermalpz[count]>>thermalid[count]>>mid>>thermalcharge[count];
			int ss1=findarray(pdgid,lenght, thermalid[count]);
			int ss2=1; 
			double ptmeson2=thermalpx[count]*thermalpx[count]+thermalpy[count]*thermalpy[count];
			if((abs(thermalid[count])<1000)&&(ptmeson2>ptcutm*ptcutm)){ss2=0;}
                        if((abs(thermalid[count])>1000)&&(ptmeson2>4.50*4.50)){ss2=0;}
			int sss=ss1*ss2;
                	if(sss== 1){count++;}//else{cout<<"dddddddd"<<endl;}
                }
                */
		if(Nfrag1>0){
                for(int ij=1;ij!=Nfrag1+1;ij++){
                        //fscanf(infilef1,"%d %lf %lf %lf %lf %lf %lf %lf  %lf %lf\n", &id, &px, &py, &pz, &energy, &mass,&xx, &yy, &zz, &tt);
                        //if(isnan(id)||isnan(px)|| isnan(py)||isnan(pz)||isnan(energy)||isnan(xx)||isnan(yy)||isnan(zz)||isnan(tt)){
                        //xx=40.0;yy=40.0; zz=40.0;tt=40.0;
                        //}
                        infilef1>>id>>px>>py>>pz>>energy>>mass>>xx>>yy>>zz>>tt;
                        int sss=findarray(pdgid,lenght, id);
                	if(sss== 1){
                        if(abs(id)==130){
                                        double ram=random(1000)/1000.00;
                                        if(ram<0.5){
                                                energy=sqrt(px*px+py*py+pz*pz+massK*massK);
                                                thermalid[count] = 321;
                                                thermalEng[count] = energy;
                                                thermalmass[count] = massK;
                                        }
                                        if((0.5<=ram)&(ram<1.0)){
                                                energy=sqrt(px*px+py*py+pz*pz+massK*massK);
                                                thermalid[count] = -321;
                                                thermalEng[count] = energy;
                                                thermalmass[count] = massK;
                                        }
                                        if((1.0<=ram)&(ram<1.1)){
                                                energy=sqrt(px*px+py*py+pz*pz+massK0*massK0);
                                                thermalid[count] = -311;
                                                thermalEng[count] = energy;
                                                thermalmass[count] = massK0;
                                        }
                                        if(1.0<=ram){
                                                energy=sqrt(px*px+py*py+pz*pz+massK0*massK0);
                                                thermalid[count] = 311;
                                                thermalEng[count] = energy;
                                                thermalmass[count] = massK0;
                                        }
                        }else{
                                thermalid[count] = id;
                                thermalEng[count] = energy;
                                thermalmass[count] = mass;
                        }
                        thermalpx[count] = px;
                        thermalpy[count] = py;
                        thermalpz[count] = pz;
                        thermalx[count] = xx;
                        thermaly[count] = yy;
                        thermalz[count] = zz;
                        thermalt[count] = tt;
                        count++;
		}
		}
                }
                if(Nfrag2>0){
                for(int ij=1;ij!=Nfrag2+1;ij++){
                //        fscanf(infilef2,"%d %lf %lf %lf %lf %lf %lf %lf  %lf %lf\n", &id, &px, &py, &pz, &energy, &mass,&xx, &yy, &zz, &tt);
                //        if(isnan(id)||isnan(px)|| isnan(py)||isnan(pz)||isnan(energy)||isnan(xx)||isnan(yy)||isnan(zz)||isnan(tt)){
                //        xx=40.0;yy=40.0; zz=40.0;tt=40.0;
                //        }
                        infilef2>>id>>px>>py>>pz>>energy>>mass>>xx>>yy>>zz>>tt;
                        int sss=findarray(pdgid,lenght, id);
                	if(sss== 1){
                        if(abs(id)==130){
                                        double ram=random(1000)/1000.00;
                                        if(ram<0.25){
                                                energy=sqrt(px*px+py*py+pz*pz+massK*massK);
                                                thermalid[count] = 321;
                                                thermalEng[count] = energy;
                                                thermalmass[count] = massK;
                                        }
                                        if((0.25<=ram)&(ram<0.50)){
                                                energy=sqrt(px*px+py*py+pz*pz+massK*massK);
                                                thermalid[count] = -321;
                                                thermalEng[count] = energy;
                                                thermalmass[count] = massK;
                                        }
                                        if((0.50<=ram)&(ram<0.750)){
                                                energy=sqrt(px*px+py*py+pz*pz+massK0*massK0);
                                                thermalid[count] = -311;
                                                thermalEng[count] = energy;
                                                thermalmass[count] = massK0;
                                        }
                                        if(0.750<=ram){
                                                energy=sqrt(px*px+py*py+pz*pz+massK0*massK0);
                                                thermalid[count] = 311;
                                                thermalEng[count] = energy;
                                                thermalmass[count] = massK0;
                                        }
                        }else{
                                thermalid[count] = id;
                                thermalEng[count] = energy;
                                thermalmass[count] = mass;
                        }
                        thermalpx[count] = px;
                        thermalpy[count] = py;
                        thermalpz[count] = pz;
                        thermalx[count] = xx;
                        thermaly[count] = yy;
                        thermalz[count] = zz;
                        thermalt[count] = tt;
                        count++;
                }
                }
                }
                if(Nfrag3>0){
                for(int ij=1;ij!=Nfrag3+1;ij++){
                //        fscanf(infilef3,"%d %lf %lf %lf %lf %lf %lf %lf  %lf %lf\n", &id, &px, &py, &pz, &energy, &mass,&xx, &yy, &zz, &tt);
                //        if(isnan(id)||isnan(px)|| isnan(py)||isnan(pz)||isnan(energy)||isnan(xx)||isnan(yy)||isnan(zz)||isnan(tt)){
                //        xx=40.0;yy=40.0; zz=40.0;tt=40.0;
                //        }
                        infilef3>>id>>px>>py>>pz>>energy>>mass>>xx>>yy>>zz>>tt;
                        int sss=findarray(pdgid,lenght, id);
                	if(sss== 1){
                        if(abs(id)==130){
                                        double ram=random(1000)/1000.00;
                                        if(ram<0.25){
                                                energy=sqrt(px*px+py*py+pz*pz+massK*massK);
                                                thermalid[count] = 321;
                                                thermalEng[count] = energy;
                                                thermalmass[count] = massK;
                                        }
                                        if((0.25<=ram)&(ram<0.50)){
                                                energy=sqrt(px*px+py*py+pz*pz+massK*massK);
                                                thermalid[count] = -321;
                                                thermalEng[count] = energy;
                                                thermalmass[count] = massK;
                                        }
                                        if((0.50<=ram)&(ram<0.750)){
                                                energy=sqrt(px*px+py*py+pz*pz+massK0*massK0);
                                                thermalid[count] = -311;
                                                thermalEng[count] = energy;
                                                thermalmass[count] = massK0;
                                        }
                                        if(0.750<=ram){
                                                energy=sqrt(px*px+py*py+pz*pz+massK0*massK0);
                                                thermalid[count] = 311;
                                                thermalEng[count] = energy;
                                                thermalmass[count] = massK0;
                                        }
                        }else{
                                thermalid[count] = id;
                                thermalEng[count] = energy;
                                thermalmass[count] = mass;
                        }
                        thermalpx[count] = px;
                        thermalpy[count] = py;
                        thermalpz[count] = pz;
                        thermalx[count] = xx;
                        thermaly[count] = yy;
                        thermalz[count] = zz;
                        thermalt[count] = tt;
                        count++;
                }
                }
                }
                if(Nfrag4>0){
                for(int ij=1;ij!=Nfrag4+1;ij++){
                        //fscanf(infilef4,"%d %lf %lf %lf %lf %lf %lf %lf  %lf %lf\n", &id, &px, &py, &pz, &energy, &mass,&xx, &yy, &zz, &tt);
                        //if(isnan(id)||isnan(px)|| isnan(py)||isnan(pz)||isnan(energy)||isnan(xx)||isnan(yy)||isnan(zz)||isnan(tt)){
                        //xx=40.0;yy=40.0; zz=40.0;tt=40.0;
                        // }
                        infilef4>>id>>px>>py>>pz>>energy>>mass>>xx>>yy>>zz>>tt;
                        int sss=findarray(pdgid,lenght, id);
                	if(sss== 1){
                        if(abs(id)==130){
                                        double ram=random(1000)/1000.00;
                                        if(ram<0.25){
                                                energy=sqrt(px*px+py*py+pz*pz+massK*massK);
                                                thermalid[count] = 321;
                                                thermalEng[count] = energy;
                                                thermalmass[count] = massK;
                                        }
                                        if((0.25<=ram)&(ram<0.50)){
                                                energy=sqrt(px*px+py*py+pz*pz+massK*massK);
                                                thermalid[count] = -321;
                                                thermalEng[count] = energy;
                                                thermalmass[count] = massK;
                                        }
                                        if((0.50<=ram)&(ram<0.750)){
                                                energy=sqrt(px*px+py*py+pz*pz+massK0*massK0);
                                                thermalid[count] = -311;
                                                thermalEng[count] = energy;
                                                thermalmass[count] = massK0;
                                        }
                                        if(0.750<=ram){
                                                energy=sqrt(px*px+py*py+pz*pz+massK0*massK0);
                                                thermalid[count] = 311;
                                                thermalEng[count] = energy;
                                                thermalmass[count] = massK0;
                                        }
                        }else{
                                thermalid[count] = id;
                                thermalEng[count] = energy;
                                thermalmass[count] = mass;
                        }
                        thermalpx[count] = px;
                        thermalpy[count] = py;
                        thermalpz[count] = pz;
                        thermalx[count] = xx;
                        thermaly[count] = yy;
                        thermalz[count] = zz;
                        thermalt[count] = tt;
                        count++;
                }
                }
                }
                if(Nfrag5>0){
                for(int ij=1;ij!=Nfrag5+1;ij++){
                        //fscanf(infilef5,"%d %lf %lf %lf %lf %lf %lf %lf  %lf %lf\n", &id, &px, &py, &pz, &energy, &mass,&xx, &yy, &zz, &tt);
                        //if(isnan(id)||isnan(px)|| isnan(py)||isnan(pz)||isnan(energy)||isnan(xx)||isnan(yy)||isnan(zz)||isnan(tt)){
                        //xx=40.0;yy=40.0; zz=40.0;tt=40.0;
                        // }
                       infilef5>>id>>px>>py>>pz>>energy>>mass>>xx>>yy>>zz>>tt;
                        int sss=findarray(pdgid,lenght, id);
                	if(sss== 1){
                        if(abs(id)==130){
                                        double ram=random(1000)/1000.00;
                                        if(ram<0.25){
                                                energy=sqrt(px*px+py*py+pz*pz+massK*massK);
                                                thermalid[count] = 321;
                                                thermalEng[count] = energy;
                                                thermalmass[count] = massK;
                                        }
                                        if((0.25<=ram)&(ram<0.50)){
                                                energy=sqrt(px*px+py*py+pz*pz+massK*massK);
                                                thermalid[count] = -321;
                                                thermalEng[count] = energy;
                                                thermalmass[count] = massK;
                                        }
                                        if((0.50<=ram)&(ram<0.750)){
                                                energy=sqrt(px*px+py*py+pz*pz+massK0*massK0);
                                                thermalid[count] = -311;
                                                thermalEng[count] = energy;
                                                thermalmass[count] = massK0;
                                        }
                                        if(0.750<=ram){
                                                energy=sqrt(px*px+py*py+pz*pz+massK0*massK0);
                                                thermalid[count] = 311;
                                                thermalEng[count] = energy;
                                                thermalmass[count] = massK0;
                                        }
                        }else{
                                thermalid[count] = id;
                                thermalEng[count] = energy;
                                thermalmass[count] = mass;
                        }
                        thermalpx[count] = px;
                        thermalpy[count] = py;
                        thermalpz[count] = pz;
                        thermalx[count] = xx;
                        thermaly[count] = yy;
                        thermalz[count] = zz;
                        thermalt[count] = tt;
                        count++;
                }
                }
                }

                if(Nfrag6>0){
                for(int ij=1;ij!=Nfrag6+1;ij++){
                        //fscanf(infilef5,"%d %lf %lf %lf %lf %lf %lf %lf  %lf %lf\n", &id, &px, &py, &pz, &energy, &mass,&xx, &yy, &zz, &tt);
                        //if(isnan(id)||isnan(px)|| isnan(py)||isnan(pz)||isnan(energy)||isnan(xx)||isnan(yy)||isnan(zz)||isnan(tt)){
                        //xx=40.0;yy=40.0; zz=40.0;tt=40.0;
                        // }
                       infilef6>>id>>px>>py>>pz>>energy>>mass>>xx>>yy>>zz>>tt;
                        int sss=findarray(pdgid,lenght, id);
                	if(sss== 1){
                        if(abs(id)==130){
                                        double ram=random(1000)/1000.00;
                                        if(ram<0.25){
                                                energy=sqrt(px*px+py*py+pz*pz+massK*massK);
                                                thermalid[count] = 321;
                                                thermalEng[count] = energy;
                                                thermalmass[count] = massK;
                                        }
                                        if((0.25<=ram)&(ram<0.50)){
                                                energy=sqrt(px*px+py*py+pz*pz+massK*massK);
                                                thermalid[count] = -321;
                                                thermalEng[count] = energy;
                                                thermalmass[count] = massK;
                                        }
                                        if((0.50<=ram)&(ram<0.750)){
                                                energy=sqrt(px*px+py*py+pz*pz+massK0*massK0);
                                                thermalid[count] = -311;
                                                thermalEng[count] = energy;
                                                thermalmass[count] = massK0;
                                        }
                                        if(0.750<=ram){
                                                energy=sqrt(px*px+py*py+pz*pz+massK0*massK0);
                                                thermalid[count] = 311;
                                                thermalEng[count] = energy;
                                                thermalmass[count] = massK0;
                                        }
                        }else{
                                thermalid[count] = id;
                                thermalEng[count] = energy;
                                thermalmass[count] = mass;
                        }
                        thermalpx[count] = px;
                        thermalpy[count] = py;
                        thermalpz[count] = pz;
                        thermalx[count] = xx;
                        thermaly[count] = yy;
                        thermalz[count] = zz;
                        thermalt[count] = tt;
                        count++;
                }
                }
                }

                if(Nfrag7>0){
                for(int ij=1;ij!=Nfrag7+1;ij++){
                        //fscanf(infilef5,"%d %lf %lf %lf %lf %lf %lf %lf  %lf %lf\n", &id, &px, &py, &pz, &energy, &mass,&xx, &yy, &zz, &tt);
                        //if(isnan(id)||isnan(px)|| isnan(py)||isnan(pz)||isnan(energy)||isnan(xx)||isnan(yy)||isnan(zz)||isnan(tt)){
                        //xx=40.0;yy=40.0; zz=40.0;tt=40.0;
                        // }
                       infilef7>>id>>px>>py>>pz>>energy>>mass>>xx>>yy>>zz>>tt;
                        int sss=findarray(pdgid,lenght, id);
                	if(sss== 1){
                        if(abs(id)==130){
                                        double ram=random(1000)/1000.00;
                                        if(ram<0.25){
                                                energy=sqrt(px*px+py*py+pz*pz+massK*massK);
                                                thermalid[count] = 321;
                                                thermalEng[count] = energy;
                                                thermalmass[count] = massK;
                                        }
                                        if((0.25<=ram)&(ram<0.50)){
                                                energy=sqrt(px*px+py*py+pz*pz+massK*massK);
                                                thermalid[count] = -321;
                                                thermalEng[count] = energy;
                                                thermalmass[count] = massK;
                                        }
                                        if((0.50<=ram)&(ram<0.750)){
                                                energy=sqrt(px*px+py*py+pz*pz+massK0*massK0);
                                                thermalid[count] = -311;
                                                thermalEng[count] = energy;
                                                thermalmass[count] = massK0;
                                        }
                                        if(0.750<=ram){
                                                energy=sqrt(px*px+py*py+pz*pz+massK0*massK0);
                                                thermalid[count] = 311;
                                                thermalEng[count] = energy;
                                                thermalmass[count] = massK0;
                                        }
                        }else{
                                thermalid[count] = id;
                                thermalEng[count] = energy;
                                thermalmass[count] = mass;
                        }
                        thermalpx[count] = px;
                        thermalpy[count] = py;
                        thermalpz[count] = pz;
                        thermalx[count] = xx;
                        thermaly[count] = yy;
                        thermalz[count] = zz;
                        thermalt[count] = tt;
                        count++;
                }
                }
                }


                srand((double)clock());
                srand(unsigned(time(0)));

		//coalescence hadron
		for(int ij=1;ij!=Nejet+1;ij++){
			fscanf(infilec,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %d\n", &mid, &id, &px, &py, &pz, &energy, &xx, &yy, &zz, &tt,&mark);
//			infile2>>mid>>id>>px>>py>>pz>>energy>>xx>>yy>>zz>>tt;
			if(isnan(id)||isnan(px)|| isnan(py)||isnan(pz)||isnan(energy)||isnan(xx)||isnan(yy)||isnan(zz)||isnan(tt)||isnan(mark)){
				xx=40.0;yy=40.0; zz=40.0;tt=40.0;
			}
			// get the real hadron id of jet hadron and coalescence hadron 
			//pion
			//if(mark!=0){//CCCCCCCut
				if((abs(id)==1)||(abs(id)==11)){	
					double ram=random(1000)/1000.00;
	//				cout<<"ram0= "<<ram<<endl;
					if(ram<0.333){
						energy=sqrt(px*px+py*py+pz*pz+masspi*masspi);
						thermalid[count] = 211;
						thermalEng[count] = energy;
						thermalmass[count] = masspi;

						thermalpx[count] = px;
						thermalpy[count] = py;
						thermalpz[count] = pz;
						thermalx[count] = xx;
						thermaly[count] = yy;
						thermalz[count] = zz;
						thermalt[count] = tt;
						count++;
					}
					if((0.333<=ram) & (ram<0.666)){
						energy=sqrt(px*px+py*py+pz*pz+masspi*masspi);
						thermalid[count] = -211;
						thermalEng[count] = energy;
						thermalmass[count] = masspi;

						thermalpx[count] = px;
						thermalpy[count] = py;
						thermalpz[count] = pz;
						thermalx[count] = xx;
						thermaly[count] = yy;
						thermalz[count] = zz;
						thermalt[count] = tt;
						count++;
					}
					if(0.666<=ram){
						energy=sqrt(px*px+py*py+pz*pz+masspi0*masspi0);
						thermalid[count] = 111;
						thermalEng[count] = energy;
						thermalmass[count] = masspi0;

						thermalpx[count] = px;
						thermalpy[count] = py;
						thermalpz[count] = pz;
						thermalx[count] = xx;
						thermaly[count] = yy;
						thermalz[count] = zz;
						thermalt[count] = tt;
						count++;
					}
				}
				//kaon
				if((abs(id)==2)||(abs(id)==22)){	
					double ram=random(1000)/1000.00;
	//				cout<<"ram1= "<<ram<<endl;
					if(ram<0.25){
						energy=sqrt(px*px+py*py+pz*pz+massK*massK);
						thermalid[count] = 321;
						thermalEng[count] = energy;
						thermalmass[count] = massK;

						thermalpx[count] = px;
						thermalpy[count] = py;
						thermalpz[count] = pz;
						thermalx[count] = xx;
						thermaly[count] = yy;
						thermalz[count] = zz;
						thermalt[count] = tt;
						count++;
					}
					if((0.25<=ram)&(ram<0.50)){
						energy=sqrt(px*px+py*py+pz*pz+massK*massK);
						thermalid[count] = -321;
						thermalEng[count] = energy;
						thermalmass[count] = massK;

						thermalpx[count] = px;
						thermalpy[count] = py;
						thermalpz[count] = pz;
						thermalx[count] = xx;
						thermaly[count] = yy;
						thermalz[count] = zz;
						thermalt[count] = tt;
						count++;
					}
					if((0.50<=ram)&(ram<0.750)){
						energy=sqrt(px*px+py*py+pz*pz+massK0*massK0);
						thermalid[count] = -311;
						thermalEng[count] = energy;
						thermalmass[count] = massK0;

						thermalpx[count] = px;
						thermalpy[count] = py;
						thermalpz[count] = pz;
						thermalx[count] = xx;
						thermaly[count] = yy;
						thermalz[count] = zz;
						thermalt[count] = tt;
						count++;
					}
					if(0.750<=ram){
						energy=sqrt(px*px+py*py+pz*pz+massK0*massK0);
						thermalid[count] = 311;
						thermalEng[count] = energy;
						thermalmass[count] = massK0;

						thermalpx[count] = px;
						thermalpy[count] = py;
						thermalpz[count] = pz;
						thermalx[count] = xx;
						thermaly[count] = yy;
						thermalz[count] = zz;
						thermalt[count] = tt;
						count++;
					}
				}
				//nucleon(p,n)
				if((abs(id)==3)||(abs(id)==33)){	
					double ram=random(1000)/1000.00;
	//				cout<<"ram2= "<<ram<<endl;
					if(ram<0.25){
						energy=sqrt(px*px+py*py+pz*pz+massp*massp);
						thermalid[count] = 2212;
						thermalEng[count] = energy;
						thermalmass[count] = massp;

						thermalpx[count] = px;
						thermalpy[count] = py;
						thermalpz[count] = pz;
						thermalx[count] = xx;
						thermaly[count] = yy;
						thermalz[count] = zz;
						thermalt[count] = tt;
						count++;
					}
					if((0.25<=ram)&(ram<0.50)){
						energy=sqrt(px*px+py*py+pz*pz+massp*massp);
						thermalid[count] = -2212;
						thermalEng[count] = energy;
						thermalmass[count] = massp;

						thermalpx[count] = px;
						thermalpy[count] = py;
						thermalpz[count] = pz;
						thermalx[count] = xx;
						thermaly[count] = yy;
						thermalz[count] = zz;
						thermalt[count] = tt;
						count++;
					}
					if((0.50<=ram)&(ram<0.750)){
						energy=sqrt(px*px+py*py+pz*pz+massn*massn);
						thermalid[count] = -2112;
						thermalEng[count] = energy;
						thermalmass[count] = massn;

						thermalpx[count] = px;
						thermalpy[count] = py;
						thermalpz[count] = pz;
						thermalx[count] = xx;
						thermaly[count] = yy;
						thermalz[count] = zz;
						thermalt[count] = tt;
						count++;
					}
					if(0.750<=ram){
						energy=sqrt(px*px+py*py+pz*pz+massn*massn);
						thermalid[count] = 2112;
						thermalEng[count] = energy;
						thermalmass[count] = massn;

						thermalpx[count] = px;
						thermalpy[count] = py;
						thermalpz[count] = pz;
						thermalx[count] = xx;
						thermaly[count] = yy;
						thermalz[count] = zz;
						thermalt[count] = tt;
						count++;
					}
				}
				//Lambda
				if((abs(id)==4)||(abs(id)==44)){	
					double ram=random(1000)/1000.00;
					if(ram<0.5){
						energy=sqrt(px*px+py*py+pz*pz+massLambda*massLambda);
						thermalid[count] = 3122;
						thermalEng[count] = energy;
						thermalmass[count] = massLambda;

						thermalpx[count] = px;
						thermalpy[count] = py;
						thermalpz[count] = pz;
						thermalx[count] = xx;
						thermaly[count] = yy;
						thermalz[count] = zz;
						thermalt[count] = tt;
						count++;
					}
					if(0.5<=ram){
						energy=sqrt(px*px+py*py+pz*pz+massLambda*massLambda);
						thermalid[count] = -3122;
						thermalEng[count] = energy;
						thermalmass[count] = massLambda;

						thermalpx[count] = px;
						thermalpy[count] = py;
						thermalpz[count] = pz;
						thermalx[count] = xx;
						thermaly[count] = yy;
						thermalz[count] = zz;
						thermalt[count] = tt;
						count++;
					}
				}
				//phi
				if((abs(id)==5)||(abs(id)==55)){	
	//				double ram=random(1000)/1000.00;
						energy=sqrt(px*px+py*py+pz*pz+massphi*massphi);
						thermalid[count] = 333;
						thermalEng[count] = energy;
						thermalmass[count] = massphi;

						thermalpx[count] = px;
						thermalpy[count] = py;
						thermalpz[count] = pz;
						thermalx[count] = xx;
						thermaly[count] = yy;
						thermalz[count] = zz;
						thermalt[count] = tt;
						count++;
				}
				//Omega
				if((abs(id)==6)||(abs(id)==66)){	
					double ram=random(1000)/1000.00;
					if(ram<0.5){
						energy=sqrt(px*px+py*py+pz*pz+massOmega*massOmega);
						thermalid[count] = 3334;
						thermalEng[count] = energy;
						thermalmass[count] = massOmega;

						thermalpx[count] = px;
						thermalpy[count] = py;
						thermalpz[count] = pz;
						thermalx[count] = xx;
						thermaly[count] = yy;
						thermalz[count] = zz;
						thermalt[count] = tt;
						count++;
					}
					if(0.5<=ram){
						energy=sqrt(px*px+py*py+pz*pz+massOmega*massOmega);
						thermalid[count] = -3334;
						thermalEng[count] = energy;
						thermalmass[count] = massOmega;

						thermalpx[count] = px;
						thermalpy[count] = py;
						thermalpz[count] = pz;
						thermalx[count] = xx;
						thermaly[count] = yy;
						thermalz[count] = zz;
						thermalt[count] = tt;
						count++;
					}
				}
				// from fragment
				//pion
				if(abs(id)==111){	
						thermalmass[count] = masspi0;
						energy=sqrt(px*px+py*py+pz*pz+thermalmass[count]*thermalmass[count]);
						thermalid[count] = id;
						thermalEng[count] = energy;
						thermalpx[count] = px;
						thermalpy[count] = py;
						thermalpz[count] = pz;
						thermalx[count] = xx;
						thermaly[count] = yy;
						thermalz[count] = zz;
						thermalt[count] = tt;
						count++;
				}
				//pion
				if(abs(id)==211){	
						thermalmass[count] = masspi;
						energy=sqrt(px*px+py*py+pz*pz+thermalmass[count]*thermalmass[count]);
						thermalid[count] = id;
						thermalEng[count] = energy;
						thermalpx[count] = px;
						thermalpy[count] = py;
						thermalpz[count] = pz;
						thermalx[count] = xx;
						thermaly[count] = yy;
						thermalz[count] = zz;
						thermalt[count] = tt;
						count++;
				}
				//kaon0
				if(abs(id)==311){	
						thermalmass[count] = massK0;
						energy=sqrt(px*px+py*py+pz*pz+thermalmass[count]*thermalmass[count]);
						thermalid[count] = id;
						thermalEng[count] = energy;
						thermalpx[count] = px;
						thermalpy[count] = py;
						thermalpz[count] = pz;
						thermalx[count] = xx;
						thermaly[count] = yy;
						thermalz[count] = zz;
						thermalt[count] = tt;
						count++;
				}
				//kaon
				if(abs(id)==321){	
						thermalmass[count] = massK;
						energy=sqrt(px*px+py*py+pz*pz+thermalmass[count]*thermalmass[count]);
						thermalid[count] = id;
						thermalEng[count] = energy;
						thermalpx[count] = px;
						thermalpy[count] = py;
						thermalpz[count] = pz;
						thermalx[count] = xx;
						thermaly[count] = yy;
						thermalz[count] = zz;
						thermalt[count] = tt;
						count++;
				}
				//p
				if(abs(id)==2212){	
						thermalmass[count] = massp;
						energy=sqrt(px*px+py*py+pz*pz+thermalmass[count]*thermalmass[count]);
						thermalid[count] = id;
						thermalEng[count] = energy;
						thermalpx[count] = px;
						thermalpy[count] = py;
						thermalpz[count] = pz;
						thermalx[count] = xx;
						thermaly[count] = yy;
						thermalz[count] = zz;
						thermalt[count] = tt;
						count++;
				}
				//n
				if(abs(id)==2112){	
						thermalmass[count] = massn;
						energy=sqrt(px*px+py*py+pz*pz+thermalmass[count]*thermalmass[count]);
						thermalid[count] = id;
						thermalEng[count] = energy;
						thermalpx[count] = px;
						thermalpy[count] = py;
						thermalpz[count] = pz;
						thermalx[count] = xx;
						thermaly[count] = yy;
						thermalz[count] = zz;
						thermalt[count] = tt;
						count++;
				}
				//Lambda
				if(abs(id)==3122){	
						thermalmass[count] = massLambda;
						energy=sqrt(px*px+py*py+pz*pz+thermalmass[count]*thermalmass[count]);
						thermalid[count] = id;
						thermalEng[count] = energy;
						thermalpx[count] = px;
						thermalpy[count] = py;
						thermalpz[count] = pz;
						thermalx[count] = xx;
						thermaly[count] = yy;
						thermalz[count] = zz;
						thermalt[count] = tt;
						count++;
				}
				//phi
				if(abs(id)==333){	
						thermalmass[count] = massphi;
						energy=sqrt(px*px+py*py+pz*pz+thermalmass[count]*thermalmass[count]);
						thermalid[count] = id;
						thermalEng[count] = energy;
						thermalpx[count] = px;
						thermalpy[count] = py;
						thermalpz[count] = pz;
						thermalx[count] = xx;
						thermaly[count] = yy;
						thermalz[count] = zz;
						thermalt[count] = tt;
						count++;
				}
				//Omega
				if(abs(id)==3334){	
						thermalmass[count] = massOmega;
						energy=sqrt(px*px+py*py+pz*pz+thermalmass[count]*thermalmass[count]);
						thermalid[count] = id;
						thermalEng[count] = energy;
						thermalpx[count] = px;
						thermalpy[count] = py;
						thermalpz[count] = pz;
						thermalx[count] = xx;
						thermaly[count] = yy;
						thermalz[count] = zz;
						thermalt[count] = tt;
						count++;
				}
		}

		outfile<<ie<<"        "<<count<<"        "<<0<<"        "<<0<<endl;
		for(int ij=0;ij!=count;ij++){
			outfile<<ij<<"  "<<thermalid[ij]<<"   "<<thermalpx[ij]<<"   "<<thermalpy[ij]<<"   "<<thermalpz[ij]<<"   "<<thermalEng[ij]<<"   "<<thermalmass[ij]<<"   "<<thermalx[ij]<<"   "<<thermaly[ij]<<"   "<<thermalz[ij]<<"   "<<thermalt[ij]<<endl;
		}
	/*	
	delete thermalEng;
	delete thermalpx;
	delete thermalpy;
	delete thermalpz;
	delete thermalid;
	delete thermalmass;
	delete thermalx;
	delete thermaly;
	delete thermalz;
	delete thermalt;
	*/

	}
	outfile.close();
	return 0;
}

// check particle ID
int findarray(int*p, int len,int val)
{
	//p[len] = val;
	int i;
	int ret=0;
	for (i = 0; i!= len; i++)
	{
		if (p[i] == val)
		{ret=1;break;}
	}
	return ret;
} 
	
