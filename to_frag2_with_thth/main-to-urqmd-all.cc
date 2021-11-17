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
#define cut1 20000000.0
#define cut2 20000001.0
#define cut3 20000002.350
#define cut4 20000003.350
#define cut5 20000004.0
#define cut6 20000005.0


using namespace std;

int findarray(int*p,int len,int val);
int main( int argc, char* argv[])
{
	int nrun = atoi(argv[1]);
	//string path1=string(argv[2]);
	//string path2=string(argv[3]);

	char infilefrag[128];
	char outfiles[128];
        char outfiles1[128];
        char outfiles2[128];
        char outfiles3[128];
        char outfiles4[128];
        char outfiles5[128];
        char outfiles6[128];
        char outfiles7[128];

	int id, mid,thermalid[50000]={0},thermalcharge[50000]={0};

	int Netherml,Nejet,mark,Nfrag;
	double thermalpx[50000]={0.0},thermalpy[50000]={0.0},thermalpz[50000]={0.0},thermalEng[50000]={0.0},thermalmass[50000]={0.0},thermalx[50000]={0.0},thermaly[50000]={0.0},thermalz[50000]={0.0},thermalt[50000]={0.0},thermaleatp[50000]={0.0},etap[50000]={0.0},Qscale[50000]={0.0};
       // double px2[50000]={0.0},py2[50000]={0.0},pz2[50000]={0.0},energy2[50000]={0.0},mass2[50000]={0.0},xx2[50000]={0.0},yy2[50000]={0.0},zz2[50000]={0.0},tt2[50000]={0.0};
        double px,py,pz,energy,mass,xx,yy,zz,tt,Q0;

        double masspi=0.13957,masspi0=0.13498;
        double massK=0.49368,massK0=0.49765;
        double massp=0.93827,massn=0.93957;
        double massLambda=1.11568,massphi=1.01946,massOmega=1.67243;

        sprintf(outfiles1,"jet_parton1.dat");
	ofstream outfile1(outfiles1);
        sprintf(outfiles2,"jet_parton2.dat");
        ofstream outfile2(outfiles2);
        sprintf(outfiles3,"jet_parton3.dat");
        ofstream outfile3(outfiles3);
        sprintf(outfiles4,"jet_parton4.dat");
        ofstream outfile4(outfiles4);
        sprintf(outfiles5,"jet_parton5.dat");
        ofstream outfile5(outfiles5);
        sprintf(outfiles6,"jet_parton6.dat");
        ofstream outfile6(outfiles6);
        sprintf(outfiles7,"jet_parton7.dat");
        ofstream outfile7(outfiles7);


	
	string ipput_filename1;
	ipput_filename1 = "./remnant_jet_parton.dat";// oooooooooooooutput files of final hadrons

        sprintf(infilefrag,  ipput_filename1.c_str());
        FILE* infilef2;
        infilef2 = fopen(infilefrag,"r");
        /*
        string ipput_filename2;
        ipput_filename2 = path2+"thermal_jet.dat";// oooooooooooooutput files of final hadrons

        sprintf(infilejet, ipput_filename2.c_str());
        FILE* infile2;
        infile2 = fopen(infilejet,"r");
       
	outfile<<"OSC1997A"<<endl;
	utfile<<"final_id_p_x"<<endl;
	outfile<<" 3DHydro       1.1  (197,    79)+(197,    79)  eqsp  0.1000E+03         1"<<endl;
	*/
	int totalnum=0;
	for(int ie=0;ie!=nrun;ie++){
                fscanf(infilef2,"%d %d\n",&mid,&Nfrag);
		cout<<Nfrag<<endl;
                //getline(inhy,stempa);
                int c1=0,c2=0,c3=0,c4=0,c5=0,c6=0,c7=0;
                for(int ij=1;ij!=Nfrag+1;ij++){
                        fscanf(infilef2,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf\n", &mid,&id, &px, &py, &pz, &energy,&xx, &yy, &zz, &tt);
                        Q0=0.0;
                        if(isnan(id)||isnan(px)|| isnan(py)||isnan(pz)||isnan(energy)||isnan(xx)||isnan(yy)||isnan(zz)||isnan(tt)){
                        xx=40.0;yy=40.0; zz=40.0;tt=40.0;
                        }
			thermalid[ij]=id;
                        thermalpx[ij] = px;
                        thermalpy[ij] = py;
                        thermalpz[ij] = pz;
			thermalEng[ij]=energy;
                        thermalx[ij] = xx;
                        thermaly[ij] = yy;
                        thermalz[ij] = zz;
                        thermalt[ij] = tt;
			Qscale[ij]   = Q0;
			etap[ij]=0.5*log((energy+pz)/(energy-pz));
			if(abs(etap[ij])<=cut1)c1++;
                        if(cut1<etap[ij]&&etap[ij]<=cut2)c2++;
                        if(cut2<etap[ij]&&etap[ij]<=cut3)c3++;
                        if(cut3<etap[ij]&&etap[ij]<=cut4)c4++;
                        if(cut4<etap[ij]&&etap[ij]<=cut5)c5++;
                        if(cut5<etap[ij]&&etap[ij]<=cut6)c6++;
                        if(cut6<etap[ij])c7++;

		}

                outfile1<<mid<<"        "<<c1<<endl;
                outfile2<<mid<<"        "<<c2<<endl;
                outfile3<<mid<<"        "<<c3<<endl;
                outfile4<<mid<<"        "<<c4<<endl;
                outfile5<<mid<<"        "<<c5<<endl;
                outfile6<<mid<<"        "<<c6<<endl;
                outfile7<<mid<<"        "<<c7<<endl;

		//cout<<c1<<" "<<c2<<" "<<c3<<" "<<c4<<" "<<c5<<" "<<c1+c2+c3+c4+c5<<endl;
		for(int ij=1;ij!=Nfrag+1;ij++){
			if(abs(etap[ij])<=cut1)		
			outfile1<<mid<<"  "<<thermalid[ij]<<"   "<<thermalpx[ij]<<"   "<<thermalpy[ij]<<"   "<<thermalpz[ij]<<"   "<<thermalEng[ij]<<"   "<<"   "<<thermalx[ij]<<"   "<<thermaly[ij]<<"   "<<thermalz[ij]<<"   "<<thermalt[ij]<<" "<<Qscale[ij]<<endl;

                        if((cut1<etap[ij])&&(etap[ij]<=cut2))
                        outfile2<<mid<<"  "<<thermalid[ij]<<"   "<<thermalpx[ij]<<"   "<<thermalpy[ij]<<"   "<<thermalpz[ij]<<"   "<<thermalEng[ij]<<"   "<<"   "<<thermalx[ij]<<"   "<<thermaly[ij]<<"   "<<thermalz[ij]<<"   "<<thermalt[ij]<<endl;

                        if((cut2<etap[ij])&&(etap[ij]<=cut3))
                        outfile3<<mid<<"  "<<thermalid[ij]<<"   "<<thermalpx[ij]<<"   "<<thermalpy[ij]<<"   "<<thermalpz[ij]<<"   "<<thermalEng[ij]<<"   "<<"   "<<thermalx[ij]<<"   "<<thermaly[ij]<<"   "<<thermalz[ij]<<"   "<<thermalt[ij]<<endl;

                        if((cut3<etap[ij])&&(etap[ij]<=cut4))
                        outfile4<<mid<<"  "<<thermalid[ij]<<"   "<<thermalpx[ij]<<"   "<<thermalpy[ij]<<"   "<<thermalpz[ij]<<"   "<<thermalEng[ij]<<"   "<<"   "<<thermalx[ij]<<"   "<<thermaly[ij]<<"   "<<thermalz[ij]<<"   "<<thermalt[ij]<<endl;

                        if((cut4<etap[ij])&&(etap[ij]<=cut5))
                        outfile5<<mid<<"  "<<thermalid[ij]<<"   "<<thermalpx[ij]<<"   "<<thermalpy[ij]<<"   "<<thermalpz[ij]<<"   "<<thermalEng[ij]<<"   "<<"   "<<thermalx[ij]<<"   "<<thermaly[ij]<<"   "<<thermalz[ij]<<"   "<<thermalt[ij]<<endl;

                        if((cut5<etap[ij])&&(etap[ij]<=cut6))
                        outfile6<<mid<<"  "<<thermalid[ij]<<"   "<<thermalpx[ij]<<"   "<<thermalpy[ij]<<"   "<<thermalpz[ij]<<"   "<<thermalEng[ij]<<"   "<<"   "<<thermalx[ij]<<"   "<<thermaly[ij]<<"   "<<thermalz[ij]<<"   "<<thermalt[ij]<<endl;

                        if(cut6<etap[ij])
                        outfile7<<mid<<"  "<<thermalid[ij]<<"   "<<thermalpx[ij]<<"   "<<thermalpy[ij]<<"   "<<thermalpz[ij]<<"   "<<thermalEng[ij]<<"   "<<"   "<<thermalx[ij]<<"   "<<thermaly[ij]<<"   "<<thermalz[ij]<<"   "<<thermalt[ij]<<endl;


		}

	}
	outfile1.close();
        outfile2.close();
        outfile3.close();
        outfile4.close();
        outfile5.close();
        outfile6.close();
        outfile7.close();

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
	
