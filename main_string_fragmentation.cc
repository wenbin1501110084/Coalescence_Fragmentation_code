#include <assert.h>
#include <time.h>
#include <sstream>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string.h>
#include <algorithm>
#include <vector>
#include "Pythia8/Pythia.h" 
#define PI 3.1415926
#define QQ0 1.0

//*****************************
//  The Pythia 8 string fragmentation
//  copywrite Wenbin Zhao, 2020, 12.24
//  When you use this code, please cite our papers:
//      W.~Zhao, etc. al.[arXiv:2103.14657 [hep-ph]].                                               
//      W.~Zhao, etc. al.Phys. Rev. Lett.125, (2020) no.7, 072301.
//*****************************

using namespace Pythia8;
int searchmin(double*p, int *q,int*s,int len);
int searchmax(double*p, int *q,int len);
int searchmin2(double*p,int len);

char infiles[128];

//===========================================================================================
int main(int argv, char* argc[])
{
    string random_str = string(argc[1]);
    int n_event=atoi(argc[1]);
    string path=string(argc[2]);
    string path2=string(argc[3]);
    string output_filename2;
    output_filename2 = path+"hadrons_frag1.dat";// output files of final hadrons
    cout << output_filename2 << endl;
    string ramdomseed_str = "Random:seed = "+random_str;
    ofstream output2(output_filename2.c_str()); 
    if (!output2.is_open() ) {
        cout << "cannot open output file:"<< endl
         << output_filename2 << endl;
        return -1;
    }
    string ipput_filename2;
    ipput_filename2 = path2+"remnant_jet_parton.dat";// output files of final hadrons
    sprintf(infiles, ipput_filename2.c_str()); //input parton file 
    FILE* infile1;
    infile1 = fopen(infiles,"r");

    Pythia pythia;
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = 0");  
    pythia.readString("Beams:eCM = 5020");
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212"); 
    // Standard settings
    pythia.readString("HardQCD:all = on");
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString("ColourReconnection:reconnect = on");
    pythia.readString("ColourReconnection:mode = 1");
    pythia.readString("MultipartonInteractions:pT0Ref = 2.15");     
    pythia.readString("ColourReconnection:allowDoubleJunRem = off");
    pythia.readString("ColourReconnection:junctionCorrection = 1.15");
    pythia.readString("ColourReconnection:timeDilationMode = 3");    
    pythia.readString("ColourReconnection:timeDilationPar = 0.18");  
    pythia.readString("ProcessLevel:all = off"); // The trick!
    pythia.readString(" StringFlav:probStoUD=0.50");
    //pythia.readString("StringFlav:BtoMratio=0.5");
    pythia.readString("StringFlav:probQQtoQ=0.34");
    pythia.readString("HadronLevel:Decay = off");
    pythia.readString("HadronLevel:Hadronize = on");
    pythia.init();

    double hbarc = 0.19732;
    double c_px, c_py ,c_pz, c_e, c_m,c_x,c_y,c_z,c_t;
    int c_id,tt;
    int acol_ip[5000]={0};
    int pp_collision=0;    //used for total cross section
    int  ccbar_num=0;
    double cmeson_px,cmeson_py,cmeson_pz, cmeson_P,tempp;
    double cbar_meson_px,cbar_meson_py,cbar_meson_pz, cbar_meson_P,pt_square,cbar_meson_energy;
    double pt,mass,amid,Qmid;
    int mid,ie,II,status,col,acol,Npart,NN,Ntotal,PAosi,simble,mmaxindex;
    int idp[10000]={0},idpo[10000]={0};
    double pxpo[10000]={0.0},pypo[10000]={0.0},pzpo[10000]={0.0},epo[10000]={0.0},ptpo[10000]={0.0},xxpo[10000]={0.0},yypo[10000]={0.0},zzpo[10000]={0.0},ttpo[10000]={0.0},phio[10000]={0.0},etao[10000]={0.0};//,mass[10000]={0.0};
    double pxp[10000]={0.0},pyp[10000]={0.0},pzp[10000]={0.0},ep[10000]={0.0},ptp[10000]={0.0},xxp[10000]={0.0},yyp[10000]={0.0},zzp[10000]={0.0},ttp[10000]={0.0},phi[10000]={0.0},distance[10000]={0.0},dsting[100][1000]={0.0},Qscale[10000]={0.0};//,mass[10000]={0.0};
    int nncol[10000]={0},aacol[10000]={0},index[10000]={0},qindex[10000]={0},aqindex[10000]={0},used[10000]={0},pair[1000]={0},apair[1000]={0},gindex[1000]={0},strings[100][1000]={0},Snum[1000]={0};//,nncol_mid[10000]={0},aacol_mid[10000]={0};
    // event loop
    for (int iEvent=0; iEvent<n_event; iEvent++) {     
        fscanf(infile1,"%d %d\n",&mid,&Npart);
        if (Npart==0 ) {output2 << iEvent+1<<" "<<0 << endl;continue;} 
        int Nquark=0;
        int Naquark =0;
        int Ngluon=0;
        int Npair=0;
        for (int ll=0;ll<Npart;ll++) {
                fscanf(infile1,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf\n",&mid,&c_id,&c_px,&c_py,&c_pz,&c_e,&c_x, &c_y, &c_z,&c_t);// format of input partons
                Qmid = 0.0; // The scale for the parton shower. 
                idpo[ll]=c_id;
                if(c_id==21){gindex[Ngluon]=ll;Ngluon++;}
                pxpo[ll]=c_px;
                pypo[ll]=c_py;
                pzpo[ll]=c_pz;
                ptpo[ll]=c_px*c_px+c_py*c_py;
                double pmg=sqrt(c_px*c_px+c_py*c_py+c_pz*c_pz);
                epo[ll]=c_e;
                xxpo[ll]=c_x;
                yypo[ll]=c_y;
                zzpo[ll]=c_z;
                ttpo[ll]=c_t;
                Qscale[ll]=sqrt(Qmid);
                double aamid=0.5*log((pmg+c_pz)/(pmg-c_pz));
                etao[ll]=aamid;
                phio[ll]=atan2(c_py,c_px);
                index[ll]=ll;
                used[ll]=0;
        }
        mmaxindex=searchmax(ptpo,index,Npart);// get the leading parton
        //**** all partons are gluon ****
        // add fictive quark anti-quark
        if (Npart>0&&Ngluon==Npart) {
            idpo[Npart]=-3;pxpo[Npart]=0.10;pypo[Npart]=0.20;pzpo[Npart]=100000.0;
            epo[Npart]=sqrt(pzpo[Npart]*pzpo[Npart]+pypo[Npart]*pypo[Npart]+pxpo[Npart]*pxpo[Npart]+0.5*0.5);
            xxpo[Npart]=0.0;yypo[Npart]=0.0;zzpo[Npart]=0.0;ttpo[Npart]=0.0;phio[Npart]=atan2(pypo[Npart],pxpo[Npart]);index[Npart]=Npart;
            used[Npart]=0;
            etao[Npart]=10000.0;
            Npart=Npart+1;
            idpo[Npart]=3;pxpo[Npart]=0.10;pypo[Npart]=0.20;pzpo[Npart]=-100000.0;
            epo[Npart]=sqrt(pzpo[Npart]*pzpo[Npart]+pypo[Npart]*pypo[Npart]+pxpo[Npart]*pxpo[Npart]+0.5*0.5);
            xxpo[Npart]=0.0;yypo[Npart]=0.0;zzpo[Npart]=0.0;ttpo[Npart]=0.0;phio[Npart]=atan2(pypo[Npart],pxpo[Npart]);index[Npart]=Npart;
            used[Npart]=0;
            etao[Npart]=-10000.0;
            Npart=Npart+1;
        }

        // ****************** give the col and acol to gluon and quarks *********************//
        int Nmax=101;
        int add=0;
        int jj,Nmin;
        //*** find the quark pairs ***
        for (int ll=0;ll<Npart;ll++) {
            if ((idpo[ll]<=5)&((idpo[ll])>0)) {
                qindex[Nquark]=ll;Nquark++;
            }
            if ((idpo[ll]>=-5)&((idpo[ll])<0)) {
                aqindex[Naquark]=ll;Naquark++;
            }
        }
        int excess=abs(Nquark-Naquark);
        if (excess>0) {
            int ppp=excess % 2;
            int ecc=excess/2;   
            if (Nquark>Naquark) {
                if (ecc>0) {
                    for (int gg=0; gg<ecc; gg++) {
                        idpo[qindex[Nquark-1]]=-1*idpo[qindex[Nquark-1]];
                        aqindex[Naquark]=qindex[Nquark-1];
                        Naquark++;
                        Nquark--;
                    }
                }
                if (ppp==1) {
                    idpo[Npart]=-3; pxpo[Npart]=0.20; pypo[Npart]=0.20;
                    pzpo[Npart]=10000.20;epo[Npart]=sqrt(pzpo[Npart]*pzpo[Npart]+pypo[Npart]*pypo[Npart]+pxpo[Npart]*pxpo[Npart]+0.5*0.5);
                    xxpo[Npart]=0.0;yypo[Npart]=0.0;zzpo[Npart]=0.0;ttpo[Npart]=0.0;phio[Npart]=atan2(pypo[Npart],pxpo[Npart]);index[Npart]=Npart; 
                    etao[Npart]=1.0;
                    used[Npart]=0;
                    aqindex[Naquark]=Npart;
                    Npart++;
                    Naquark++;
                }
            }
            if (Naquark>Nquark) {
                if (ecc>0) {
                    for (int gg=0; gg<ecc; gg++) {
                        idpo[aqindex[Naquark-1]]=-1*idpo[aqindex[Naquark-1]];
                        qindex[Nquark]=aqindex[Naquark-1];
                        Nquark++;
                        Naquark--;
                    }
                }
                if (ppp==1) {
                    idpo[Npart]=3; pxpo[Npart]=0.20; pypo[Npart]=0.20;
                    pzpo[Npart]=100000.20;epo[Npart]=sqrt(pzpo[Npart]*pzpo[Npart]+pypo[Npart]*pypo[Npart]+pxpo[Npart]*pxpo[Npart]+0.5*0.5);
                    xxpo[Npart]=0.0;yypo[Npart]=0.0;zzpo[Npart]=0.0;ttpo[Npart]=0.0;phio[Npart]=atan2(pypo[Npart],pxpo[Npart]);index[Npart]=Npart;
                    etao[Npart]=1.0;
                    used[Npart]=0;
                    qindex[Nquark]=Npart;
                    Npart++;
                    Nquark++;
                }
            }
        }
        // ************* find the quark-anti-quark pairs with smallest distance *************
       //Nquark>Naquark>0
        if (Nquark>Naquark) {
            // if the only one (anti-)quark
            if(Naquark==1){
                int used_q[1000]={0};
                for (int pp=0;pp<Nquark;pp++) {
                    distance[pp]=(phio[aqindex[0]]-phio[qindex[pp]])*(phio[aqindex[0]]-phio[qindex[pp]])+(etao[aqindex[0]]-etao[qindex[pp]])*(etao[aqindex[0]]-etao[qindex[pp]]);
                    used_q[pp]=used[qindex[pp]];
                }
                int minindex=searchmin(distance,qindex,used_q,Nquark);
                apair[Npair]=aqindex[0]; used[aqindex[0]]=1; pair[Npair]=minindex; used[minindex]=1; Npair++;// first pair
                for(int pp=0;pp<Nquark;pp++){
                    int iindex=qindex[pp];
                    if(used[iindex]==1)continue;
                    idpo[Npart]=-3;pxpo[Npart]=0.0;pypo[Npart]=0.0;pzpo[Npart]=100000.0;epo[Npart]=sqrt(pzpo[Npart]*pzpo[Npart]+pypo[Npart]*pypo[Npart]+pxpo[Npart]*pxpo[Npart]+0.5*0.5);
                    apair[Npair]=Npart; used[Npart]=1; pair[Npair]=iindex; used[iindex]=1; Npair++;//  pair
                    Npart++;
                }
            }else{
                //sort the input partons according to pT
                int temindex[1000]={0};
                for(int yy=0;yy<Naquark;yy++){
                     temindex[yy]=aqindex[yy];
                 }
                for (int i = 0; i < Naquark-1; i++){
                    for (int j = 0; j < Naquark - 1 - i; j++){
                        if (ptpo[temindex[j]] > ptpo[temindex[j + 1]]) {
                            tempp = ptpo[temindex[j]];
                            ptpo[temindex[j]] = ptpo[temindex[j + 1]];
                            ptpo[temindex[j+1]] = tempp;
                            tt=aqindex[j];
                            aqindex[j] = aqindex[j+1];
                            aqindex[j+1] = tt; 
                        }
                    }
                }
                // Assign the pairs
                for (int ii=0; ii<Naquark; ii++) {
                    int used_q[1000]={0};
                    for (int pp=0;pp<Nquark;pp++) {
                        distance[pp]=(phio[aqindex[ii]]-phio[qindex[pp]])*(phio[aqindex[ii]]-phio[qindex[pp]])+(etao[aqindex[ii]]-etao[qindex[pp]])*(etao[aqindex[ii]]-etao[qindex[pp]]);
                        used_q[pp]=used[qindex[pp]];
                    }
                    int minindex=searchmin(distance,qindex,used_q,Nquark);
                    apair[Npair]=aqindex[ii]; used[aqindex[ii]]=1; pair[Npair]=minindex; used[minindex]=1; Npair++;//  pair
                }
            }
            // the remnant quarks
            for(int yy=0;yy<Nquark;yy++){
                if(used[qindex[yy]]==1)continue;
                idpo[Npart]=-3;pxpo[Npart]=0.0;pypo[Npart]=0.0;pzpo[Npart]=100000.0;
                epo[Npart]=sqrt(pzpo[Npart]*pzpo[Npart]+pypo[Npart]*pypo[Npart]+pxpo[Npart]*pxpo[Npart]+0.5*0.5);
                xxpo[Npart]=0.0;yypo[Npart]=0.0;zzpo[Npart]=0.0;ttpo[Npart]=0.0;phio[Npart]=atan2(pypo[Npart],pxpo[Npart]);index[Npart]=Npart;
                apair[Npair]=Npart; used[Npart]=1; pair[Npair]=qindex[yy]; used[qindex[yy]]=1; Npair++;//  pair
                Npart++;
            }
        }
    
        //Naquark>Nquark>0
        if(Naquark>Nquark){
            // if the only one (anti-)quark
            if(Nquark==1){
                int used_q[1000]={0};
                for(int pp=0;pp<Naquark;pp++){
                    distance[pp]=(phio[qindex[0]]-phio[aqindex[pp]])*(phio[qindex[0]]-phio[aqindex[pp]])+(etao[qindex[0]]-etao[aqindex[pp]])*(etao[qindex[0]]-etao[aqindex[pp]]);
                    used_q[pp]=used[aqindex[pp]];
                }
                int minindex=searchmin(distance,aqindex,used_q,Naquark);
                pair[Npair]=qindex[0]; used[qindex[0]]=1; apair[Npair]=minindex; used[minindex]=1; Npair++;// first pair
                for(int pp=0;pp<Naquark;pp++){
                    int iindex=aqindex[pp];
                    if(used[iindex]==1)continue;
                    idpo[Npart]=3;pxpo[Npart]=0.0;pypo[Npart]=0.0;pzpo[Npart]=100000.0;epo[Npart]=sqrt(pzpo[Npart]*pzpo[Npart]+pypo[Npart]*pypo[Npart]+pxpo[Npart]*pxpo[Npart]+0.5*0.5);
                    pair[Npair]=Npart; used[Npart]=1; apair[Npair]=iindex; used[iindex]=1; Npair++;//  pair
                    Npart++;
                }
            }else{
                //sort the input partons according to pT// track the index
                int temindex[1000]={0};
                for(int yy=0;yy<Nquark;yy++){
                     temindex[yy]=qindex[yy];
                 }
                for (int i = 0; i < Nquark-1; i++){
                    for (int j = 0; j < Nquark - 1 - i; j++){
                        if (ptpo[temindex[j]] > ptpo[temindex[j + 1]]) {
                            tempp = ptpo[temindex[j]];
                            ptpo[temindex[j]] = ptpo[temindex[j + 1]];
                            ptpo[temindex[j+1]] = tempp;
                            tt=qindex[j];
                            qindex[j] = qindex[j+1];
                            qindex[j+1] = tt; 
                        }
                    }
                }
                // Assign the pairs
                for(int ii=0;ii<Nquark;ii++){
                    int used_q[1000]={0};
                    for(int pp=0;pp<Naquark;pp++){
                        distance[pp]=(phio[qindex[ii]]-phio[aqindex[pp]])*(phio[qindex[ii]]-phio[aqindex[pp]])+(etao[qindex[ii]]-etao[aqindex[pp]])*(etao[qindex[ii]]-etao[aqindex[pp]]);
                        used_q[pp]=used[aqindex[pp]];
                    }
                    int minindex=searchmin(distance,aqindex,used_q,Nquark);
                    pair[Npair]=qindex[ii]; used[qindex[ii]]=1; apair[Npair]=minindex; used[minindex]=1; Npair++;//  pair
                }
            }
            // the remnant ani-quarks
            for(int yy=0;yy<Naquark;yy++){
                if(used[aqindex[yy]]==1)continue;
                idpo[Npart]=3;pxpo[Npart]=0.0;pypo[Npart]=0.0;pzpo[Npart]=100000.0;
                epo[Npart]=sqrt(pzpo[Npart]*pzpo[Npart]+pypo[Npart]*pypo[Npart]+pxpo[Npart]*pxpo[Npart]+0.5*0.5);
                xxpo[Npart]=0.0;yypo[Npart]=0.0;zzpo[Npart]=0.0;ttpo[Npart]=0.0;phio[Npart]=atan2(pypo[Npart],pxpo[Npart]);index[Npart]=Npart;
                pair[Npair]=Npart; used[Npart]=1; apair[Npair]=aqindex[yy]; used[aqindex[yy]]=1; Npair++;//  pair
                Npart++;
            }
        }
    
        //Naquark=Nquark>0
        if(Naquark==Nquark){
            if(Nquark>1){cout<<"df";
                //sort the input partons according to Delta R
                for(int www=0;www<Nquark;www++){
                    double disquark[200][200]={0.0};
                    // calculate the distance between quark pairs
                    for(int ii=0;ii<Nquark;ii++){
                        for(int pp=0;pp<Naquark;pp++){
                            if((used[qindex[ii]]==0)&&(used[aqindex[pp]]==0)){
                                double sphi=phio[qindex[ii]]-phio[aqindex[pp]];
                                if(abs(sphi)>PI){sphi=2*PI-abs(sphi);}
                                disquark[ii][pp]=sphi*sphi+(etao[qindex[ii]]-etao[aqindex[pp]])*(etao[qindex[ii]]-etao[aqindex[pp]]);
                            }else{
                                disquark[ii][pp]=1000000000000000.0;
                            }
                        }
                    }
                    // get the index of the minum value of distance
                    int qnm=0,aqnm=0;
                    double mm=disquark[0][0];
                    for (int i = 0; i < Nquark; i++){
                        for(int j=0;j<Naquark;j++){
                            if (mm > disquark[i][j]){
                                mm = disquark[i][j];
                                qnm = i;
                                aqnm=j;
                            }
                        }
                    }
                    pair[Npair]=qindex[qnm]; used[qindex[qnm]]=1; apair[Npair]=aqindex[aqnm]; used[aqindex[aqnm]]=1; Npair++;//  pair
                }
            }
            if(Naquark==1){
                pair[0]=qindex[0];used[qindex[0]]=1; apair[0]=aqindex[0]; used[aqindex[0]]=1; Npair++;//  pair
            }
        }
        // *** connect gluon the quark anti-quark pairs
        for(int zz=0;zz<Npair;zz++){
            Snum[zz]=0;
        }
        if(Ngluon>0){
            for(int pp=0;pp<Ngluon;pp++){
                double disg[1000]={0.0},disq[1000]={0.0};
                for(int nn=0;nn<Npair;nn++){
                    double sphi=phio[gindex[pp]]-phio[pair[nn]];
                    if(abs(sphi)>PI){sphi=2*PI-abs(sphi);}
                    double d1=sphi*sphi+(etao[gindex[pp]]-etao[pair[nn]])*(etao[gindex[pp]]-etao[pair[nn]]);
                    double sphi2=phio[gindex[pp]]-phio[apair[nn]];
                    if(abs(sphi2)>PI){sphi2=2*PI-abs(sphi2);}
                    double d2=sphi2*sphi2+(etao[gindex[pp]]-etao[apair[nn]])*(etao[gindex[pp]]-etao[apair[nn]]);
                    disg[nn]=sqrt(d1)+sqrt(d2); disq[nn]=d1;
                }
                int minindex=searchmin2(disg,Npair);//minindex belongs to Npair
                    strings[minindex][Snum[minindex]]=gindex[pp];
                    dsting[minindex][Snum[minindex]]=disq[minindex];
                    Snum[minindex]++;
            }
        }
        // sort the gluons // track the index
        for(int kkk=0;kkk<Npair;kkk++){
            if(Snum[kkk]>=2){
                int temindex[1000]={0};
                     for(int yy=0;yy<Snum[kkk];yy++){
                         temindex[yy]=strings[kkk][yy];
                     }
                for (int i = 0; i < Snum[kkk]-1; i++){
                    for (int j = 0; j < Snum[kkk] - 1 - i; j++){
                        if (dsting[kkk][j] > dsting[kkk][j+1]) {
                            tempp = dsting[kkk][j];
                            dsting[kkk][j] = dsting[kkk][j+1];
                            dsting[kkk][j+1] = tempp;
                            tt=strings[kkk][j];
                            strings[kkk][j] = strings[kkk][j+1];
                            strings[kkk][j+1] = tt; 
                        }
                    }
                }
            }
        }
        //***** connect the partons by strings *****
        int Cocon=102; 
        cout<<"=== The string connections ==="<<endl;
        for(int pp=0;pp<Npair;pp++){
            int Aocon=Cocon+Snum[pp];
            int Cocon2=Cocon;
            cout<<"====   "<< Snum[pp] <<"   ======"<<endl;
            nncol[pair[pp]]=Cocon;  aacol[pair[pp]]=0;// quark
            cout<< nncol[pair[pp]]<<" "<<aacol[pair[pp]]<<" "<<idpo[pair[pp]]<<endl;
            for(int yy=0;yy<Snum[pp];yy++){
                int I1=strings[pp][yy]; 
                aacol[I1]=Cocon; Cocon++;
                nncol[I1]=Cocon; cout<< nncol[I1]<<" "<<idpo[I1]<<" ";cout<< aacol[I1]<<" "<<idpo[I1]<<endl;
            }
            nncol[apair[pp]]=0; aacol[apair[pp]]=Cocon;Cocon++;// anti-quark
            cout<< nncol[apair[pp]]<<" "<<aacol[apair[pp]]<<" "<<idpo[apair[pp]]<<endl;
        }
        cout<<"================================"<<endl;

// ****************** append the partons into pythia event *********************//
        double m_str=0.0, x_str=0.0,y_str=0.0,z_str=0.0,t_str=0.0;
        double x_hadron,y_hadron,z_hadron,t_hadron,hmt;
        pythia.event.clear();
        double maxQ0 = 2.0;//maxQ0>=QQ0; wenbin 
        double minQ0 = 0.4;//minum Q0;
        pythia.event.reset();
        for (int tt=0;tt<Npart;tt++){
            if(idpo[tt]==21){mass=0.0;}
            else{
                //if(abs(idp[tt])<=2)mass=0.330;
                //if(abs(idp[tt])==3)mass=0.50;
                mass=sqrt(abs(epo[tt]*epo[tt]-pzpo[tt]*pzpo[tt]-pypo[tt]*pypo[tt]-pxpo[tt]*pxpo[tt]));
            }
            pythia.event.append(idpo[tt],62,nncol[tt],aacol[tt],pxpo[tt],pypo[tt],pzpo[tt],epo[tt],mass);
            if(maxQ0<Qscale[tt])Qscale[tt]=maxQ0;
            //if(minQ0>Qscale[tt])Qscale[tt]=minQ0;
            //pythia.event[tt].scale(Qscale[tt]);//QQ0 the initial scale of input partons 
            // get the center of mass of the strings and corresponding posistion 
            m_str=m_str+mass;
            x_str=x_str+mass*xxpo[tt];
            y_str=y_str+mass*yypo[tt];
            z_str=z_str+mass*zzpo[tt];
            t_str=t_str+mass*ttpo[tt];
        }
        if(m_str==0)m_str=0.10;
//****** fragment the remnant partons ***********
        //pythia.forceTimeShower(1,Npart,maxQ0);//Continue the FSR to the defaulted scale 
        pythia.forceHadronLevel();
        int simble = 0;
        for(int i=0; i<pythia.event.size();i++) {
            if (pythia.event[i].isFinal() ) {
                c_id = pythia.event[i].id();
                if( (abs(c_id) !=22)&&(abs(c_id) !=11)) {
                    simble=simble+1;
                }
            }
        }
        if(simble==0){output2 << iEvent+1<<" "<<simble << endl;}
        if(simble>0){
            output2 << iEvent+1<<" "<<simble << endl;
            for(int i=0; i<pythia.event.size();i++)
                {
                if (pythia.event[i].isFinal() ){
                    c_id = pythia.event[i].id();
                    if ( (abs(c_id) !=22)&&(abs(c_id) !=11)) {
                        cbar_meson_px = pythia.event[i].px();
                        cbar_meson_py = pythia.event[i].py();
                        cbar_meson_pz = pythia.event[i].pz();
                        cbar_meson_energy = pythia.event[i].e();
                        // get the posistion of the final hadrons
                        hmt=sqrt(pythia.event[i].m()*pythia.event[i].m()+cbar_meson_px*cbar_meson_px+cbar_meson_py*cbar_meson_py);
                        x_hadron=x_str/m_str+hbarc*cbar_meson_px/hmt;
                        y_hadron=y_str/m_str+hbarc*cbar_meson_py/hmt;
                        z_hadron=z_str/m_str+hbarc*cbar_meson_pz/hmt;
                        t_hadron=t_str/m_str+hbarc*cbar_meson_energy/hmt;
                        output2 << c_id << " "<<cbar_meson_px<<"  "<<cbar_meson_py<<"  "<<cbar_meson_pz <<" "<<cbar_meson_energy << "  "<< pythia.event[i].m()<<" "<< x_hadron<<" "<<y_hadron<<" "<<z_hadron<<" "<<t_hadron<<" "<<endl;
                    }
                }
            }
        }
        pythia.next();

}
output2.close();


  return 0;
}

int searchmin(double*p,int*q,int*s,int len)
{
    double m = 100000000.0;
    int k;
    for (int i = 0; i < len; ++i)
    {
        if ((m > p[i])&&(s[i]==0))
        {
            m = p[i];
            k = i;
        }
    }
    return q[k];
} 


int searchmax(double*p,int*q,int len)
{
    double m = 0.0;
    int k;
    for (int i = 0; i < len; ++i)
    {
        if (m < p[i])
        {
            m = p[i];
            k = i;
        }
    }
    return q[k];
} 

    
int searchmin2(double*p,int len)
{
    double m = 100000000.0;
    int k;
    for (int i = 0; i < len; ++i)
    {
        if (m > p[i])
        {
            m = p[i];
            k = i;
        }
    }
    return k;
} 
    
