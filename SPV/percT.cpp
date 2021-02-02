#include "std_include.h"

#include "Simulation.h"

#include "voronoiQuadraticEnergy.h"
#include "selfPropelledParticleDynamics.h"
#include "DatabaseNetCDFSPV.h"

#include <cstring>

using namespace std;

int main(){

    int ini_size=21;
    double p0_ini[]={3.75, 3.76, 3.77, 3.78, 3.79, 3.80, 3.81, 3.82, 3.83, 3.84, 3.85, 3.86, 3.87, 3.88, 3.89, 3.90, 3.91, 3.92, 3.93, 3.94, 3.95};
    double tension_avg_ini[]={0.0243023,0.0195984,0.016296,0.0140956,0.0125365,0.012417,0.0126778,0.0120944,0.0103144,0.00970038,0.00877825,0.00717078,0.00699281,0.00611137 ,0.00498653 ,0.00391426,0.00351072,0.00265684,0.00232919,0.00213616,0.00117537};
    double tension_sigma_ini[]={0.0204401,0.0208826,0.0209951,0.0206429,0.0214414,0.0212116,0.0218593,0.0213172,0.0219175,0.0218288,0.0219368,0.0222094,0.0222121,0.022528,0.0220646,0.0226943,0.0224833,0.022345,0.0224517,0.0219839,0.0222303};

    //double tension_avg_ini[]={0.0485758, 0.0399299, 0.0328189, 0.0284635, 0.0253197, 0.02551, 0.0257305, 0.024717, 0.0208015, 0.0193719, 0.0175794, 0.0143416, 0.0136164, 0.0121144, 0.00944348, 0.00769726, 0.00650248, 0.00530765, 0.00440967, 0.00354915, 0.00184046};
    //double tension_sigma_ini[]={0.027757, 0.0276927, 0.0280816, 0.0279209, 0.0290757, 0.0288898, 0.0300414, 0.0296319, 0.02998, 0.0301359, 0.0308032, 0.0306167, 0.0310135, 0.0310541, 0.0306512, 0.0311785, 0.0310794, 0.0312943, 0.0313861, 0.0316016, 0.0312311};

    int cont=0;
    int last=60;
    int scriptsIni=16;
    for(int start=1; start<=last; start++){

      int v1, v2, v3, v4, v10, v11;
      double v5, v6, v7, v8, v9;
      //string access_dados="/home/diogo/MEGA/cenas/CFTC_Cluster/scripts35/dados.txt";
      //string access_dados="/media/hdd/Results/scripts39/dados.txt";
      //string access_dados="/media/hdd/Results/GPU/scripts5/dados.txt";
      //string access_dados="/home/diogo/MEGA/cenas/GPU/Results_glass/scripts"+to_string(scriptsIni)+"/dados.txt";
      //string access_dados="/media/hdd/Results_GRF/scripts"+to_string(scriptsIni)+"/dados.txt";
      string access_dados="/mnt/d/Results/GRF/scripts"+to_string(scriptsIni)+"/dados.txt";
      ifstream file(access_dados);

      for(int line=0; line<start; line++){
	file>> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7 >> v8 >> v9 >> v10 >> v11;
      }
      file.close();

      char data[1000000]; 
      //sprintf(dataname, "/home/diogo/MEGA/cenas/CFTC_Cluster/scripts1/Job_%d/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", avg, v9, v10, v11, v6, v7, v8);
      //sprintf(data, "/media/hdd/Results/scripts39/Job_%d/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", (start-1)*v14+1, v9, v10, v11, v6, v7, v8);
      //sprintf(data, "/media/hdd/Results/GPU/scripts5/Job_%d/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", (start-1)*v14+1, v9, v10, v11, v6, v7, v8);
      //sprintf(data, "/home/diogo/MEGA/cenas/GPU/Results_glass/scripts%d/Job_%d/test_voronoi.nc", scriptsIni, (start-1)*v11+1);
      //sprintf(data, "/media/hdd/Results_GRF/scripts%d/Job_%d/test_voronoi.nc", scriptsIni, (start-1)*v11+1);
      sprintf(data, "/mnt/d/Results/GRF/scripts%d/Job_%d/test_voronoi.nc", scriptsIni, (start-1)*v11+1);
 
      SPVDatabaseNetCDF nc(v1,data,NcFile::ReadOnly);

        int N=v1;
        int a=N;
	double L=sqrt(N);
        int b=nc.GetNumRecs();
        int total_size=v2;
        int amostras=v11;
        //int start=1;
        int fim=start;
        vector <double> time(b,0);
	vector <double> cluster_size_avg_p(b, 0);
        vector <double> perc_avg_p(b, 0);
	vector <double> phi_avg_p(b, 0);
	vector <double> phi_avg_p_2(b, 0);
	int numBins=100;
	double cdistSize=(double(N+1))/double(numBins);
	vector <double> cluster_dist_p(numBins, 0);
        //double max=v6;
	//double min=max-v11;
        string filename;
        double dt=v5;
	int LA=last*amostras;
	//double pThresh=3.8;
	double pThresh=0.08;

	int idx=0;
	for(int i=0; i<ini_size; i++)
	{
		if(abs(p0_ini[i]-v6)<0.00001)idx=i;
	}
	pThresh=3*tension_sigma_ini[idx];

        vector <int> cluster_p(N,0);
        vector <int> cluster_list_p(N,-1);
        vector <int> cluster_perc_p(N,0);
        vector <int> cluster_size_p(N,0);

        shared_ptr<VoronoiQuadraticEnergy> spv  = make_shared<VoronoiQuadraticEnergy>(N,1.0,4.0,true);
        FILE *f;
 
        for(int avg=(start-1)*amostras+1; avg<=fim*amostras; avg++){
            cont=cont+1;
            int curr_line=(avg-1)/amostras+1;           
	    file.open(access_dados);
            for(int line=0; line<curr_line; line++){
                file>> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7 >> v8 >> v9 >> v10 >> v11;
            }
            file.close();
            string add_JOB_num=to_string(avg);
            cout<<add_JOB_num<<endl;
           
            char dataname[1000000]; 
            //sprintf(dataname, "/home/diogo/MEGA/cenas/CFTC_Cluster/scripts1/Job_%d/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", avg, v9, v10, v11, v6, v7, v8);
	    //sprintf(dataname, "/media/hdd/Results/scripts39/Job_%d/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", avg, v9, v10, v11, v6, v7, v8);
	    //sprintf(dataname, "/media/hdd/Results/GPU/scripts5/Job_%d/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", avg, v9, v10, v11, v6, v7, v8);
            //sprintf(dataname, "/home/diogo/MEGA/cenas/GPU/Results_glass/scripts%d/Job_%d/test_voronoi.nc", scriptsIni, avg);
            //sprintf(dataname, "/media/hdd/Results_GRF/scripts%d/Job_%d/test_voronoi.nc", scriptsIni, avg);
	    sprintf(dataname, "/mnt/d/Results/GRF/scripts%d/Job_%d/test_voronoi.nc", scriptsIni, avg);
	    
            SPVDatabaseNetCDF ncdat(N,dataname,NcFile::ReadOnly);
            //Check if the file exists in the output folder. if it does then do the scan
            if ((f = fopen(dataname, "r")) == NULL){
                printf("Error! opening file\n");
                printf("%s\n",dataname);
                exit (911);
                return -1;
            }
            else{
                fclose(f);
                ncdat.ReadState(spv,0,false);
                b=ncdat.GetNumRecs();
                //printf("b=%u\n", b);
            }
           
            //ArrayHandle<double2> h_p(spv->cellPositions,access_location::host,access_mode::read);
            //ArrayHandle<double2> h_app(spv->AreaPeriPreferences,access_location::host,access_mode::read);
            //ArrayHandle<int> h_nn(spv->neighborNum,access_location::host,access_mode::read);
            //ArrayHandle<int> h_n(spv->neighbors,access_location::host,access_mode::read);
            //ArrayHandle<double2> h_ap(spv->returnAreaPeri(),access_location::host,access_mode::read);

            double x11,x12,x21,x22;
            spv->Box->getBoxDims(x11,x12,x21,x22);
            L=x11;
	    
            for(int nn = 0; nn<b; nn++){

	        ncdat.ReadState(spv,nn,true);
		time[nn]= spv->currentTime;
                ArrayHandle<double2> h_p(spv->cellPositions,access_location::host,access_mode::read);
                ArrayHandle<double2> h_app(spv->AreaPeriPreferences,access_location::host,access_mode::read);
                ArrayHandle<int> h_nn(spv->neighborNum,access_location::host,access_mode::read);
                ArrayHandle<int> h_n(spv->neighbors,access_location::host,access_mode::read);
                ArrayHandle<double2> h_ap(spv->returnAreaPeri(),access_location::host,access_mode::read);

		fill(cluster_p.begin(), cluster_p.end(), 0);
		fill(cluster_list_p.begin(), cluster_list_p.end(), -1);
		fill(cluster_perc_p.begin(), cluster_perc_p.end(), 0);
		fill(cluster_size_p.begin(), cluster_size_p.end(), 0);

		int cluster_val_p=0;
		int cluster_array_size_p=0;
		double maxp=-10000;
		int CValp=0;
		vector <int> cleft(0, 0);
		vector <int> ctop(0, 0);
                vector <int> cright(0, 0);
                vector <int> cbottom(0, 0);
		vector <int> cBoundaries(N, 0);
		for(int i=0; i<N; i++)
		{
			cluster_array_size_p=0;
			if(cBoundaries[i]==0)
			{
				for(int k=0; k<h_nn.data[i]; k++)
				{
					int nsite=h_n.data[spv->n_idx(k, i)];
					double x=sqrt((h_p.data[i].x-h_p.data[nsite].x)*(h_p.data[i].x-h_p.data[nsite].x));
					double y=sqrt((h_p.data[i].y-h_p.data[nsite].y)*(h_p.data[i].y-h_p.data[nsite].y));
					if(x>L/2)
					{
						if(h_p.data[i].x>h_p.data[nsite].x && cBoundaries[i]!=1 && cBoundaries[i]<5)
						{
							cright.push_back(i);
							if(cBoundaries[i]==0)cBoundaries[i]=1;
							else
                                                        {
                                                                if(cBoundaries[i]==2)cBoundaries[i]=5;
                                                                else cBoundaries[i]=8;
                                                        }

						}
						else if (h_p.data[i].x<h_p.data[nsite].x && cBoundaries[i]!=3 && cBoundaries[i]<5)
						{
							cleft.push_back(i);
							if(cBoundaries[i]==0)cBoundaries[i]=3;
							else
                                                        {
                                                                if(cBoundaries[i]==2)cBoundaries[i]=6;
                                                                else cBoundaries[i]=7;
                                                        }
						}
					}
					if(y>L/2)
					{
                                                if(h_p.data[i].y>h_p.data[nsite].y && cBoundaries[i]!=2 && cBoundaries[i]<5)
                                                {
                                                        ctop.push_back(i);
							if(cBoundaries[i]==0)cBoundaries[i]=2;
							else
							{
								if(cBoundaries[i]==1)cBoundaries[i]=5;
								else cBoundaries[i]=6;
							}
                                                }
                                                else if(h_p.data[i].y<h_p.data[nsite].y && cBoundaries[i]!=4 && cBoundaries[i]<5)
                                                {
                                                        cbottom.push_back(i);
							if(cBoundaries[i]==0)cBoundaries[i]=4;
                                                        else 
							{
								if(cBoundaries[i]==1)cBoundaries[i]=8;
								else cBoundaries[i]=7;
							}
                                                }
					}
				}
			}

			//percolation of tensions cluster
			//if(h_app.data[i].y<pThresh && cluster_list_p[i]<0)
			double tensionVal=h_ap.data[i].y-h_app.data[i].y;
                        if(abs(tension_avg_ini[idx]-tensionVal)>pThresh && cluster_list_p[i]<0)
                        {
                                cluster_val_p++;
				cluster_array_size_p++;
                                cluster_p[cluster_array_size_p-1]=i;
                                cluster_perc_p[cluster_val_p-1]=0;
                                cluster_size_p[cluster_val_p-1]=0;
				cluster_list_p[i]=cluster_val_p;

                                while(cluster_array_size_p>0)
                                {
                                        int pos=cluster_p[cluster_array_size_p-1];
					cluster_array_size_p--;
                                        cluster_size_p[cluster_val_p-1]++;
                                        int size=h_nn.data[pos];
                                        int neigh=h_n.data[spv->n_idx(0,pos)];
                                        for(int j=0; j<size; j++)
                                        {
                                                int next=h_n.data[spv->n_idx((j+1)%size,pos)];
						//if(h_app.data[neigh].y<pThresh && cluster_list_p[neigh]<0)
                                                if(abs(tension_avg_ini[idx]-(h_ap.data[neigh].y-h_app.data[neigh].y))>pThresh && cluster_list_p[neigh]<0)
						{
							cluster_array_size_p++;
							cluster_p[cluster_array_size_p-1]=neigh;
							cluster_list_p[neigh]=cluster_val_p-1;
						}
                                                neigh=next;
                                        }
                                }
				if(cluster_size_p[cluster_val_p-1]>maxp)
                                {
                                        maxp=cluster_size_p[cluster_val_p-1];
                                        CValp=cluster_val_p-1;
                                }
                        }
		}

		if(maxp>0)
		{
			phi_avg_p[nn]+=maxp/N;
			phi_avg_p_2[nn]+=(maxp/N)*(maxp/N);
		}
		else
		{
			phi_avg_p[nn]+=0;
                        phi_avg_p_2[nn]+=0;
		}
		
		//check if largest cluster percolates horizontally
		bool percCluster=false;
		cluster_array_size_p=0;
		fill(cluster_p.begin(), cluster_p.end(), 0);
		vector <double> disp(N, 0);
		double maxd=-100000;
		int maxSite=0;
		for(int i=0; i<cleft.size(); i++)
		{
			int pos=cleft[i];
			bool rightC=false;
			if(cluster_list_p[pos]==CValp)
			{
				for(int k=0; k<cright.size(); k++)
				{
					int ps=cright[k];
					if(cluster_list_p[ps]==CValp)
					{
						rightC=true;
						break;
					}
				}

				if(rightC==true)
				{
					cluster_array_size_p=1;
					cluster_p[cluster_array_size_p-1]=pos;
					vector <int> newlist(N, 0);
					newlist[pos]=1;
					while(cluster_array_size_p>0 && percCluster==false)
					{
						int site=cluster_p[cluster_array_size_p-1];
						cluster_array_size_p--;
						for(int l=0; l<h_nn.data[site]; l++)
						{
							int nsite=h_n.data[spv->n_idx(l, site)];
							double xx=(h_p.data[site].x-h_p.data[nsite].x)*(-1);
							double x=sqrt(xx*xx);
							if(x<L/2 && cluster_list_p[nsite]==CValp && newlist[nsite]==0)
							{
								cluster_array_size_p++;
								newlist[nsite]=1;
								cluster_p[cluster_array_size_p-1]=nsite;
								double dispNext=xx+disp[site];
								disp[nsite]=dispNext;
								if(dispNext>maxd)
								{
									maxd=dispNext;
									maxSite=nsite;
								}
								if(cBoundaries[nsite]==1 || cBoundaries[nsite]==5 || cBoundaries[nsite]==8)percCluster=true;
							}
							if(percCluster==true)break;
						}
					}
				}
			}

			if(percCluster==false && rightC==true)
			{
				vector <double> newposx(N, 0);
				double delta=L-maxd-0.01;
				for(int k=0; k<N; k++)
				{
					newposx[k]=h_p.data[k].x+delta;
					if(newposx[k]>=L)newposx[k]-=L;
				}
				vector <int> newleft(N,0);
				for(int k=0; k<N; k++)
				{
					for(int j=0; j<h_nn.data[k]; j++)
					{
						int site=h_n.data[spv->n_idx(j, k)];
						if(cluster_list_p[site]==CValp)
						{
							double xx=(newposx[k]-newposx[site]);
							double x=sqrt(xx*xx);
							if(x>L/2 && newposx[k]<newposx[site])
							{
								newleft[k]=1;
								break;
							}
						}
					}
				}
				cluster_array_size_p=1;
				cluster_p[cluster_array_size_p-1]=maxSite;
				vector <int> newlist(N, 0);
				newlist[maxSite]=1;
				while(cluster_array_size_p>0 && percCluster==false)
				{
					int site=cluster_p[cluster_array_size_p-1];
					cluster_array_size_p--;
					for(int l=0; l<h_nn.data[site]; l++)
					{
						int nsite=h_n.data[spv->n_idx(l, site)];
						double xx=(newposx[site]-newposx[nsite]);
						double x=sqrt(xx*xx);
						if(x<L/2 && cluster_list_p[nsite]==CValp && newlist[nsite]==0)
						{
							cluster_array_size_p++;
							newlist[nsite]=1;
							cluster_p[cluster_array_size_p-1]=nsite;
							if(newleft[nsite]==1)percCluster=true;
						}
						if(percCluster==true)break;
					}
				}
			}
			if(percCluster==true)break;
		}

		//check if largest cluster percolates vertically
		cluster_array_size_p=0;
		fill(cluster_p.begin(), cluster_p.end(), 0);
		fill(disp.begin(), disp.end(), 0);
		maxd=-100000;
		maxSite=0;
		for(int i=0; i<cbottom.size(); i++)
		{
			if(percCluster==true)break;
			int pos=cbottom[i];
			bool topC=false;
			if(cluster_list_p[pos]==CValp)
			{
				for(int k=0; k<ctop.size(); k++)
				{
					int ps=ctop[k];
					if(cluster_list_p[ps]==CValp)
					{
						topC=true;
						break;
					}
				}

				if(topC==true)
				{
					cluster_array_size_p=1;
					cluster_p[cluster_array_size_p-1]=pos;
					vector <int> newlist(N, 0);
					newlist[pos]=1;
					while(cluster_array_size_p>0 && percCluster==false)
					{
						int site=cluster_p[cluster_array_size_p-1];
						cluster_array_size_p--;
						for(int l=0; l<h_nn.data[site]; l++)
						{
							int nsite=h_n.data[spv->n_idx(l, site)];
							double yy=(h_p.data[site].y-h_p.data[nsite].y)*(-1);
							double y=sqrt(yy*yy);
							if(y<L/2 && cluster_list_p[nsite]==CValp && newlist[nsite]==0)
							{
								cluster_array_size_p++;
								newlist[nsite]=1;
								cluster_p[cluster_array_size_p-1]=nsite;
								double dispNext=yy+disp[site];
								disp[nsite]=dispNext;
								if(dispNext>maxd)
								{
									maxd=dispNext;
									maxSite=nsite;
								}
								if(cBoundaries[nsite]==2 || cBoundaries[nsite]==5 || cBoundaries[nsite]==6)percCluster=true;
							}
							if(percCluster==true)break;
						}
					}
				}
			}
			if(percCluster==false && topC==true)
			{
				vector <double> newposy(N, 0);
				double delta=L-maxd-0.01;
				for(int k=0; k<N; k++)
				{
					newposy[k]=h_p.data[k].y+delta;
					if(newposy[k]>=L)newposy[k]-=L;
				}
				vector <int> newbottom(N,0);
				for(int k=0; k<N; k++)
				{
					for(int j=0; j<h_nn.data[k]; j++)
					{
						int site=h_n.data[spv->n_idx(j, k)];
						if(cluster_list_p[site]==CValp)
						{
							double yy=(newposy[k]-newposy[site]);
							double y=sqrt(yy*yy);
							if(y>L/2 && newposy[k]<newposy[site])
							{
								newbottom[k]=1;
								break;
							}
						}
					}
				}
				cluster_array_size_p=1;
				cluster_p[cluster_array_size_p-1]=maxSite;
				vector <int> newlist(N, 0);
				newlist[maxSite]=1;
				while(cluster_array_size_p>0 && percCluster==false)
				{
					int site=cluster_p[cluster_array_size_p-1];
					cluster_array_size_p--;
					for(int l=0; l<h_nn.data[site]; l++)
					{
						int nsite=h_n.data[spv->n_idx(l, site)];
						double yy=(newposy[site]-newposy[nsite]);
						double y=sqrt(yy*yy);
						if(y<L/2 && cluster_list_p[nsite]==CValp && newlist[nsite]==0)
						{
							cluster_array_size_p++;
							newlist[nsite]=1;
							cluster_p[cluster_array_size_p-1]=nsite;
							if(newbottom[nsite]==1)percCluster=true;
						}
						if(percCluster==true)break;
					}
				}
			}
			if(percCluster==true)break;
		}

		if(percCluster==true)perc_avg_p[nn]++;
                for(int k=0; k<cluster_val_p; k++)cluster_size_avg_p[nn]+=cluster_size_p[k]/cluster_val_p;

		if(nn==b-1)
		{
			for(int k=0; k<cluster_val_p; k++)
			{
				int csize=cluster_size_p[k]/cdistSize;
				cluster_dist_p[csize]+=1/cluster_val_p;
			}
		}
            }

            if(cont%amostras==0 && cont>0){
		//filename="/home/diogo/MEGA/cenas/CFTC_Cluster/Results_tissues/Stripe/s39/MSD_"+to_string(start) +".txt";
		//filename="/home/diogo/MEGA/cenas/CFTC_Cluster/Results_tissues/Glass/GPU/s1/MSD_"+to_string(start)+".txt";
                //filename="/home/diogo/MEGA/cenas/CFTC_Cluster/Results_tissues/GRF/s"+to_string(scriptsIni)+"/T_Perc_"+to_string(start)+".txt";
		filename="/mnt/c/Users/Diogop/Documents/MEGA/cenas/CFTC_Cluster/Results_tissues/GRF/s"+to_string(scriptsIni)+"/T_Perc_"+to_string(start)+".txt";
                ofstream file1;
                file1.open (filename);
		int last=0;
                for(int n = 1; n<b; n++){
                    file1 << time[n] <<" "<< perc_avg_p[n]/amostras <<" "<< cluster_size_avg_p[n]/amostras <<" "<< phi_avg_p[n]/amostras<<" "<< (phi_avg_p_2[n]/amostras)-(phi_avg_p[n]/amostras)*(phi_avg_p[n]/amostras) <<endl;
                }
                file1.close();

                //filename="/home/destevao/CellGPUSubstrate/Results/"+to_string(ini)+"/data_distribution_"+to_string(value)+".txt";
		//filename="/home/diogo/MEGA/cenas/CFTC_Cluster/Results_tissues/GRF/s"+to_string(scriptsIni)+"/T_PercHist_"+to_string(start) +".txt";
		filename="/mnt/c/Users/Diogop/Documents/MEGA/cenas/CFTC_Cluster/Results_tissues/GRF/s"+to_string(scriptsIni)+"/T_PercHist_"+to_string(start) +".txt";
                file1.open (filename);
		for(int i=0; i<numBins; i++)
		{
		    file1 << i*cdistSize <<" "<< cluster_dist_p[i]/(amostras) << endl;
		}
                file1.close();

                //filename="/home/destevao/CellGPUSubstrate/Results/"+to_string(ini)+"/data_distribution_"+to_string(value)+".txt";
                //filename="/home/diogo/MEGA/cenas/CFTC_Cluster/Results_tissues/GRF/s"+to_string(scriptsIni)+"/T_PercFinal.txt";
		filename="/mnt/c/Users/Diogop/Documents/MEGA/cenas/CFTC_Cluster/Results_tissues/GRF/s"+to_string(scriptsIni)+"/T_PercFinal.txt";
		if(start==1)
		{
			file1.open (filename);
			file1.close();
		}
                file1.open (filename, std::ios_base::app);
                file1 << perc_avg_p[b-1]/amostras  <<" "<< cluster_size_avg_p[b-1]/amostras <<" "<< phi_avg_p[b-1]/amostras <<" "<< (phi_avg_p_2[b-1]/amostras)-(phi_avg_p[b-1]/amostras)*(phi_avg_p[b-1]/amostras)  << " " << v1 << " "<< v6 <<" "<< v8 <<" "<< v9 <<" "<< v10*sqrt(v1)/1024 <<" "<<tension_avg_ini[idx]<<" "<<tension_sigma_ini[idx]<<  endl;
                file1.close();
            }
        }
    }
    return 0;
}
