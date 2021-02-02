#include "std_include.h"

#include "Simulation.h"

#include "voronoiQuadraticEnergy.h"
#include "selfPropelledParticleDynamics.h"
#include "DatabaseNetCDFSPV.h"

#include <cstring>

using namespace std;

int main(){

    int cont=0;
    int last=126;
    int scriptsIni=2;
    vector <int> gather(0,0);
    for(int start=1; start<=last; start++){

      int v1, v2, v3, v4, v10, v11;
      double v5, v6, v7, v8, v9;
      //string access_dados="/home/diogo/MEGA/cenas/CFTC_Cluster/scripts35/dados.txt";
      //string access_dados="/media/hdd/Results/scripts39/dados.txt";
      //string access_dados="/media/hdd/Results/GPU/scripts5/dados.txt";
      //string access_dados="/home/diogo/MEGA/cenas/GPU/Results_glass/scripts"+to_string(scriptsIni)+"/dados.txt";
      string access_dados="/media/hdd/Results_GRF/scripts"+to_string(scriptsIni)+"/dados.txt";
      //string access_dados="/mnt/d/Results/GRF/scripts"+to_string(scriptsIni)+"/dados.txt";
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
      sprintf(data, "/media/hdd/Results_GRF/scripts%d/Job_%d/test_voronoi.nc", scriptsIni, (start-1)*v11+1);
      //sprintf(data, "/mnt/d/Results/GRF/scripts%d/Job_%d/test_voronoi.nc", scriptsIni, (start-1)*v11+1);
 
      SPVDatabaseNetCDF nc(v1,data,NcFile::ReadOnly);

        int N=v1;
        int a=N;
        int b=nc.GetNumRecs();
        int total_size=v2;
        int amostras=v11;
        //int start=1;
        int fim=start;
        vector <double> time(b,0);
	vector <double> tensions(b,0);
	vector <double> tension2(b*N*amostras, 0);
	vector <double> force(b,0);
	vector <double> fluid(b,0);
	vector <double> solid(b,0);
	vector <double> tensionsalt(b,0);
	vector <double> edgesalt(b,0);
        //double max=v6;
	//double min=max-v11;
        string filename;
        double dt=v5;
	int LA=last*amostras;
	int bins=100;
	vector <double> tHist(bins, 0);
	vector <double> taltHist(bins, 0);
	vector <double> pHist(bins, 0);
	vector <double> fHist(bins, 0);
	double numbinsp, numbinsp0, numbinsf, numbinstalt;
        double deltap, deltap0, deltaf, deltatalt;
	double p0Thresh=3.812;

        shared_ptr<VoronoiQuadraticEnergy> spv  = make_shared<VoronoiQuadraticEnergy>(N,1.0,4.0,true);
        FILE *f;
 
	double maxp=-1000000, maxp0=maxp, maxf=maxp, maxtalt=maxp, minp0=1000000, minp=minp0, minf=minp, mintalt=minp;
	int sample=0;
        for(int avg=(start-1)*amostras+1; avg<=fim*amostras; avg++){
	    sample++;
            cont=cont+1;
            double L;

            int curr_line=(avg-1)/amostras+1;           
            //int v1, v2, v3, v4, v9, v12, v14;
            //double v5, v6, v7, v8, v10, v11, v13;
            //string access_dados="/home/diogo/MEGA/cenas/CFTC_Cluster/scripts1/dados.txt";
	    //string access_dados="/media/hdd/Results/scripts25/dados.txt";
	    //string access_dados="/media/diogo/Elements/Tissues_results/scripts/dados2.txt";
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
            sprintf(dataname, "/media/hdd/Results_GRF/scripts%d/Job_%d/test_voronoi.nc", scriptsIni, avg);
	    //sprintf(dataname, "/mnt/d/Results/GRF/scripts%d/Job_%d/test_voronoi.nc", scriptsIni, avg);

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
		for(int ii=0; ii<N; ii++)spv->computeVoronoiForceCPU(ii);
		ArrayHandle<double2> h_f(spv->returnForces(),access_location::host,access_mode::read);
                ArrayHandle<double2> h_p(spv->cellPositions,access_location::host,access_mode::read);
                ArrayHandle<double2> h_app(spv->AreaPeriPreferences,access_location::host,access_mode::read);
                ArrayHandle<int> h_nn(spv->neighborNum,access_location::host,access_mode::read);
                ArrayHandle<int> h_n(spv->neighbors,access_location::host,access_mode::read);
                ArrayHandle<double2> h_ap(spv->returnAreaPeri(),access_location::host,access_mode::read);

		double tau=0;
		double fd=0;
		double sd=0;
		double fi=0;
		double talt=0;
		double edges=0;
		for(int i=0; i<N; i++)
		{
			for(int j=0; j<h_nn.data[i]; j++)
			{
				//cout<< h_nn.data[i]<<" "<<h_n.data[spv->n_idx(j,i)]<<" "<< spv->n_idx(j,i)  <<endl;
				//if(i<spv->n_idx(j,i))tau+=(h_ap.data[i].y-h_app.data[i].y)+(h_ap.data[h_n.data[spv->n_idx(j,i)]].y-h_app.data[h_n.data[spv->n_idx(j,i)]].y);
				int z=h_n.data[spv->n_idx(j,i)];
				if(z>i)
				{
					talt+=(h_ap.data[i].y-h_app.data[i].y)+(h_ap.data[z].y-h_app.data[z].y);
					edges++;
				}			
			}
			tau+=(h_ap.data[i].y-h_app.data[i].y);
			tension2[i+(sample-1)*N+nn*amostras*N]=(h_ap.data[i].y-h_app.data[i].y);
			fi+=sqrt(h_f.data[i].x*h_f.data[i].x+h_f.data[i].y*h_f.data[i].y);
			if(h_app.data[i].y>p0Thresh)fd++;
			else sd++;

			if(nn==b-1)
			{
				for(int j=0; j<h_nn.data[i]; j++)
				{
					int z=h_n.data[spv->n_idx(j,i)];
					if(z>i)
					{
						//double tt=(h_ap.data[i].y-h_app.data[i].y)+(h_ap.data[h_n.data[spv->n_idx(j,i)]].y-h_app.data[h_n.data[spv->n_idx(j,i)]].y);
						//if(tt>maxp)maxp=tt;
						//if(tt<minp)minp=tt;
		                                double ttalt=(h_ap.data[i].y-h_app.data[i].y)+(h_ap.data[z].y-h_app.data[z].y);
						if(ttalt>maxtalt)maxtalt=ttalt;
		                                if(ttalt<mintalt)mintalt=ttalt;
					}
				}
				double tt=(h_ap.data[i].y-h_app.data[i].y);
				if(tt>maxp)maxp=tt;
                                if(tt<minp)minp=tt;

				if(h_app.data[i].y>maxp0)maxp0=h_app.data[i].y;
				if(h_app.data[i].y<minp0)minp0=h_app.data[i].y;

				double ff=sqrt(h_f.data[i].x*h_f.data[i].x+h_f.data[i].y*h_f.data[i].y);
				if(ff>maxf)maxf=ff;
				if(ff<minf)minf=ff;
			}
		}
		tensions[nn]+=tau/N;
		fluid[nn]+=fd/N;
		solid[nn]+=sd/N;
		force[nn]+=fi/N;
		tensionsalt[nn]+=talt;
		edgesalt[nn]+=edges;
            }

            if(cont%amostras==0 && cont>0){
		//filename="/home/diogo/MEGA/cenas/CFTC_Cluster/Results_tissues/Stripe/s39/MSD_"+to_string(start) +".txt";
		//filename="/home/diogo/MEGA/cenas/CFTC_Cluster/Results_tissues/Glass/GPU/s1/MSD_"+to_string(start)+".txt";
                filename="/home/diogo/MEGA/cenas/CFTC_Cluster/Results_tissues/GRF/s"+to_string(scriptsIni)+"/network_"+to_string(start)+".txt";
		//filename="/mnt/c/Users/Diogop/Documents/MEGA/cenas/CFTC_Cluster/Results_tissues/GRF/s"+to_string(scriptsIni)+"/network_"+to_string(start)+".txt";
                ofstream file1;
                file1.open (filename);
		int last=0;

                for(int n = 1; n<b; n++)
		{	
			double tension_2=0;
			double tensionsFinal=tensions[n]/amostras;
			for(int kk=0; kk<amostras; kk++)
			{
				for(int ll=0; ll<N; ll++)
				{
					tension_2+=(tension2[ll+kk*N+n*N*amostras]-tensionsFinal)*(tension2[ll+kk*N+n*N*amostras]-tensionsFinal);
				}
			}
			tension_2/=(N*amostras-1);

    			file1 << time[n] <<" "<< tensionsFinal <<" "<< sqrt(tension_2) <<" "<< fluid[n]/amostras <<" "<<  solid[n]/amostras <<" "<< force[n]/amostras <<" "<< tensionsalt[n]/(amostras*edgesalt[n]) <<endl;
                }
                file1.close();

                deltap=maxp-minp;
                deltap0=maxp0-minp0;
		deltaf=maxf-minf;
		deltatalt=maxtalt-mintalt;
                numbinsp=deltap/bins;
                numbinsp0=deltap0/bins;
		numbinsf=deltaf/bins;
		numbinstalt=deltatalt/bins;
		if(numbinsp==0)numbinsp=1;
		if(numbinsp0==0)numbinsp0=1;
		if(numbinsf==0)numbinsf=1;
		if(numbinstalt==0)numbinstalt=1;
                fill(tHist.begin(), tHist.end(), 0);
                fill(pHist.begin(), pHist.end(), 0);
		fill(fHist.begin(), fHist.end(), 0);
		fill(taltHist.begin(), taltHist.end(), 0);
                int dp, dt, df, dtalt;
		double totalStd=0;
		double tensionAvg=tensions[b-1]/amostras;
		for(int aa=(start-1)*amostras+1; aa<=fim*amostras; aa++)
		{
                sprintf(dataname, "/media/hdd/Results_GRF/scripts%d/Job_%d/test_voronoi.nc", scriptsIni, aa);
		//sprintf(dataname, "/mnt/d/Results/GRF/scripts%d/Job_%d/test_voronoi.nc", scriptsIni, aa);
                SPVDatabaseNetCDF nc(N,dataname,NcFile::ReadOnly);
		b=nc.GetNumRecs();
		nc.ReadState(spv,b-1,true);
		for(int ii=0; ii<N; ii++)spv->computeVoronoiForceCPU(ii);
                ArrayHandle<double2> h_f(spv->returnForces(),access_location::host,access_mode::read);
                ArrayHandle<double2> h_p(spv->cellPositions,access_location::host,access_mode::read);
                ArrayHandle<double2> h_app(spv->AreaPeriPreferences,access_location::host,access_mode::read);
                ArrayHandle<int> h_nn(spv->neighborNum,access_location::host,access_mode::read);
                ArrayHandle<int> h_n(spv->neighbors,access_location::host,access_mode::read);
                ArrayHandle<double2> h_ap(spv->returnAreaPeri(),access_location::host,access_mode::read);

		double std=0;
                for(int i=0; i<N; i++)
                {
                        dp=(h_app.data[i].y-minp0)/numbinsp0;
			pHist[dp]++;
			for(int j=0; j<h_nn.data[i]; j++)
                        {
				int z=h_n.data[spv->n_idx(j,i)];
				if(z>i)
				{
					//dt=((h_ap.data[i].y-h_app.data[i].y)+(h_ap.data[h_n.data[spv->n_idx(j,i)]].y-h_app.data[h_n.data[spv->n_idx(j,i)]].y)-minp)/numbinsp;
					//tHist[dt]++;
					dtalt=((h_ap.data[i].y-h_app.data[i].y)+(h_ap.data[z].y-h_app.data[z].y)-mintalt)/numbinstalt;
					taltHist[dtalt]++;
				}
                        }
			double tensionVal=(h_ap.data[i].y-h_app.data[i].y);
			dt=(tensionVal-minp)/numbinsp;
			tHist[dt]++;
			std+=(tensionVal-tensionAvg)*(tensionVal-tensionAvg);

			double ff=sqrt(h_f.data[i].x*h_f.data[i].x+h_f.data[i].y*h_f.data[i].y);
			df=(ff-minf)/numbinsf;
			fHist[df]++;
                }
		std=std/(N-1);
                std=sqrt(std);
		totalStd+=std;
		}

                //filename="/home/destevao/CellGPUSubstrate/Results/"+to_string(ini)+"/data_distribution_"+to_string(value)+".txt";
                //ofstream file1
		filename="/home/diogo/MEGA/cenas/CFTC_Cluster/Results_tissues/GRF/s"+to_string(scriptsIni)+"/hist_"+to_string(start) +".txt";
                //filename="/home/diogo/CellGPU/SPV/Results/s"+to_string(scriptsIni)+"/Diff.txt";
		//filename="/mnt/c/Users/Diogop/Documents/MEGA/cenas/CFTC_Cluster/Results_tissues/GRF/s"+to_string(scriptsIni)+"/hist_"+to_string(start) +".txt";
                file1.open (filename);
		for(int i=0; i<bins; i++)
		{
		    file1 << i*numbinsp+minp  <<" "<< tHist[i]/(amostras*N) << " " << i*numbinsp0+minp0 << " "<< pHist[i]/(amostras*N) <<" "<< i*numbinsf+minf <<" "<< fHist[i]/(amostras*N) <<" "<< i*numbinstalt+mintalt <<" "<< taltHist[i]/(amostras*3*N) << endl;
		}
                file1.close();

		filename="/home/diogo/MEGA/cenas/CFTC_Cluster/Results_tissues/GRF/s"+to_string(scriptsIni)+"/std.txt";
                if(start==1)
                {
                        file1.open (filename);
                        file1.close();
                }
                file1.open (filename, std::ios_base::app);
		file1<<totalStd/amostras<<" "<<v9<<" "<< tensionAvg<<endl;
		file1.close();

            }
        }
    }

    cout<<"not finished:"<<endl;
    for(int i=0; i<gather.size(); i++)cout<<gather[i]<<endl;
    return 0;
}
