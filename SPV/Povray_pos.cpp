#include "std_include.h"

#include "Simulation.h"

#include "voronoiQuadraticEnergy.h"
#include "selfPropelledParticleDynamics.h"
#include "DatabaseNetCDFSPV.h"

#include <cstring>

using namespace std;


int main(){

  int last=126;
  int scriptsIni=2;
    for(int start=1; start<=last; start++){

      int v1, v2, v3, v4, v10, v11;
      double v5, v6, v7, v8, v9;
      //string access_dados="/home/diogo/MEGA/cenas/GPU/Results/scripts27/dados.txt";
      //string access_dados="/mnt/d/Results/GRF/scripts"+to_string(scriptsIni)+"/dados.txt";
      string access_dados="/media/hdd/Results_GRF/scripts"+to_string(scriptsIni)+"/dados.txt";
      //string access_dados="/media/hdd/Results/scripts38/dados.txt";
      //string access_dados="/media/hdd/Results/GPU/scripts11/dados.txt";
      //string access_dados="/home/diogo/MEGA/cenas/CellGPU_Substrate2/dados.txt";

      ifstream file(access_dados);

      for(int line=0; line<start; line++){
	file>> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7 >> v8 >> v9 >> v10 >> v11;
      }

      file.close();

      char data[1000000]; 
      //sprintf(data, "/mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV/build2/test_voronoi.nc");
      //sprintf(data, "/mnt/d/Results/GRF/scripts%d/Job_%d/test_voronoi.nc",  scriptsIni, (start-1)*v11+1);
      sprintf(data, "/media/hdd/Results_GRF/scripts%d/Job_%d/test_voronoi.nc",  scriptsIni, (start-1)*v11+1);
      //sprintf(data, "/media/hdd/Results/scripts38/Job_%d/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", (start-1)*v14+1, v9, v10, v11, v6, v7, v8);
      //sprintf(data, "/media/hdd/Results/GPU/scripts11/Job_%d/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", (start-1)*v14+1, v9, v10, v11, v6, v7, v8);
      //sprintf(data, "/home/diogo/MEGA/cenas/CellGPU_Substrate2/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", v9, v10, v11, v6, v7, v8);

      SPVDatabaseNetCDF nc(v1,data,NcFile::ReadOnly);
      
        int N=v1;
        int a=N;
        int b=nc.GetNumRecs();
	
        int amostras=1;//v14;
        //int start=2;
        int fim=start;
	int time_total=b;
	int start_time=b-1;
        string filename;
	int val=1;
        double p0Thresh=3.8;
        double pThresh=0.08;

	EOMPtr spp = make_shared<selfPropelledParticleDynamics>(N);
        shared_ptr<VoronoiQuadraticEnergy> spv  = make_shared<VoronoiQuadraticEnergy>(N, false);

        FILE *f;
	    	        	
	for(int nn = start_time; nn<time_total; nn++){

	  cout<<"NEW TIME----------------"<<" "<<nn<<endl;
	    
	  for(int avg=(start-1)*amostras+1; avg<=(start-1)*amostras+1; avg++){
            double L;
	    
            int curr_line=(avg-1)/amostras+1;           
            //int v1, v2, v3, v4, v9, v12, v14;
            //double v5, v6, v7, v8, v10, v11, v13;
            //string access_dados="/home/diogo/MEGA/cenas/CFTC_Cluster/scripts20/dados.txt";
	    file.open(access_dados);

            for(int line=0; line<curr_line; line++){
	      file>> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7 >> v8 >> v9 >> v10 >> v11;
            }

            file.close();

            string add_JOB_num=to_string(avg);
            //if(avg<10)add_JOB_num="00"+to_string(avg);
            //else if(avg<100)add_JOB_num="0"+to_string(avg);
            //else add_JOB_num=to_string(avg);
           
            char dataname[1000000]; 
            //sprintf(dataname, "/mnt/c/Users/Diogop/Documents/MEGA/cenas/CellGPU/SPV/build2/test_voronoi.nc");
	    //sprintf(dataname, "/mnt/d/Results/GRF/scripts%d/Job_%d/test_voronoi.nc",  scriptsIni, (start-1)*v11+1);
            sprintf(dataname, "/media/hdd/Results_GRF/scripts%d/Job_%d/test_voronoi.nc",  scriptsIni, (start-1)*v11+1);
	    //sprintf(dataname, "/media/hdd/Results/scripts38/Job_%d/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", avg, v9, v10, v11, v6, v7, v8);
	    //sprintf(dataname, "/media/hdd/Results/GPU/scripts10/Job_%d/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", avg, v9, v10, v11, v6, v7, v8);
            //sprintf(dataname, "/home/diogo/MEGA/cenas/CellGPU_Substrate2/test_voronoi_St%d_T%g_Iv%g_P%g_A%g_V%g.nc", v9, v10, v11, v6, v7, v8);
	    
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
            }

	    ncdat.ReadState(spv, nn, true);
	
	    //filename="/home/diogo/MEGA/cenas/CFTC_Cluster/Results_tissues/Stripe/s32/POV/snap00"+to_string(val)+".pov";
	    //filename="/home/diogo/MEGA/cenas/CellGPU_Substrate/snap00"+to_string(val)+".pov";
	    //filename="snap00"+to_string(val)+".pov";
            //filename="/mnt/c/Users/Diogop/Documents/MEGA/cenas/CFTC_cluster/Results_tissues/GRF/s"+to_string(scriptsIni)+"/Images/snap00"+to_string(start)+".pov";
	    filename="/home/diogo/MEGA/cenas/CFTC_Cluster/Results_tissues/GRF/s"+to_string(scriptsIni)+"/Images/snap00"+to_string(start)+".pov";
	    ofstream file1;
	    file1.open (filename);
	    //filename="/mnt/c/Users/Diogop/Documents/MEGA/cenas/CFTC_cluster/Results_tissues/GRF/s"+to_string(scriptsIni)+"/Images/Csnap00"+to_string(start)+".pov";
	    filename="/home/diogo/MEGA/cenas/CFTC_Cluster/Results_tissues/GRF/s"+to_string(scriptsIni)+"/Images/Csnap00"+to_string(start)+".pov";
            ofstream file2;
            file2.open (filename);
            //filename="/mnt/c/Users/Diogop/Documents/MEGA/cenas/CFTC_cluster/Results_tissues/GRF/s"+to_string(scriptsIni)+"/Images/Tsnap00"+to_string(start)+".pov";
	    filename="/home/diogo/MEGA/cenas/CFTC_Cluster/Results_tissues/GRF/s"+to_string(scriptsIni)+"/Images/Tsnap00"+to_string(start)+".pov";
	    ofstream file3;
	    file3.open (filename);

            ArrayHandle<double2> h_p(spv->cellPositions,access_location::host,access_mode::read);
            ArrayHandle<double2> h_app(spv->AreaPeriPreferences,access_location::host,access_mode::read);
	    ArrayHandle<int> h_nn(spv->neighborNum,access_location::host,access_mode::read);
	    ArrayHandle<int> h_n(spv->neighbors,access_location::host,access_mode::read);
	    ArrayHandle<double2> h_ap(spv->returnAreaPeri(),access_location::host,access_mode::read);

            double delta=0;
	    double deltap=0;
            double minp0=100000000000;
            double maxp0=0;
	    double minp=10000000000;
	    double maxp=0;
            int Ncells=N;
            for(int i=0; i<Ncells; i++)
            {
                if(h_app.data[i].y>maxp0)maxp0=h_app.data[i].y;
                if(h_app.data[i].y<minp0)minp0=h_app.data[i].y;
		if(abs(h_ap.data[i].y-h_app.data[i].y)>maxp)maxp=abs(h_ap.data[i].y-h_app.data[i].y);
                if(abs(h_ap.data[i].y-h_app.data[i].y)<minp)minp=abs(h_ap.data[i].y-h_app.data[i].y);
            }
            delta=maxp0-minp0;
	    deltap=maxp-minp;

	    for (int i = 0; i < Ncells; ++i)
	      {
		//get Delaunay neighbors of the cell
		int neigh = h_nn.data[i];
		vector<int> ns(neigh);
		for (int nn = 0; nn < neigh; ++nn)
		  {
		    ns[nn]=h_n.data[spv->n_idx(nn,i)];
		  };

		double2 circumcent;
		double2 nnextp,nlastp;
		double2 pi = h_p.data[i];
		double2 rij, rik;
		double2 voro;
		double2 voron;
		double2 voroi;
		double color1=0.01, color2=(h_app.data[i].y-minp0)/delta, color3=0.01;
		if(delta<=1e-12)color2=0.01;
		double color4=0, color5=0, color6=0;
		double diff=abs(h_ap.data[i].y-h_app.data[i].y);
		double color7=(diff-minp)/deltap;
		if(deltap<=1e-12)color7=0.01;
	        double color8=color7,color9=color7;
		//color2=0.1+((h_ap.data[i].y-(3.9-v11))/v11)*(0.175-0.1);

		if(h_app.data[i].y>p0Thresh)
		  {
		    //color1=0.898039216;
		    //color2=0.988235294;
		    //color3=0.760784314;
		    color4=0.99;
		    color5=0.99;
		    color6=0.99;

		    color1=0.99;
                    color2=0.99;
                    color3=0.99;
		  }
		else
		  {
		    //color1=0.615686275;
		    //color2=0.878431373;
		    //color3=0.678431373;
		    color4=0.01;
		    color5=0.01;
		    color6=0.01;

		    color1=0.01;
                    color2=0.01;
                    color3=0.01;
		  }

		if(h_ap.data[i].y>p0Thresh)
		{
                    color7=0.99;
                    color8=0.99;
                    color9=0.99;

		    if(color1<0.1)
		    {
			    color1=0.01;
			    color2=0.5;
			    color3=0.01;
		    }
		    else
		    {
			    color1=0.5;
                            color2=0.01;
                            color3=0.01;
		    }
		}
		else
		{
                    color7=0.01;
                    color8=0.01;
                    color9=0.01;
		}

		nlastp = h_p.data[ns[ns.size()-1]];
		spv->Box->minDist(nlastp,pi,rij);
		for (int zz = 0; zz < neigh; ++zz)
		  {
		    if(zz>0)voron=voro;
		    nnextp = h_p.data[ns[zz]];
		    spv->Box->minDist(nnextp,pi,rik);
		    Circumcenter(rij,rik,circumcent);
		    voro = circumcent;
		    rij=rik;

		    if(zz>0){

		      //if(h_p.data[i].x<sqrt(N)/2+5 && h_p.data[i].x>sqrt(N)/2-5 && h_p.data[i].y>0 && h_p.data[i].y<5)
		      //{
		      file1 << setprecision(15) << "Voro(" << h_p.data[i].x+voro.x << "," << h_p.data[i].y+voro.y << "," << h_p.data[i].x << "," << h_p.data[i].y << "," << h_p.data[i].x+voron.x << "," << h_p.data[i].y+voron.y << "," << color1 << "," << color2<< "," << color3  << ")"  << endl;
                      file2 << setprecision(15) << "Voro(" << h_p.data[i].x+voro.x << "," << h_p.data[i].y+voro.y << "," << h_p.data[i].x << "," << h_p.data[i].y << "," << h_p.data[i].x+voron.x << "," << h_p.data[i].y+voron.y << "," << color4 << "," << color5<< "," << color6  << ")"  << endl;
                      file3 << setprecision(15) << "Voro(" << h_p.data[i].x+voro.x << "," << h_p.data[i].y+voro.y << "," << h_p.data[i].x << "," << h_p.data[i].y << "," << h_p.data[i].x+voron.x << "," << h_p.data[i].y+voron.y << "," << color7 << "," << color8<< "," << color9  << ")"  << endl;
		      //}
		    }
		    else voroi=voro;
		  }
		//if(h_p.data[i].x<sqrt(N)/2+5 && h_p.data[i].x>sqrt(N)/2-5 && h_p.data[i].y>0 && h_p.data[i].y<5)
		//{
		/*
		  if(i==15){
		  color1=0.01;
		  color2=0.01;
		  color3=0.25;}*/

		file1 << setprecision(15) << "Voro(" << h_p.data[i].x+voro.x << "," << h_p.data[i].y+voro.y << "," << h_p.data[i].x << "," << h_p.data[i].y << "," << h_p.data[i].x+voroi.x << "," << h_p.data[i].y+voroi.y << "," << color1 << "," << color2 << "," << color3 << ")"  << endl;
		file2 << setprecision(15) << "Voro(" << h_p.data[i].x+voro.x << "," << h_p.data[i].y+voro.y << "," << h_p.data[i].x << "," << h_p.data[i].y << "," << h_p.data[i].x+voroi.x << "," << h_p.data[i].y+voroi.y << "," << color4 << "," << color5 << "," << color6 << ")"  << endl;
		file3 << setprecision(15) << "Voro(" << h_p.data[i].x+voro.x << "," << h_p.data[i].y+voro.y << "," << h_p.data[i].x << "," << h_p.data[i].y << "," << h_p.data[i].x+voroi.x << "," << h_p.data[i].y+voroi.y << "," << color7 << "," << color8 << "," << color9 << ")"  << endl;
		//}
	      }
	    val++;

	    file1.close();
	    file2.close();
	    file3.close();
	  }
	}	
    }
    return 0;
}
