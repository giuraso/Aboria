/*
 * sphdem.cpp
 * 
 * Copyright 2014 Martin Robinson
 *
 * This file is part of Aboria.
 *
 * Aboria is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Aboria is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Aboria.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 7 Feb 2014
 *      Author: robinsonm
 */

#include "sphdem.h"
#include "Visualisation.h"

#include <fstream>
#include <string>
#include <sstream>

#include <vtkFloatArray.h>

std::string make_output_filename(size_t timeout) {
   std::ostringstream ss;
   ss << "output_" << timeout << ".data";
   return ss.str();
}

int calc_part_num(double diam, double L){
	int num = int(L/diam);
	return num;
}

int main(int argc, char **argv) {
	ofstream demfile;
    demfile.open ("velocity.txt");

	ofstream sphfile;
    sphfile.open ("sphcheck.txt");

    ofstream timefile;
    timefile.open ("time.txt");
  
	auto dem = DemType::New();
	auto sph = SphType::New();
	auto params = ptr<Params>(new Params());


	const int timesteps = 20000;
	const int nout = 2000;
	const int timesteps_per_out = timesteps/nout;
	const double Lx = 6/1000.0;
	const int nx = 30;
	const double Lz = 6/1000.0;
	const int nz = 30;
	
	 /* dem parameters
	 */
	params->dem_diameter = 0.0001;
	params->dem_gamma = 0.0;
	params->dem_k = 10.0;
	params->dem_vol = (1.0/6.0)*PI*pow(params->dem_diameter,3);
	const double dem_dens = 2500.0;
	params->dem_mass = params->dem_vol*dem_dens;
	const double dem_min_reduced_mass = 0.5*params->dem_mass;
	params->dem_dt = (1.0/50.0)*PI/sqrt(params->dem_k/dem_min_reduced_mass-pow(0.5*params->dem_gamma/dem_min_reduced_mass,2));

	/*

	/*
	 * sph parameters
	 */
	params->sph_hfac = 1.5;
	params->sph_visc = 8.9e-07;
	params->sph_refd = 1000.0;
	params->sph_dens = 1000.0;
	params->sph_gamma = 7;
	const double VMAX = 2.0*sqrt(2*9.81*Lz);
	const double CSFAC = 10.0;
	params->sph_spsound = CSFAC*VMAX;
	params->sph_prb = pow(params->sph_refd/params->sph_dens,params->sph_gamma-1.0)*pow(params->sph_spsound,2)*params->sph_refd/params->sph_gamma;
	const double psep = Lx/nx;
	params->sph_dt = std::min(0.25*params->sph_hfac*psep/params->sph_spsound,0.125*pow(params->sph_hfac*psep,2)/params->sph_visc);
	params->sph_mass = params->sph_dens*pow(psep,NDIM);

	std::cout << "sph dt= "<<params->sph_dt<<std::endl;
	std::cout << "dem dt= "<<params->dem_dt<<std::endl;

	const int wait_dem=8000;

	params->dem_time_drop = params->sph_dt*wait_dem;

	params->sph_time_damping = params->sph_dt*(800);	
	params->time = 0;


	/*
	 * define domain / geometry
	 */
	auto dem_geometry = [params](DemType::Value& i) {
		Vect3d acceleration;
		acceleration << 0,0,-9.81;
		REGISTER_DEM_PARTICLE(i);
		const double dem_diameter = params->dem_diameter;
		const double dem_k = params->dem_k;
		const double dem_gamma = params->dem_gamma;
		const double dem_mass = params->dem_mass;
		const double dem_vol = params->dem_vol;

		const double overlap = dem_diameter/2.0-r[2];
		if (overlap>0) {
			const double overlap_dot = -v[2];
			const Vect3d normal(0,0,1);
			acceleration += (dem_k*overlap + dem_gamma*overlap_dot)*normal/dem_mass;
		}
		
		return acceleration;
	};

	auto sph_geometry = [](SphType::Value& i) {
		Vect3d acceleration;
		acceleration << 0,0,-9.81;
		return acceleration;
	};

	const Vect3d min(0,0,-3.0*psep);
	const Vect3d max(Lx,Lx,Lz);
	const Vect3b periodic(true,true,false);

	/*
	 * create sph and dem particles
	 */
	sph->create_particles_grid(min,max,Vect3i(nx,nx,nz+3),[psep,params](SphType::Value& i) {
		REGISTER_SPH_PARTICLE(i);
		h = params->sph_hfac*psep;
		omega = 1.0;
		v << 0,0,0;
		v0 << 0,0,0;
		dddt = 0;
		e = 1;
		rho = params->sph_dens;
		f << 0,0,0;
		f0 << 0,0,0;

		if (r[2]<0) {
			fixed = true;
		} else {
			fixed = false;
		}
	});

	

	int nside=34;
	int nheight=13;
	double npart=nside*nside*nheight;
	double separ=(Lx-params->dem_diameter)/(nside-1);
	double maxz=Lz/2+(nheight-1)*separ+params->dem_diameter;
	double voidage=(Lx*Lx*((nheight-1)*separ+params->dem_diameter)-npart*params->dem_vol)/(Lx*Lx*((nheight-1)*separ+params->dem_diameter));
	std::cout<<params->dem_vol<<" "<<voidage<<endl;
	double minx=params->dem_diameter/2;
	double miny=params->dem_diameter/2;
	double minz=Lz/2-params->dem_diameter/2;
	double maxx=Lx-params->dem_diameter/2;
	double maxy=Lx-params->dem_diameter/2;
	
	const Vect3d minC(minx,miny,minz);
	const Vect3d maxC(maxx,maxy,maxz);
	std::cout<<maxx<<" "<<maxy<<" "<<maxz<<endl;
	std::cout<<minx<<" "<<miny<<" "<<minz<<endl;

dem->create_particles_grid(minC,maxC,Vect3i(nside,nside,nheight),[params](DemType::Value& i) {
		REGISTER_DEM_PARTICLE(i);
		v << 0,0,0;
		v0 << 0,0,0;
		f << 0,0,0;
		f0 << 0,0,0;
	});


	/*
* setup output stuff
*

auto sph_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
auto dem_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
sph->copy_to_vtk_grid(sph_grid);
dem->copy_to_vtk_grid(dem_grid);
//Visualisation::vtkWriteGrid("vis/at_start_sph",0,sph_grid);
//Visualisation::vtkWriteGrid("vis/at_start_dem",0,dem_grid);


//	/Visualisation vis(min,max);
//	vis.glyph_points(sph_grid);
//	vis.start_render_loop();
	vtkSmartPointer<vtkFloatArray> vis_v = vtkSmartPointer<vtkFloatArray>::New();
	vtkSmartPointer<vtkFloatArray> vis_f = vtkSmartPointer<vtkFloatArray>::New();
	vtkSmartPointer<vtkFloatArray> vis_rho = vtkSmartPointer<vtkFloatArray>::New();
	vtkSmartPointer<vtkFloatArray> vis_eps = vtkSmartPointer<vtkFloatArray>::New();

	vtkSmartPointer<vtkFloatArray> vis_v_dem = vtkSmartPointer<vtkFloatArray>::New();
	vtkSmartPointer<vtkFloatArray> vis_f_dem = vtkSmartPointer<vtkFloatArray>::New();

	vis_v->SetName("v");
	vis_f->SetName("f");
	vis_rho->SetName("rho");
	vis_eps->SetName("eps");

	vis_v->SetNumberOfComponents(3);
	vis_f->SetNumberOfComponents(3);

	vis_v_dem->SetName("vdem");
	vis_f_dem->SetName("fdem");

	vis_v_dem->SetNumberOfComponents(3);
	vis_f_dem->SetNumberOfComponents(3);

	sph_grid->GetPointData()->AddArray(vis_v);
	sph_grid->GetPointData()->AddArray(vis_f);
	sph_grid->GetPointData()->AddArray(vis_rho);
	sph_grid->GetPointData()->AddArray(vis_eps);

	dem_grid->GetPointData()->AddArray(vis_v_dem);
	dem_grid->GetPointData()->AddArray(vis_f_dem);
	

std::cout << "starting...."<<std::endl;
sph->init_neighbour_search(min,max,2*params->sph_hfac*psep,periodic);
dem->init_neighbour_search(min,max,params->dem_diameter,periodic);
*/
int nSPHpart=nx*nx*(nz+3);
size_t outputfile=0;
	for (int i = 0; i < nout; ++i) {
		for (int k = 0; k < timesteps_per_out; ++k) {
			//std::this_thread::sleep_for(std::chrono::seconds(1));
	
			sphdem(sph,dem,params,sph_geometry,dem_geometry);
			std::cout <<"run"<<std::endl;
		}
		timefile<<params->time<<" "<<i<<" "<<params->sph_dt<<endl;
		std::cout <<"iteration ok"<<i<<std::endl;
		//vis.stop_render_loop();
/*		sph->copy_to_vtk_grid(sph_grid);
		vis_v->SetNumberOfTuples(sph->size());
		vis_f->SetNumberOfTuples(sph->size());
		vis_rho->SetNumberOfTuples(sph->size());
		vis_eps->SetNumberOfTuples(sph->size());		

		vis_v_dem->SetNumberOfTuples(dem->size());
		vis_f_dem->SetNumberOfTuples(dem->size());

*/		
		int ii = 0;
		for (SphType::iterator i = sph->begin(); i != sph->end(); i++,ii++) {
			REGISTER_SPH_PARTICLE((*i));
/*			vis_v->SetTuple3(ii,v[0],v[1],v[2]);
			vis_f->SetTuple3(ii,r[0],r[1],r[2]);
			vis_rho->SetValue(ii,rho);
			vis_eps->SetValue(ii,e);
*/		}
			if ((i>1100) && ((i%300)==0)){
			outputfile++;
			
			FILE *file = fopen(make_output_filename(outputfile).c_str(), "w");
			fclose(file);
					
			ofstream outfile;
			outfile.open(make_output_filename(outputfile).c_str());
			outfile<<nSPHpart<< " "<<params->time<< " "<< min[0]<< " "<< min[1]<< " "<< min[2]<< " "<< max[0]<< " "<< max[1]<< " "<< max[2]<<endl;
			
			int ii = 0;
			for (SphType::iterator i = sph->begin(); i != sph->end(); i++,ii++) {
				REGISTER_SPH_PARTICLE((*i));
			outfile<<r[0]<<" "<<r[1]<< " "<<r[2]<<" ";
			outfile<<v[0]<<" "<<v[1]<< " "<<v[2]<< " ";				
			outfile<<psep/2<<" ";
			outfile<<0<< " "<<0<<" "<<0<<" ";
			outfile<<0<<" "<<0<<" "<<0<<" ";
			outfile<<e<<endl;
			}
				
			outfile.close();
			}
		
		ii = 0;
		for (DemType::iterator i = dem->begin(); i != dem->end(); i++,ii++) {
			REGISTER_DEM_PARTICLE((*i));
//			vis_v_dem->SetTuple3(ii,v[0],v[1],v[2]);
//			vis_f_dem->SetTuple3(ii,f[0],f[1],f[2]);
                if (ii==1){
		                demfile << v[2] <<" "<< f0[0] <<" "<< f0[1] <<" "<< f0[2];
				                demfile <<" " << r[0] <<" "<< r[1] <<" "<< r[2] << endl;}
						                

}
		
//sph->copy_to_vtk_grid(sph_grid);
//dem->copy_to_vtk_grid(dem_grid);

//Visualisation::vtkWriteGrid("vis/sph",i,sph_grid);
//Visualisation::vtkWriteGrid("vis/dem",i,dem_grid);
		//vis.restart_render_loop();
	

		//vis.stop_render_loop();
		//vis.restart_render_loop();
	}
/*
auto sph_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
auto dem_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
sph->copy_to_vtk_grid(sph_grid);
dem->copy_to_vtk_grid(dem_grid);
Visualisation::vtkWriteGrid("vis/at_start_sph",0,sph_grid);
Visualisation::vtkWriteGrid("vis/at_start_dem",0,dem_grid);


std::cout << "starting...."<<std::endl;
sph->init_neighbour_search(min,max,2*params->sph_hfac*psep,periodic);
dem->init_neighbour_search(min,max,params->dem_diameter,periodic);

for (int i = 0; i < nout; ++i) {
for (int k = 0; k < timesteps_per_out; ++k) {
//std::this_thread::sleep_for(std::chrono::seconds(1));
sphdem(sph,dem,params,sph_geometry,dem_geometry);

}
std::cout <<"iteration "<<i<<std::endl;
		myfile << params->dem_diameter << endl;
sph->copy_to_vtk_grid(sph_grid);
dem->copy_to_vtk_grid(dem_grid);
Visualisation::vtkWriteGrid("vis/sph",i,sph_grid);
Visualisation::vtkWriteGrid("vis/dem",i,dem_grid);
}
*/

	demfile.close();
	sphfile.close();
	timefile.close();
	ofstream parfile;
    parfile.open ("params.txt");
		parfile << params->sph_dt <<" "<< params->time <<" "<< timesteps <<" "<< nout<<endl;
		parfile <<params->sph_hfac  <<" "<< nx <<" "<< nz <<" "<< Lx <<" "<< Lz << endl;
		parfile << params->sph_time_damping <<" "<< params->dem_time_drop << endl;
	parfile.close();

	
}
