/*
 * sphdem.h
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

#ifndef SPHDEM_H_
#define SPHDEM_H_

#include "Aboria.h"
using namespace Aboria;

#include "sph_common.h"
#include <tuple>



enum {DEM_FORCE, DEM_VELOCITY, DEM_VELOCITY0};
typedef std::tuple<Vect3d,Vect3d,Vect3d> DemTuple;
typedef Particles<DemTuple> DemType;

#define GET_TUPLE(type,name,position,particle) type& name = std::get<position>(particle.get_data())
#define REGISTER_DEM_PARTICLE(particle) \
				const Vect3d& r = particle.get_position(); \
				GET_TUPLE(Vect3d,f,DEM_FORCE,particle); \
				GET_TUPLE(Vect3d,v,DEM_VELOCITY,particle); \
				GET_TUPLE(Vect3d,v0,DEM_VELOCITY0,particle)
#define REGISTER_NEIGHBOUR_DEM_PARTICLE(tuple) \
				const Vect3d& dx = std::get<1>(tuple); \
				const DemType::Value& j = std::get<0>(tuple); \
				const Vect3d& rj = j.get_position(); \
				const GET_TUPLE(Vect3d,fj,DEM_FORCE,j); \
				const GET_TUPLE(Vect3d,vj,DEM_VELOCITY,j); \
				const GET_TUPLE(Vect3d,v0j,DEM_VELOCITY0,j)


enum {SPH_FORCE, SPH_VELOCITY, SPH_VELOCITY0, SPH_DENS, SPH_POROSITY, SPH_H, SPH_DDDT,SPH_PDR2,SPH_OMEGA,SPH_FIXED};
typedef std::tuple<Vect3d,Vect3d,Vect3d,double,double,double,double,double,double,bool> SphTuple;
typedef Particles<SphTuple> SphType;

#define REGISTER_SPH_PARTICLE(particle) \
				const Vect3d& r = particle.get_position(); \
				GET_TUPLE(bool,fixed,SPH_FIXED,particle); \
				GET_TUPLE(Vect3d,f,SPH_FORCE,particle); \
				GET_TUPLE(Vect3d,v,SPH_VELOCITY,particle); \
				GET_TUPLE(Vect3d,v0,SPH_VELOCITY0,particle); \
				GET_TUPLE(double,rho,SPH_DENS,particle); \
				GET_TUPLE(double,e,SPH_POROSITY,particle); \
				GET_TUPLE(double,h,SPH_H,particle); \
				GET_TUPLE(double,dddt,SPH_DDDT,particle); \
				GET_TUPLE(double,pdr2,SPH_PDR2,particle); \
				GET_TUPLE(double,omega,SPH_OMEGA,particle)

#define REGISTER_NEIGHBOUR_SPH_PARTICLE(tuple) \
				const Vect3d& dx = std::get<1>(tuple); \
				const SphType::Value& j = std::get<0>(tuple); \
				const Vect3d& rj = j.get_position(); \
				const GET_TUPLE(bool,fixedj,SPH_FIXED,j); \
				const GET_TUPLE(Vect3d,fj,SPH_FORCE,j); \
				const GET_TUPLE(Vect3d,vj,SPH_VELOCITY,j); \
				const GET_TUPLE(Vect3d,v0j,SPH_VELOCITY0,j); \
				const GET_TUPLE(double,rhoj,SPH_DENS,j); \
				const GET_TUPLE(double,ej,SPH_POROSITY,j); \
				const GET_TUPLE(double,hj,SPH_H,j); \
				const GET_TUPLE(double,dddtj,SPH_DDDT,j); \
				const GET_TUPLE(double,pdr2j,SPH_PDR2,j); \
				const GET_TUPLE(double,omegaj,SPH_OMEGA,j)

struct Params {
	double sph_dt,sph_mass,sph_hfac,sph_visc,sph_refd,sph_gamma,sph_spsound,sph_prb;
	double dem_dt,dem_mass,dem_diameter,dem_k,dem_gamma;
};


template<typename GeometryType>
void dem_start(ptr<DemType> dem,
		ptr<Params> params,
		GeometryType geometry) {

	const double dt = params->dem_dt;
	const double dem_diameter = params->dem_diameter;
	const double dem_k = params->dem_k;
	const double dem_gamma = params->dem_gamma;
	const double dem_mass = params->dem_mass;


	dem->update_positions(dem->begin(),dem->end(),[dt](DemType::Value& i) {
		REGISTER_DEM_PARTICLE(i);

		v0 = v + dt/2*f;
		v += dt * f;
		return r + dt * v0;
	});

	std::for_each(dem->begin(),dem->end(),[&geometry,dem,dem_k,dem_gamma,dem_mass,dem_diameter](DemType::Value& i) {
		REGISTER_DEM_PARTICLE(i);

		f << 0,0,0;
		f = f + geometry(i);

		for (auto tpl: i.get_neighbours(dem)) {
			REGISTER_NEIGHBOUR_DEM_PARTICLE(tpl);

			if (i.get_id()==j.get_id()) continue;

			const double r = dx.norm();
			const double overlap = dem_diameter-r;
			if (overlap>0) {
				const Vect3d normal = dx/r;
				const Vect3d dv = v-vj;
				const double overlap_dot = -dv.dot(normal);
				f += (dem_k*overlap + dem_gamma*overlap_dot)*normal/dem_mass;
			}
		}

	});

}

template<typename GeometryType>
void dem_end(ptr<DemType> dem,
		ptr<Params> params,
		GeometryType geometry) {

	const double dt = params->dem_dt;
	const double dem_diameter = params->dem_diameter;
	const double dem_k = params->dem_k;
	const double dem_gamma = params->dem_gamma;
	const double dem_mass = params->dem_mass;

	std::for_each(dem->begin(),dem->end(),[dt](DemType::Value& i) {
		REGISTER_DEM_PARTICLE(i);

		v = v0 + dt/2 * f;
	});
}

template<typename GeometryType>
void integrate_dem(const double dt,ptr<DemType> dem,
		ptr<Params> params,
		GeometryType geometry) {
	const int num_it = floor(dt/params->dem_dt);
	for (int i = 0; i < num_it; ++i) {
		dem_start(dem,params,geometry);
		dem_end(dem,params,geometry);
	}
	const double save_dem_dt = params->dem_dt;
	params->dem_dt = dt - num_it*params->dem_dt;
	if (params->dem_dt>0) {
		dem_start(dem,params,geometry);
		dem_end(dem,params,geometry);
	}
	params->dem_dt = save_dem_dt;
}

template<typename SphGeometryType,typename DemGeometryType>
void sphdem(ptr<SphType> sph,ptr<DemType> dem,
		ptr<Params> params,
		SphGeometryType sph_geometry,DemGeometryType dem_geometry) {

	const double dt = params->sph_dt;
	const double sph_mass = params->sph_mass;
	const double sph_prb = params->sph_prb;
	const double sph_refd = params->sph_refd;
	const double sph_gamma = params->sph_gamma;
	const double sph_visc = params->sph_visc;



	/*
	 * 0 -> 1/2 step
	 */
	sph->update_positions(sph->begin(),sph->end(),[dt,sph_mass](SphType::Value& i) {
		REGISTER_SPH_PARTICLE(i);

		if (!fixed) v += dt/2 * f;
		//h = pow(sph_mass/rho,1.0/NDIM);
		return r + dt/2 * v;
	});



//	std::for_each(sph->begin(),sph->end(),[dem](SphType::Value& i) {
//		REGISTER_SPH_PARTICLE(i);
//
//		e = 0;
//		for (auto tpl: i.get_neighbours(dem)) {
//			REGISTER_NEIGHBOUR_DEM_PARTICLE(tpl);
//			const double r = dx.norm();
//
//			e += dem_vol*W(r/h,h);
//		}
//	});
//
//	std::for_each(dem->begin(),dem->end(),[sph](DemType::Value& i) {
//		REGISTER_SPH_PARTICLE(i);
//
//		for (auto tpl: i.get_neighbours(sph)) {
//			REGISTER_NEIGHBOUR_SPH_PARTICLE(tpl);
//			const double r = dx.norm();
//
//			//interpolate needed sph values here
//		}
//	});
//
//
//	integrate_dem(dt,dem,params,dem_geometry);
//
//
//	std::for_each(sph->begin(),sph->end(),[dem](SphType::Value& i) {
//		REGISTER_SPH_PARTICLE(i);
//
//		e = 0;
//
//		for (auto tpl: i.get_neighbours(dem)) {
//			REGISTER_NEIGHBOUR_DEM_PARTICLE(tpl);
//
//			const double r = dx.norm();
//			e += dem_vol*W(r/h,h);
//		}
//	});
//
//	/*
//	 * Calculate omega
//	 */
//	std::for_each(sph->begin(),sph->end(),[&sph_geometry,sph_mass](SphType::Value& i) {
//		REGISTER_SPH_PARTICLE(i);
//
//		omega << 0,0,0;
//		for (auto tpl: i.get_neighbours(sph)) {
//			REGISTER_NEIGHBOUR_SPH_PARTICLE(tpl);
//
//			const double r = dx.norm();
//			omega -= sph_mass*pow(r,2)*F(r/h,h);
//		}
//		omega *= 1.0/(rho*NDIM);
//	});

	/*
	 * Calculate change in density
	 */
	std::for_each(sph->begin(),sph->end(),[sph,&sph_geometry,sph_mass](SphType::Value& i) {
		REGISTER_SPH_PARTICLE(i);

		dddt = 0;
		for (auto tpl: i.get_neighbours(sph)) {
			REGISTER_NEIGHBOUR_SPH_PARTICLE(tpl);

			const double r2 = dx.squaredNorm();
			if (r2 > 4.0*h*h) continue;
			if (r2 == 0) continue;
			const double r = sqrt(r2);

			dddt += sph_mass*(v-vj).dot(dx*F(r/h,h));
		}
		dddt *= 1.0/omega;

	});

	/*
	 * 1/2 -> 1 step
	 */
	sph->update_positions(sph->begin(),sph->end(),[dt,sph_mass,sph_prb,sph_refd,sph_gamma](SphType::Value& i) {
		REGISTER_SPH_PARTICLE(i);

		rho += dt * dddt;
		//h = pow(sph_mass/rho,1.0/NDIM);
		v0 = v;
		if (!fixed) v += dt/2 * f;
		const double press = sph_prb*(pow(rho/sph_refd,sph_gamma) - 1.0);
		pdr2 = press/pow(rho,2);
		return r + dt/2 * v0;
	});

	/*
	 * acceleration calculation
	 */
	std::for_each(sph->begin(),sph->end(),[sph,&sph_geometry,sph_mass,sph_visc](SphType::Value& i) {
		REGISTER_SPH_PARTICLE(i);

		f << 0,0,0;
		f = f + sph_geometry(i);

		for (auto tpl: i.get_neighbours(sph)) {
			REGISTER_NEIGHBOUR_SPH_PARTICLE(tpl);
			const double r2 = dx.squaredNorm();
			if (r2 > 4.0*h*h) continue;
			if (r2 == 0) continue;
			const double r = sqrt(r2);
			const double fdash = F(r/h,h);
			const double fdashj = F(r/hj,hj);

			/*
			 * pressure gradient
			 */
			f += -sph_mass*((1.0/omega)*pdr2*fdash + (1.0/omegaj)*pdr2j*fdashj)*dx;

			const Vect3d dv = v-vj;
			const double visc = sph_visc*(rho+rhoj)/(rho*rhoj);

			/*
			 * viscosity (morris)
			 */
			f += dv*sph_mass*visc*fdash;
		}

	});


	/*
	 * 1/2 -> 1 step for velocity
	 */
	std::for_each(sph->begin(),sph->end(),[dt](SphType::Value& i) {
		REGISTER_SPH_PARTICLE(i);

		if (!fixed) v = v0 + dt/2 * f;
	});

}


#endif /* SPHDEM_H_ */
