/*
 * BucketSort.h
 *
 * Copyright 2012 Martin Robinson
 *
 * This file is part of RD_3D.
 *
 * RD_3D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * RD_3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with RD_3D.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 22 Oct 2012
 *      Author: robinsonm
 */

#ifndef BUCKETSORT_H_
#define BUCKETSORT_H_

#include <boost/iterator/iterator_facade.hpp>
#include "Vector.h"
#include "Constants.h"
#include "Log.h"
#include <vector>
#include <iostream>


namespace Aboria {

const int CELL_EMPTY = -1;


class BucketSort {
public:
	class const_iterator
	  : public boost::iterator_facade<
	        const_iterator
	      , const int
	      , boost::forward_traversal_tag
	    >
	{
	 public:
		  const_iterator()
	      : m_node(),my_index(-1),self(false),radius2(0) {
			  cell_empty.push_back(CELL_EMPTY);
			  m_node = cell_empty.begin();
		  }

	    void go_to_next_candidate() {
	    	m_node = bucket_sort->linked_list.begin() + *m_node;
	    	if (*m_node == CELL_EMPTY) {
	    		if (surrounding_cell_offset_i == surrounding_cell_offset_end) return;
	    		surrounding_cell_offset_i++;
	    	}
	    	if (surrounding_cell_offset_i != surrounding_cell_offset_end) {
	    		m_node = cell_i + *surrounding_cell_offset_i;
	    	} else if (self) {
	    		m_node = cell_i;
	    		while (*m_node != my_index) {
	    	    	m_node = bucket_sort->linked_list.begin() + *m_node;
	    		}
	    	}
	    }

	    explicit const_iterator(const BucketSort* bucket_sort,
	    		const Vect3d centre, const double radius,
	    		const int my_index = -1, const bool self = false)
	    : bucket_sort(bucket_sort),
	      cell_i(bucket_sort->cells.begin() + bucket_sort->find_cell_index(centre)),
	      my_index(my_index),self(self),
	      centre(centre),radius2(radius*radius) {
	    	surrounding_cell_offset_i = bucket_sort->surrounding_cell_offsets.begin();
	    	if (self) {
	    		surrounding_cell_offset_end = surrounding_cell_offset_i
	    						+ (bucket_sort->surrounding_cell_offsets.size()-1)/2;
	    	} else {
	    		surrounding_cell_offset_end = bucket_sort->surrounding_cell_offsets.end();
	    	}
	    	m_node = cell_i + *surrounding_cell_offset_i;
	    	while ((centre-bucket_sort->correct_position_for_periodicity(centre, bucket_sort->positions[*m_node])).squaredNorm()
	    			> radius2) {
	    		go_to_next_candidate();
	    	}
	    }

	 private:
	    friend class boost::iterator_core_access;

	    bool equal(const_iterator const& other) const {
	        return *m_node == *(other.m_node);
	    }

	    void increment() {
	    	while ((centre-bucket_sort->correct_position_for_periodicity(centre, bucket_sort->positions[*m_node])).squaredNorm()
	    			> radius2) {
	    		go_to_next_candidate();
	    	}
	    }


	    const int& dereference() const
	    { return *m_node; }

	    const BucketSort* bucket_sort;
	    std::vector<int>::const_iterator m_node;
	    std::vector<int>::const_iterator cell_i;
	    //Value* const linked_list;
	    const int my_index;
	    const bool self;
	    const Vect3d centre;
	    const double radius2;
	    std::vector<int> cell_empty;
	//    std::vector<Vect3d>::const_iterator positions;
	//    std::vector<Value>::const_iterator linked_list;
	    std::vector<int>::const_iterator surrounding_cell_offset_i,surrounding_cell_offset_end;
	};

	BucketSort(Vect3d low, Vect3d high, Vect3b periodic):
		low(low),high(high),domain_size(high-low),periodic(periodic),
		empty_cell(CELL_EMPTY) {
		LOG(2,"Creating bucketsort data structure with lower corner = "<<low<<" and upper corner = "<<high);
		const double dx = (high-low).maxCoeff()/10.0;
		reset(low, high, dx);
	}

	void reset(const Vect3d& low, const Vect3d& high, double _max_interaction_radius);
	inline const Vect3d& get_low() {return low;}
	inline const Vect3d& get_high() {return high;}

	template<typename T, typename F>
	void embed_points(const T& begin_iterator, const T& end_iterator, const F& return_vect3d);
	const_iterator find_broadphase_neighbours(const Vect3d& r, const double radius, const int my_index, const bool self) const;
	const_iterator end() const;
	Vect3d correct_position_for_periodicity(const Vect3d& source_r, const Vect3d& to_correct_r) const;
	Vect3d correct_position_for_periodicity(const Vect3d& to_correct_r) const;

private:
	inline int vect_to_index(const Vect3i& vect) const {
		return vect[0] * num_cells_along_yz + vect[1] * num_cells_along_axes[1] + vect[2];
	}
	inline int find_cell_index(const Vect3d &r) const {
		const Vect3i celli = ((r-low).cwiseProduct(inv_cell_size) + Vect3d(1.0,1.0,1.0)).cast<int>();
		ASSERT((celli[0] >= 0) && (celli[0] < num_cells_along_axes[0]), "position is outside of x-range");
		ASSERT((celli[1] >= 0) && (celli[1] < num_cells_along_axes[1]), "position is outside of y-range");
		ASSERT((celli[2] >= 0) && (celli[2] < num_cells_along_axes[2]), "position is outside of z-range");
		//std::cout << " looking in cell " << celli <<" out of total cells " << num_cells_along_axes << " at position " << r<< std::endl;
		return vect_to_index(celli);
	}

	std::vector<Vect3d> positions;
    std::vector<int> cells;
    std::vector<std::vector<int> > ghosting_indices_pb;
    std::vector<std::pair<int,int> > ghosting_indices_cb;
    std::vector<int> dirty_cells;
	std::vector<int> linked_list;
	std::vector<int> neighbr_list;
	Vect3d low,high,domain_size;
	const Vect3b periodic;
	Vect3d cell_size,inv_cell_size;
	Vect3i num_cells_along_axes;
	int num_cells_along_yz;
	double max_interaction_radius;
	std::vector<int> surrounding_cell_offsets;
	const int empty_cell;
};

template<typename T, typename F>
void BucketSort::embed_points(const T& begin_iterator, const T& end_iterator, const F& return_vect3d) {
	const unsigned int n = std::distance(begin_iterator,end_iterator);
	linked_list.assign(n, CELL_EMPTY);
	positions.resize(n);
	const bool particle_based = dirty_cells.size() < cells.size();
	if (particle_based) {
		for (int i: dirty_cells) {
			cells[i] = CELL_EMPTY;
		}
	} else {
		cells.assign(cells.size(), CELL_EMPTY);
	}

	dirty_cells.clear();
	int i = 0;
	for (auto it = begin_iterator; it != end_iterator; ++it, ++i) {
		positions[i] = return_vect3d(it);
		const int celli = find_cell_index(positions[i]);
		const int cell_entry = cells[celli];

		// Insert into own cell
		cells[celli] = i;
		dirty_cells.push_back(celli);
		linked_list[i] = cell_entry;

		// Insert into ghosted cells
		if (particle_based) {
			for (int j: ghosting_indices_pb[celli]) {
				const int cell_entry = cells[j];
				cells[j] = i;
				dirty_cells.push_back(j);
				linked_list[i] = cell_entry;
			}
		}
	}

	if (!particle_based) {
		for (std::vector<std::pair<int,int> >::iterator index_pair = ghosting_indices_cb.begin(); index_pair != ghosting_indices_cb.end(); ++index_pair) {
			//BOOST_FOREACH(std::pair<int,int> index_pair, ghosting_indices) {
			cells[index_pair->first] = cells[index_pair->second];
		}
	}
}




}

#endif /* BUCKETSORT_H_ */