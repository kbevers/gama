/*  
   GNU Gama -- adjustment of geodetic networks
   Copyright (C) 2003  Ales Cepek <cepek@fsv.cvut.cz>

   This file is part of the GNU Gama C++ library.
   
   This library is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this library; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/* $Id: g3_parameter.h,v 1.16 2003/12/28 16:42:34 uid66336 Exp $  */

#include <cstddef>
#include <gnu_gama/model.h>
#include <gnu_gama/list.h>


#ifndef GNU_gama_____g3______parameter______h__________GNUgamag3parameterh
#define GNU_gama_____g3______parameter______h__________GNUgamag3parameterh


namespace GNU_gama { namespace g3 {


  class Parameter {
  public:
    
    Parameter() : val(0), cor(0), dif(0.05) {}
    // ~Parameter() {}
    
    double value() const { return val + cor; }

    double init_value() const { return val; }
    double correction() const { return cor; }
    double step_size () const { return dif; } 
    std::size_t index() const { return ind; }

    void set_init_value(double p) { val = p; cor = 0; }
    void set_correction(double p) { cor = p; }
    void set_step_size (double p) { dif = p; }
    void set_index(std::size_t t) { ind = t; }

    void set_unused() { state_ = unused_; }
    void set_fixed () { state_ = fixed_;  }
    void set_free  () { state_ = free_;   }
    void set_constr() { state_ = constr_; }

    bool active() const { return state_ != unused_; }
    bool unused() const { return state_ == unused_; }
    bool fixed () const { return state_ == fixed_;  }
    bool free  () const { return state_ &  free_;   }
    bool constr() const { return state_ == constr_; }

  private:
    
    double val;
    double cor;
    double dif;
    std::size_t ind;

    enum 
      {
        unused_ = 0,
        fixed_  = 1,
        free_   = 2,
        constr_ = 4 + free_
      } state_;

  };
  
}}
  
#endif
