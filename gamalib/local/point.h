/*  
    Geodesy and Mapping C++ Library (GNU GaMa / GaMaLib)
    Copyright (C) 2000  Ales Cepek <cepek@fsv.cvut.cz>

    This file is part of the GNU GaMa / GaMaLib C++ Library.
    
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

/*
 *  $Id: point.h,v 1.1 2001/12/07 12:38:37 cepek Exp $
 */

#ifndef GaMaLib___class__Point_h
#define GaMaLib___class__Point_h

#include <gamalib/float.h>
#include <gamalib/exception.h>
#include <gamalib/language.h>

namespace GaMaLib {


class Point {
public:

  struct XYZ 
  { 
    Double x, y, z; 
    XYZ(Double px, Double py, Double pz) : x(px), y(py), z(pz) {}
  };
  struct XY  
  { 
    Double x, y;    
    XY(Double px, Double py) : x(px), y(py) {}
  };
  struct ZZ   
  { 
    Double z; 
    ZZ(Double pz) : z(pz) {}
  }; 


  Point()      : bxy_(false), bz_(false), pst_(unused_) {}
  Point(XYZ p) : x_(p.x), y_(p.y), z_(p.z), bxy_(true), bz_(true),  
                 pst_(unused_) {}
  Point(XY  p) : x_(p.x), y_(p.y), bxy_(true), bz_(false), pst_(unused_) {}
  Point(ZZ  p) : z_(p.z), bxy_(false), bz_(true), pst_(unused_) {}
  

  Double y() const 
    { 
      if (!bxy_) throw GaMaLib::Exception(T_POBS_bad_data);
      return y_; 
    }
  Double x() const 
    { 
      if (!bxy_) throw GaMaLib::Exception(T_POBS_bad_data); 
      return x_; 
    }
  Double z() const 
    { 
      if (!bz_ ) throw GaMaLib::Exception(T_POBS_bad_data); 
      return z_; 
    }

  void set_xy  (Double x, Double y) { bxy_ = true; x_ = x; y_ = y; }
  void set_z   (Double z)           { bz_  = true; z_ = z; }
  void unset_xy() { bxy_ = false; }
  void unset_z () { bz_  = false; }
  bool test_xy () const { return bxy_; }
  bool test_z  () const { return bz_;  }

  int& index_y() { return iy_; }         // indexes in project equations
  int& index_x() { return ix_; }
  int& index_z() { return iz_; }
  int index_y() const { return iy_; }
  int index_x() const { return ix_; }
  int index_z() const { return iz_; }
  
  void set_fixed_xy()       { pst_ &= ~active_xy_; pst_ |= xy_fixed_;     }
  void set_free_xy()        { pst_ &= ~active_xy_; pst_ |= xy_adjusted_;  }
  void set_constrained_xy() { pst_ &= ~active_xy_; 
                              pst_ |= (xy_adjusted_ | xy_constrained_);   }
  void set_fixed_z()        { pst_ &= ~active_z_;  pst_ |= z_fixed_;      }
  void set_free_z()         { pst_ &= ~active_z_;  pst_ |= z_adjusted_;   }
  void set_constrained_z()  { pst_ &= ~active_z_;
                              pst_ |= (z_adjusted_  | z_constrained_);    }
  void unused_xy()          { pst_ &= ~active_xy_;                        }
  void unused_z()           { pst_ &= ~active_z_;                         }
  void unused()             { pst_ = unused_;                             }
  
  bool active_xy()      const { return pst_ & active_xy_;      }
  bool fixed_xy()       const { return pst_ & xy_fixed_;       }
  bool free_xy()        const { return pst_ & xy_adjusted_;    }
  bool constrained_xy() const { return pst_ & xy_constrained_; }
  bool fixed_z()        const { return pst_ & z_fixed_;        }
  bool free_z()         const { return pst_ & z_adjusted_;     }
  bool constrained_z()  const { return pst_ & z_constrained_;  }
  bool active_z()       const { return pst_ & active_z_;       }
  bool active()         const { return pst_ & active_;         }

  /*
  void set_fixed_point()        { pst_ &= ~active_xy_; pst_ |= xy_fixed_;     }
  void set_free_point()         { pst_ &= ~active_xy_; pst_ |= xy_adjusted_;  }
  void set_constrained_point()  { pst_ &= ~active_xy_; 
                                  pst_ |= (xy_adjusted_ | xy_constrained_);   }
  void set_fixed_height()       { pst_ &= ~active_z_;  pst_ |= z_fixed_;      }
  void set_free_height()        { pst_ &= ~active_z_;  pst_ |= z_adjusted_;   }
  void set_constrained_height() { pst_ &= ~active_z_;  pst_ |= z_constrained_;}
  void unset_network_point()    { pst_ &= ~active_xy_;                        }
  void unset_height_point()     { pst_ &= ~active_z_;                         }
  void unused_point()           { pst_ = unused_;                             }
  
  bool network_point()      const { return pst_ & active_xy_;      }
  bool fixed_point()        const { return pst_ & xy_fixed_;       }
  bool free_point()         const { return pst_ & xy_adjusted_;    }
  bool constrained_point()  const { return pst_ & xy_constrained_; }
  bool fixed_height()       const { return pst_ & z_fixed_;        }
  bool free_height()        const { return pst_ & z_adjusted_;     }
  bool constrained_height() const { return pst_ & z_constrained_;  }
  bool height_point()       const { return pst_ & active_z_;       }
  */
private:

  Double  x_, y_, z_;       // coordinates
  bool    bxy_, bz_;        // coordinates are/are_not defined
  int     ix_, iy_, iz_;    // indexes of unknowns in project equations
  int     pst_;             // point status in adjustment

  enum 
  { 
    unused_          =  0, 
    xy_fixed_        =  1, 
    xy_adjusted_     =  2, 
    xy_constrained_  =  4,
    z_fixed_         =  8,
    z_adjusted_      = 16,
    z_constrained_   = 32,
    active_xy_   = xy_fixed_  | xy_adjusted_ | xy_constrained_,
    active_z_    = z_fixed_   | z_adjusted_  | z_constrained_,
    active_      = active_xy_ | active_z_
  };

};

}   // namespace GaMaLib;

#endif





