/*  
    C++ Matrix/Vector templates (gMatVec 0.9.14)
    Copyright (C) 1999  Ales Cepek <cepek@fsv.cvut.cz>

    This file is part of the gMatVec C++ Matrix/Vector template library.
    
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
 *  $Id: memrep.h,v 1.1 2001/12/07 11:59:46 cepek Exp $
 */

#ifndef gMatVec_MemRep__h_
#define gMatVec_MemRep__h_

#include <cstring>
#include <gmatvec/inderr.h>

namespace gMatVec {
  
template <class Float=double, class Exc=Exception>
class MemRep {

private:

   struct Mrep {
     Float* m;
     Index  sz;
     Index  n;
     
     ~Mrep() { delete[] m; }
     Mrep() : m(0), sz(0), n(1) {}
     Mrep(Index nsz) : m(new Float[nsz]), sz(nsz), n(1) {}
     Mrep(Index nsz, const Float* p) : m(new Float[nsz]), sz(nsz), n(1)
       {
         using namespace std;
         memcpy(m, p, sz*sizeof(Float));
       }
     
     private:
     Mrep(const Mrep&);
     Mrep& operator=(const Mrep&);
   };      /* struct Mrep; */
  
   mutable Mrep* rep;

protected:

   MemRep() { rep = new Mrep; }
   MemRep(Index nsz)
      {
         if (nsz>0)
            rep = new Mrep(nsz);
         else if (nsz == 0)
            rep = new Mrep;
         else
            throw Exc(BadRank, "MemRep::MemRep(Index nsz)");
      }
   MemRep(const MemRep& x) { x.rep->n++; rep = x.rep; }
   MemRep& operator=(const MemRep& x)
      {
         x.rep->n++;                      // protect against "x = x"
         if (--rep->n == 0) delete rep;
         rep = x.rep;
         return *this;
      }
   ~MemRep() { if (--rep->n) return; delete rep; }

   void resize(Index nsz)
      {
         if (nsz == rep->sz)
            return;
         else
            {
               if (--rep->n == 0) delete rep;
               if (nsz > 0)
                  rep = new Mrep(nsz);
               else
                  rep = new Mrep;
            }
      }

   Index size() const { return rep->sz; }

public:

   typedef Float* iterator;
   iterator begin() 
     {
         if (rep->n > 1) { --rep->n; rep = new Mrep(rep->sz, rep->m); }
         return rep->m;       
     }
   iterator end() { return begin() + rep->sz; }

   typedef const Float* const_iterator;
   const_iterator begin() const { return rep->m; }
   const_iterator end()   const { return rep->m + rep->sz; }

};      /* class MemRep; */


}   // namespace gMatVec

#endif


