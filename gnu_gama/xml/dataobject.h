/*  
    Geodesy and Mapping C++ Library (GNU GaMa / GaMaLib)
    Copyright (C) 2002  Ales Cepek <cepek@fsv.cvut.cz>

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
 *  $Id: dataobject.h,v 1.4 2003/11/25 22:17:15 cepek Exp $
 */

#ifndef GaMaLib_GaMa_XML_Data_Object__object___h_
#define GaMaLib_GaMa_XML_Data_Object__object___h_

#include <string>
#include <sstream>
#include <gnu_gama/adj/adj.h>
#include <gnu_gama/g3/g3_model.h>

namespace GNU_gama { namespace DataObject {

  class Base {
  public:

    virtual ~Base() 
      {
      }
    virtual std::string xml() const = 0;

    static  std::string xml_begin();
    static  std::string xml_end();
  };


  class Text : public Base {
  public:
  
    std::string text;
  
    Text() 
      {
      }    
    Text(std::string s) : text(s) 
      {
      }    
    std::string xml() const 
      {
        if (!text.empty())
          return "\n<text>" + text + "</text>\n";

        return "";
      }
  };


  class AdjInput : public Base {
  public:
  
    GNU_gama::AdjInputData *data;
  
    AdjInput() : data(0)
      {
      }    
    AdjInput(GNU_gama::AdjInputData *d) : data(d)
      {
      }    
    std::string xml() const 
      {
        if (data) 
          {
            std::stringstream out;
            data->write_xml(out);
            
            return out.str();
          }

        return "";
      }
  };


  class g3_model : public Base {
  public:
  
    GNU_gama::g3::g3_Model *model;
  
    g3_model() : model(0)
      {
      }    
    g3_model(GNU_gama::g3::g3_Model *m) : model(m)
      {
      }    
    std::string xml() const 
      {
        if (model) 
          {
            std::stringstream out;
            model->write_xml(out);
            
            return out.str();
          }

        return "";
      }
  };

}}       // namespace GNU_gama::DataObject


#endif
















