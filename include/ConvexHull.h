/*
 **************************************************************************
 * Class: Convex Hull                                                     *
 * By Arash Partow - 2001                                                 *
 * URL: http://www.partow.net                                             *
 *                                                                        *
 * Copyright Notice:                                                      *
 * Free use of this library is permitted under the guidelines and         *
 * in accordance with the most current version of the Common Public       *
 * License.                                                               *
 * http://www.opensource.org/licenses/cpl.php                             *
 *                                                                        *
 **************************************************************************
*/


#ifndef INCLUDE_CONVEXHULL_H
#define INCLUDE_CONVEXHULL_H

#include <vector>


struct point2d
{
   point2d(double _x = 0.0 , double _y = 0.0) : x(_x), y(_y){}
   double x;
   double y;
};

class ConvexHull
{
   public:

     virtual ~ConvexHull(){};
     virtual bool operator()(const std::vector<point2d>& pnt, std::vector<point2d>& final_hull) = 0;

};


#endif
