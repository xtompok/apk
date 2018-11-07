/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package convexhull;

import java.util.List;
import java.awt.Point;
import java.awt.Polygon;
import java.awt.geom.Path2D;
import java.awt.geom.Point2D;
import static java.lang.Math.abs;
import static java.lang.Math.acos;
import static java.lang.Math.atan2;
import static java.lang.Math.sqrt;
import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.ListIterator;

/**
 *
 * @author jethro
 */
public class Algorithms {
    public static final double EPSILON = 0.0001;
    
    public enum PositionEnum {INSIDE,OUTSIDE,BOUNDARY}
    public enum OrientationEnum {CW,COLINEAR,CCW}
    
    public static OrientationEnum getOrientation(Point2D p1, Point2D p2, Point2D p3){
         double val = (p2.getY() - p1.getY()) * (p3.getX() - p2.getX()) - 
                  (p2.getX() - p1.getX()) * (p3.getY() - p2.getY()); 
         if (abs(val) < EPSILON){
             return OrientationEnum.COLINEAR;
         }else if (val > 0){
             return OrientationEnum.CW;
         }else {
             return OrientationEnum.CCW;
         }
    }
    
    public static double dotProd(double ux, double uy, double vx, double vy){
        return ux*vx + uy*vy;
    }
    
    public static double len(double ux, double uy){
        return sqrt(ux*ux + uy*uy);
    }
    
    public static double dotProdNorm(double ux, double uy, double vx, double vy){
        double prod = dotProd(ux, uy, vx, vy);
        double ulen = len(ux,uy);
        double vlen = len(vx,vy);
        return prod/(ulen*vlen);
    }
    
    public static double angle(double ux, double uy, double vx, double vy){
        return acos(dotProdNorm(ux,uy,vx,vy));
    }
    
    
    public static PositionEnum pointPolygonWinding(Point pt, Polygon poly){
        double sumAngle = 0;
        final double eps = 0.01;
        
        for (int i=0;i<poly.npoints;i++){
            double ux = poly.xpoints[i] - pt.x;
            double uy = poly.ypoints[i] - pt.y;
            
            double vx = poly.xpoints[(i+1)%poly.npoints] - pt.x;
            double vy = poly.ypoints[(i+1)%poly.npoints] - pt.y;
            
            sumAngle += angle(ux, uy, vx, vy);
        }
        
        if (abs(abs(sumAngle) - 2*Math.PI) < eps){
            return PositionEnum.INSIDE;
        }
        return PositionEnum.OUTSIDE;
        
    }
    public static Path2D jarvisScan(Point2D [] points){
        Path2D hull;
        hull = new Path2D.Double();
        
        LinkedList<Point2D> hullPoints;
        hullPoints = new LinkedList<>();
        
        Point2D miny;
        miny = points[0];
        for (Point2D pt : points){
            if (pt.getY() < miny.getY()){
                miny = pt;
            }
        }
        
        hull.moveTo(miny.getX(),miny.getY());
        hullPoints.add(new Point2D.Double(0,miny.getY()));
        hullPoints.add(miny);
        
        Point2D prevPt;
        Point2D curPt;
        curPt = miny;
        prevPt = hullPoints.getFirst();
        
        int count = 0;
        do {
            double ux;
            double uy;
            ux = prevPt.getX()- curPt.getX();
            uy = prevPt.getY()- curPt.getY();
                        
            double min;
            Point2D minPt;
            min = Double.MAX_VALUE;
            minPt = null;
            
            for (Point2D pt: points){
                if ((pt == curPt)||(pt == prevPt)){
                    continue;
                }
                double vx,vy;
                vx = curPt.getX() - pt.getX();
                vy = curPt.getY() - pt.getY();
                double m = angle(ux,uy,vx,vy);
                if (m < min){
                    min = m;
                    minPt = pt;
                }   
            }
            
            hullPoints.add(minPt);
            prevPt = curPt;
            curPt = minPt;
            
        }
        while (hullPoints.getLast()!=miny);
        
        double xs[];
        double ys[];
        xs = new double[hullPoints.size()-1];
        ys = new double[hullPoints.size()-1];
        
        for (int i=1;i<hullPoints.size();i++){
            hull.lineTo(hullPoints.get(i).getX(), hullPoints.get(i).getY());
        }
        hull.lineTo(hullPoints.get(1).getX(),hullPoints.get(1).getY());
    
        return hull;
    }
    
    public static Path2D sweepHull(Point2D [] points){
        Path2D hull;
        hull = new Path2D.Double();
        Arrays.sort(points, (Point2D p1, Point2D p2) -> Double.compare(p1.getX(), p2.getX()));
        
        List<Point2D> upperHull;
        List<Point2D> lowerHull;
        
        upperHull = new LinkedList<>();
        lowerHull = new LinkedList<>();
        
        for (Point2D pt: points){
            upperHull.add(pt);
            lowerHull.add(pt);
        }
        
        fixConvexity(upperHull, OrientationEnum.CW);
        fixConvexity(lowerHull, OrientationEnum.CCW);
        
        hull.moveTo(lowerHull.get(0).getX(), lowerHull.get(0).getY());
        
        for (Point2D pt : lowerHull){
            hull.lineTo(pt.getX(), pt.getY());
        }
        for (Point2D pt : upperHull){
            hull.lineTo(pt.getX(), pt.getY());
        }
        
        return hull;
    }
    
    public static void fixConvexity(List<Point2D> points,OrientationEnum orientation){
        if (points.size() < 3){
            return;
        }
        ListIterator<Point2D> iterator;
        iterator = points.listIterator();
        Point2D prev;
        Point2D cur;
        Point2D next;
        prev = iterator.next();
        cur = iterator.next();
        next = iterator.next();

        
        while (iterator.hasNext()){
            System.err.format("p:%h, c:%h, n:%h\n", prev,cur,next);
            // Orientace je v pořádku
            if (getOrientation(prev,cur,next) == orientation){
                prev = cur;
                cur = next;
                next = iterator.next();
                continue;
            }
            System.out.println("Orientation test failed");
            //Orientace je rozbitá => smažeme předchozí vrchol
            iterator.previous();
            iterator.remove();
            cur = iterator.next();
            next = iterator.next();
        }
        
    }
    
    public static Path2D quickHull(Point2D [] points){
        Path2D hull;
        hull = new Path2D.Double();
        
        Point2D minx;
        Point2D maxx;
        minx = points[0];
        maxx = points[0];
        for (Point2D pt : points){
            if (pt.getX() < minx.getX()){
                minx = pt;
            }
            if (pt.getX() > maxx.getX()){
                maxx = pt;
            }
        }
        
        hull.moveTo(maxx.getX(),maxx.getY());
        
    
        
        List<Point2D>[] splitted = splitPoints(minx, maxx, Arrays.asList(points));
        Point2D[] upperHull = qh(maxx,minx,splitted[0]);
        Point2D[] lowerHull = qh(minx,maxx,splitted[1]);
        
        System.out.println(upperHull.length);
        System.out.println(lowerHull.length);
        
        for (Point2D pt: lowerHull){
            hull.lineTo(pt.getX(), pt.getY());
        }

        hull.lineTo(minx.getX(), minx.getY());
        for (Point2D pt: upperHull){
            hull.lineTo(pt.getX(), pt.getY());
        }

        hull.lineTo(maxx.getX(),maxx.getY());

        return hull;
    }
    
    private static Point2D[] qh(Point2D start, Point2D end, List<Point2D> points){
        if (points.size() == 0){
            return new Point2D[0];
        }
        
        double max = -1;
        Point2D farthestPt = null;
        
        for (Point2D pt: points){
            double dist = distanceFromLine(start, end, pt);
            if (dist > max){
                max = dist;
                farthestPt = pt;
            }
        }
        
        if (farthestPt == null){
            System.out.println("FP je null");
        }
        
        List<Point2D>[] startPts = splitPoints(farthestPt,start,points);
        List<Point2D>[] endPts = splitPoints(end,farthestPt, points);
        if (startPts[0].size() + startPts[1].size() + 2 != points.size()){
            System.out.print("Nesedi soucet: ");
            System.out.format("0: %d 1: %d p:%d\n",startPts[0].size(),startPts[1].size(),points.size());
        }
        
        
        Point2D[] startHull = qh(start,farthestPt,startPts[0]);
        Point2D[] endHull = qh(farthestPt,end,endPts[0]);
        
        Point2D[] res;
        res = new Point2D[startHull.length+endHull.length+1];
        
        int residx = 0;
        for (Point2D pt: endHull){
            res[residx] = pt;
            residx++;
        }
        res[residx] = farthestPt;
        residx++;
        for (Point2D pt: startHull){
            res[residx] = pt;
            residx++;
        }
        
        return res;
        
      
    
    }
    
    public static List<Point2D> [] splitPoints(Point2D start, Point2D end, List<Point2D> points){
        double epsilon = 0.0001;
        
        List<Point2D>[] res = new List[2];
        res[0] = new LinkedList<>();
        res[1] = new LinkedList<>();
        
        for (Point2D pt: points){
            if (pt == start || pt == end){
                continue;
            }
            double det = (end.getX()-start.getX())*(pt.getY()-start.getY()) - 
                    (end.getY()-start.getY())*(pt.getX() - start.getX());
            if (det > epsilon){
                res[0].add(pt);
            }else if (det < -epsilon){
                res[1].add(pt);
            }
        }
        
        return res;
        
    
    }
    
    public static double distanceFromLine(Point2D start, Point2D end, Point2D pt){
        double nom;
        double denom;
        
        nom = abs((end.getY()-start.getY())*pt.getX() - 
                  (end.getX()-start.getX())*pt.getY() +
                  end.getX()*start.getY() - 
                  end.getY()*start.getX());
        denom = sqrt((end.getY()-start.getY())*
                     (end.getY()-start.getY())+
                     (end.getX()-start.getX())*
                     (end.getX()-start.getX()));
        return nom/denom;
    }
    
    
}
