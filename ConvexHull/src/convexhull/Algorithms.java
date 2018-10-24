/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package convexhull;

import java.awt.List;
import java.awt.Point;
import java.awt.Polygon;
import java.awt.geom.Path2D;
import java.awt.geom.Point2D;
import static java.lang.Math.abs;
import static java.lang.Math.acos;
import static java.lang.Math.atan2;
import static java.lang.Math.sqrt;
import java.util.LinkedList;

/**
 *
 * @author jethro
 */
public class Algorithms {
    public enum PositionEnum {INSIDE,OUTSIDE,BOUNDARY}
    
    public static double dotProd(double ux, double uy, double vx, double vy){
        return ux*vx + uy*vy;
    }
    
    public static double len(double ux, double uy){
        return sqrt(ux*ux + uy*uy);
    }
    
    public static double scalarProjection(double ux, double uy, double vx, double vy){
        double prod = dotProd(ux, uy, vx, vy);
        double ulen = len(ux,uy);
        double vlen = len(vx,vy);
        return prod/(ulen*vlen);
    }
    
    public static double angle(double ux, double uy, double vx, double vy){
        return acos(scalarProjection(ux,uy,vx,vy));
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
        hullPoints.add(new Point2D.Double(0,0));
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
                        
            double max;
            Point2D maxPt;
            max = Double.MAX_VALUE;
            maxPt = null;
            
            for (Point2D pt: points){
                if ((pt == curPt)||(pt == prevPt)){
                    continue;
                }
                double vx,vy;
                vx = pt.getX() - curPt.getX();
                vy = pt.getY() - curPt.getY();
                double m = scalarProjection(ux,uy,vx,vy);
                if (m < max){
                    max = m;
                    maxPt = pt;
                }   
            }
            
            hullPoints.add(maxPt);
            prevPt = curPt;
            curPt = maxPt;
            
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
    
    
}
