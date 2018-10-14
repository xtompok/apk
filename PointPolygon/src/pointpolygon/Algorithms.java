/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pointpolygon;

import java.awt.Point;
import java.awt.Polygon;
import static java.lang.Math.abs;
import static java.lang.Math.acos;
import static java.lang.Math.sqrt;

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
    
    public static double angle(double ux, double uy, double vx, double vy){
        double prod = dotProd(ux, uy, vx, vy);
        double ulen = len(ux,uy);
        double vlen = len(vx,vy);
        return acos(prod / (ulen * vlen));
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
    
    
}
