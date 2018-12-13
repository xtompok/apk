/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package setoperations;

import java.awt.Point;
import java.awt.geom.Point2D;
import static java.lang.Math.abs;
import static java.lang.Math.acos;
import static java.lang.Math.sqrt;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

/**
 *
 * @author mazurd
 */
public class Algorithms {
    
    public static IntersectionPoint calcIntersection(Edge e1, Edge e2){
        double ux,uy,vx,vy,wx,wy;
        ux = e1.end.getX() - e1.start.getX();
        uy = e1.end.getY() - e1.start.getY();
        vx = e2.end.getX() - e2.start.getX();
        vy = e2.end.getY() - e2.start.getY();
        wx = e1.start.getX() - e2.start.getX();
        wy = e1.start.getY() - e2.start.getY();
        
        double k1,k2,k3;
        k1 = vx*wy - vy*wx;
        k2 = ux*wy - uy*wx;
        k3 = vy*ux - vx*uy;
        
        double alpha, beta;
        alpha = k1/k3;
        beta = k2/k3;
        
        if (alpha > 1 || beta > 1 || alpha < 0 || beta < 0){
            return null;
        }
        
        Point2D pt;
        pt = new Point2D.Double(e1.start.getX()+alpha*ux,e1.start.getY()+alpha*uy);
           
        
        return new IntersectionPoint(pt,e1,e2,alpha,beta);   
    }
    
    public static List<IntersectionPoint> allIntersections(Polygon polyA, Polygon polyB){
        List<IntersectionPoint> ints;
        ints = new LinkedList<>();
        
        for (Edge ea : polyA.edges){
            for (Edge eb : polyB.edges){
                IntersectionPoint pt;
                pt = calcIntersection(ea, eb);
                if (pt != null){
                    ints.add(pt);
                }
            }
        }
        
        return ints;
    }
    
    
    public static List<Edge> divideEdge(Edge e, List<IntersectionPoint> points){
        List<IntersectionPoint> myInts;
        myInts = new LinkedList<>();
        
        for (IntersectionPoint pt : points){
            if (pt.e1 == e || pt.e2 == e){
                myInts.add(pt);
            }
        }
        
        Collections.sort(myInts,(IntersectionPoint t, IntersectionPoint t1) -> {
            double c1,c2;
            if (t.e1 == e){
                c1 = t.alpha;
            } else {
                c1 = t.beta;
            }
            if (t1.e1 == e){
                c2 = t1.alpha;
            } else {
                c2 = t1.beta;
            }
            if (c1 < c2){
                return -1;
            } else if (c1 > c2){
                return 1;
            }
            return 0;
        });
        
        List<Edge> divided;
        divided = new LinkedList<>();
        
        if (myInts.size() == 0){
            divided.add(e);
            return divided;
        }
        
        divided.add(new Edge(e.start,myInts.get(0).point));
        for (int i=0;i<myInts.size()-1;i++){
            divided.add(new Edge(myInts.get(i).point,myInts.get(i+1).point));
        }
        divided.add(new Edge(myInts.get(myInts.size()-1).point,e.end));
        
        return divided;  
    }
    
    public static Polygon divideAll(Polygon poly,List<IntersectionPoint> points){
        Polygon out;
        out = new Polygon();
        for (Edge e: poly.edges){
            out.edges.addAll(divideEdge(e, points));
        }
        return out;
    }
    
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
    
    public static PositionEnum pointPolygonWinding(Point2D pt, Polygon poly){
        double sumAngle = 0;
        final double eps = 0.01;
        
        for (Edge e: poly.edges){
            double ux = e.start.getX() - pt.getX();
            double uy = e.start.getY() - pt.getY();
            
            double vx = e.end.getX() - pt.getX();
            double vy = e.end.getY() - pt.getY();
            
            sumAngle += angle(ux, uy, vx, vy);
        
        }
        
        if (abs(abs(sumAngle) - 2*Math.PI) < eps){
            return PositionEnum.INSIDE;
        }
        return PositionEnum.OUTSIDE;       
    }
    
    public static void setInside(Edge e, Polygon poly){
        Point2D middle;
        middle = new Point2D.Double(
                (e.start.getX()+e.end.getX())/2,
                (e.start.getY()+e.end.getY())/2);
        
        if (pointPolygonWinding(middle, poly)==PositionEnum.INSIDE){
            e.inside = true;
        } else {
            e.inside = false;
        }
    }
    
    public static void setInside(Polygon polyA, Polygon polyB){
        for (Edge e: polyA.edges){
            setInside(e, polyB);
        }
        for (Edge e: polyB.edges){
            setInside(e, polyA);
        }
    }
    
    
}
