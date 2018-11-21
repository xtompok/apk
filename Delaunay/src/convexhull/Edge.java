/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package convexhull;

import java.awt.geom.Point2D;
import java.util.Objects;
import java.util.Set;

/**
 *
 * @author jethro
 */
public class Edge {
    public Point2D p1;
    public Point2D p2;
    
   // private static Set<Point2D[]> edges;
    
   /* public static Edge createEdge(Point2D p1, Point2D p2){
        Point2D [] pts = new 
    }*/
    
    public Edge(Point2D p1, Point2D p2){
        this.p1 = p1;
        this.p2 = p2;
    }
    
    public void swap(){
        Point2D tmp;
        tmp = p1;
        p1 = p2;
        p2 = tmp;
    }
    @Override
    public boolean equals(Object obj){
        if (!(obj instanceof Edge)){
            return false;
        }
        Edge e = (Edge)obj;
        if ((e.p1 == p1) && (e.p2 == p2)){
            return true;
        }
        return false;
    }
    
    @Override
    public int hashCode(){
        return Objects.hash(p1,p2);
    }
    
    public Edge swappedEdge(){
        return new Edge(p2,p1);
    }
    
    @Override
    public String toString(){
        return String.format("E: p1: %h, p2: %h", p1,p2);
    }
    
    
    
    
}
