/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package setoperations;

import java.awt.geom.Point2D;
import java.util.LinkedList;
import java.util.List;

/**
 *
 * @author mazurd
 */
public class Polygon {
    protected List<Edge> edges;
    
    public Polygon(){
        edges = new LinkedList<>();
    }
    
    public Polygon(Point2D[] pts){
        edges = new LinkedList<>();
        for (int i=0;i<pts.length-1;i++){
            edges.add(new Edge(pts[i],pts[i+1]));
        }
        edges.add(new Edge(pts[pts.length-1],pts[0]));
    }
}
