/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package setoperations;

import java.awt.geom.Point2D;

/**
 *
 * @author mazurd
 */
public class Edge {
    protected Point2D start;
    protected Point2D end;
    protected Boolean inside;
    
    public Edge(Point2D s, Point2D e){
        start = s;
        end = e;
        inside = null;
    }
}
