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
public class IntersectionPoint {
    protected Point2D point;
    protected Edge e1;
    protected Edge e2;
    protected double alpha; // for e1
    protected double beta; // for e2
    
    public IntersectionPoint(Point2D pt, Edge e1, Edge e2, double al, double be){
        this.point = pt;
        this.e1 = e1;
        this.e2 = e2;
        this.alpha = al;
        this.beta = be;
    }
    
}
