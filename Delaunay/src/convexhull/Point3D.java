/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package convexhull;

import java.awt.geom.Point2D;

/**
 *
 * @author jethro
 */
public class Point3D {
    private double x;
    private double y;
    private double z;
    
    public Point3D(double ax, double ay, double az){
        x = ax;
        y = ay;
        z = az;
    }
    
    public Point2D toPoint2D(){
        return new Point2D.Double(x,y);
    }
    
    
    public double getX(){
        return x;
    }
    
    public double getY(){
        return y;
    }
    
    public double getZ(){
        return z;
    }

    /**
     * @param x the x to set
     */
    public void setX(double x) {
        this.x = x;
    }

    /**
     * @param y the y to set
     */
    public void setY(double y) {
        this.y = y;
    }

    /**
     * @param z the z to set
     */
    public void setZ(double z) {
        this.z = z;
    }
    
}
