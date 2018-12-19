/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package setoperations;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.awt.geom.Path2D;
import java.util.LinkedList;
import java.util.List;
import javax.swing.BorderFactory;

/**
 *
 * @author jethro
 */
public class drawPanel extends javax.swing.JPanel {
    
    protected Polygon polyA;
    protected Polygon polyB;
    protected List<IntersectionPoint> intersec;
    protected List<Polygon> polyOut;

   
    
    public drawPanel() {
        polyA = null;
        polyB = null;
        polyOut = new LinkedList<>();
        intersec = new LinkedList<>();
        initComponents();
        setBorder(BorderFactory.createLineBorder(Color.black));
        setBackground(Color.white);
    }
    
    
    @Override
    public void paintComponent(Graphics g){
        super.paintComponent(g);
        Graphics2D gfx = (Graphics2D)g;
        
        int width = this.getWidth();
        int height = this.getHeight();
        
        if (polyA == null){
            return;
        }
        
        for (Polygon p: polyOut){
            Path2D ppath = drawSinglePoly(p);
            ppath = affineTransform(ppath, width, height);
            if (p.side == Algorithms.PositionEnum.INSIDE){
                gfx.setColor(Color.yellow);
            } else {
                gfx.setColor(Color.cyan);
            }
            gfx.fill(ppath);
        }
        
        
        Path2D[] pA = drawPoly(polyA);
        Path2D[] pB = drawPoly(polyB);
        
        Path2D pAin;
        Path2D pAout;
        Path2D pBin;
        Path2D pBout;
        
        pAin = affineTransform(pA[0], width, height);
        pBin = affineTransform(pB[0], width, height);
        pAout = affineTransform(pA[1], width, height);
        pBout = affineTransform(pB[1], width, height);
        
        gfx.setColor(Color.red);
        gfx.draw(pAin);
        gfx.setColor(Color.orange);
        gfx.draw(pAout); 
        
        gfx.setColor(Color.green);
        gfx.draw(pBin);        
        gfx.setColor(Color.blue);
        gfx.draw(pBout);
        
        
        
        
        gfx.setColor(Color.black);
        for (int i=0;i<intersec.size();i++){
            int x;
            int y;
            x = (int)(intersec.get(i).point.getX()*width);
            y = (int)((1-intersec.get(i).point.getY())*height);
            System.out.println(intersec.get(i).point.getX());
            System.out.println(x);

            gfx.drawLine(x-5, y-5, x+5, y+5);
            gfx.drawLine(x-5, y+5, x+5, y-5);

        }
    }
    
    
    private static Path2D[] drawPoly(Polygon p){
        Path2D[] paths = new Path2D[2];
        paths[0] = new Path2D.Double();
        paths[1] = new Path2D.Double();
        for (Edge e : p.edges){
            if (e.side == null || e.side == Algorithms.PositionEnum.OUTSIDE){ //outside or unknown
                paths[1].moveTo(e.start.getX(), e.start.getY());
                paths[1].lineTo(e.end.getX(), e.end.getY());
            } else { // inside
                paths[0].moveTo(e.start.getX(), e.start.getY());
                paths[0].lineTo(e.end.getX(), e.end.getY());
                System.out.println("Inside!");
            }
        }        
        return paths;
    }
    
    private static Path2D drawSinglePoly(Polygon p){
        Path2D path = new Path2D.Double();
        path.moveTo(p.edges.get(0).start.getX(), p.edges.get(0).start.getY());
        for (Edge e : p.edges){
                path.lineTo(e.end.getX(), e.end.getY());
        }        
        return path;
    }
    
    // method, which transforms and scales Path2D to canvas size
    private static Path2D affineTransform(Path2D p, int width, int height) {
        AffineTransform at = AffineTransform.getScaleInstance(width, -height);
        p.transform(at);
        at = AffineTransform.getTranslateInstance(0, height);
        p.transform(at);
        return p;
    }
    
    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 400, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 300, Short.MAX_VALUE)
        );
    }// </editor-fold>//GEN-END:initComponents


    // Variables declaration - do not modify//GEN-BEGIN:variables
    // End of variables declaration//GEN-END:variables
}
