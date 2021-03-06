/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package convexhull;

import java.awt.Point;
import java.awt.geom.Point2D;
import static java.lang.Math.cos;
import static java.lang.Math.sin;
import java.util.Random;

/**
 *
 * @author jethro
 */
public class GUI extends javax.swing.JFrame {
    

    /**
     * Creates new form GUI
     */
    public GUI() {
        initComponents();
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        drawPanel1 = new convexhull.drawPanel();
        pointsButton = new javax.swing.JButton();
        qhButton = new javax.swing.JButton();
        pointCountField = new javax.swing.JTextField();
        jarvisButton1 = new javax.swing.JButton();
        benchmarkButton = new javax.swing.JButton();
        sweepButton1 = new javax.swing.JButton();

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);

        javax.swing.GroupLayout drawPanel1Layout = new javax.swing.GroupLayout(drawPanel1);
        drawPanel1.setLayout(drawPanel1Layout);
        drawPanel1Layout.setHorizontalGroup(
            drawPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 682, Short.MAX_VALUE)
        );
        drawPanel1Layout.setVerticalGroup(
            drawPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );

        pointsButton.setText("Points");
        pointsButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                pointsButtonActionPerformed(evt);
            }
        });

        qhButton.setText("QuickHull");
        qhButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                qhButtonActionPerformed(evt);
            }
        });

        pointCountField.setText("100");

        jarvisButton1.setText("Jarvis");
        jarvisButton1.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jarvisButton1ActionPerformed(evt);
            }
        });

        benchmarkButton.setText("Benchmark!");
        benchmarkButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                benchmarkButtonActionPerformed(evt);
            }
        });

        sweepButton1.setText("SweepLine");
        sweepButton1.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                sweepButton1ActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addContainerGap()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(pointsButton, javax.swing.GroupLayout.DEFAULT_SIZE, 157, Short.MAX_VALUE)
                            .addComponent(pointCountField)))
                    .addGroup(layout.createSequentialGroup()
                        .addGap(23, 23, 23)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(benchmarkButton, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(qhButton, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(drawPanel1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap())
            .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGroup(layout.createSequentialGroup()
                    .addGap(22, 22, 22)
                    .addComponent(jarvisButton1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGap(690, 690, 690)))
            .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGroup(layout.createSequentialGroup()
                    .addGap(33, 33, 33)
                    .addComponent(sweepButton1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGap(690, 690, 690)))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(pointCountField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(32, 32, 32)
                        .addComponent(pointsButton, javax.swing.GroupLayout.PREFERRED_SIZE, 33, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(92, 92, 92)
                        .addComponent(qhButton)
                        .addGap(78, 78, 78)
                        .addComponent(benchmarkButton)
                        .addGap(0, 334, Short.MAX_VALUE))
                    .addComponent(drawPanel1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addContainerGap())
            .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGroup(layout.createSequentialGroup()
                    .addGap(145, 145, 145)
                    .addComponent(jarvisButton1)
                    .addContainerGap(492, Short.MAX_VALUE)))
            .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGroup(layout.createSequentialGroup()
                    .addGap(241, 241, 241)
                    .addComponent(sweepButton1)
                    .addContainerGap(396, Short.MAX_VALUE)))
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void pointsButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_pointsButtonActionPerformed
        int npoints;
        npoints = Integer.parseInt(pointCountField.getText());
        
        drawPanel1.points = generateRandom(npoints);
        drawPanel1.repaint();
    }//GEN-LAST:event_pointsButtonActionPerformed

    private void qhButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_qhButtonActionPerformed
        drawPanel1.hull = Algorithms.quickHull(drawPanel1.points);
        drawPanel1.repaint();
    }//GEN-LAST:event_qhButtonActionPerformed

    private void jarvisButton1ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jarvisButton1ActionPerformed
        drawPanel1.hull = Algorithms.jarvisScan(drawPanel1.points);
        drawPanel1.repaint();
    }//GEN-LAST:event_jarvisButton1ActionPerformed

    private Point2D [] generateRandom(int size){
        Point2D [] points;
        points = new Point2D[size];
        Random rnd;
        rnd = new Random();
        for (int i=0;i<size;i++){
            
            points[i] = new Point2D.Double(rnd.nextDouble(),rnd.nextDouble());
        }
        return points;
    }
    
    private Point2D[] generateCircle(int size){
        Point2D [] points;
        points = new Point2D[size];
        Random rnd;
        rnd = new Random();
        for (int i=0;i<size;i++){
            double rand = rnd.nextDouble()*2*Math.PI;
            double x = cos(rand)/2 + 0.5;
            double y = sin(rand)/2 + 0.5;
            points[i] = new Point2D.Double(x,y);
        }
        return points;
    }
    
    private Point2D [] generateGrid(int size){
        Point2D [] grid;
            grid = new Point2D.Double[size*size];
            for (int x = 0; x < size; x++){
                for (int y = 0; y < size; y++){            
                    grid[x*size+y] = new Point2D.Double(1.0/size*x,1.0/size*y);
                }
        
        }
        return grid; 
    }
    
    private void benchmarkData(Point2D[] data){
        long startTime = System.nanoTime();
        Algorithms.quickHull(data);
        long endTime = System.nanoTime();
        long quickTime = endTime - startTime;
            
        startTime = System.nanoTime();
        //Algorithms.jarvisScan(data);
        endTime = System.nanoTime();
        long jarvisTime = endTime - startTime;
            
        startTime = System.nanoTime();
        Algorithms.sweepHull(data);
        endTime = System.nanoTime();
        long sweepTime = endTime - startTime;
            
        System.out.format("%d,%d,%d,%d\n", data.length,jarvisTime/1000,quickTime/1000,sweepTime/1000);
    }
    
    
    private void benchmarkButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_benchmarkButtonActionPerformed
        final int gridfrom = 10;
        final int gridto = 500;
        
        final int sizefrom = 100;
        final int sizeto = 250_000;
        
        System.out.println("Benchmarking grid");
        for (int size = gridfrom; size < gridto; size *= 1.1){          
            Point2D[] grid = generateGrid(size);
            benchmarkData(grid);
        }
        
        System.out.println("Benchmarking random");
        for (int size = sizefrom; size < sizeto; size *=1.1){
            Point2D [] points = generateRandom(size);
            benchmarkData(points);
        }
        
        
        System.out.println("Benchmarking circle");
        for (int size = sizefrom; size < sizeto; size *=1.1){
            Point2D [] points = generateCircle(size);
            benchmarkData(points);
        }
        
        // Vygenerovat data
        // spustit algoritmy
        // vypsat výsledek
    }//GEN-LAST:event_benchmarkButtonActionPerformed

    private void sweepButton1ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_sweepButton1ActionPerformed
        drawPanel1.hull = Algorithms.sweepHull(drawPanel1.points);
        drawPanel1.repaint();
    }//GEN-LAST:event_sweepButton1ActionPerformed

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException ex) {
            java.util.logging.Logger.getLogger(GUI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(GUI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(GUI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(GUI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>

        /* Create and display the form */
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
                new GUI().setVisible(true);
            }
        });
    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton benchmarkButton;
    private convexhull.drawPanel drawPanel1;
    private javax.swing.JButton jarvisButton1;
    private javax.swing.JTextField pointCountField;
    private javax.swing.JButton pointsButton;
    private javax.swing.JButton qhButton;
    private javax.swing.JButton sweepButton1;
    // End of variables declaration//GEN-END:variables


}
