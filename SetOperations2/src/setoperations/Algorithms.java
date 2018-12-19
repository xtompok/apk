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
import java.util.ListIterator;

/**
 *
 * @author mazurd
 */
public class Algorithms {

    public static IntersectionPoint calcIntersection(Edge e1, Edge e2) {
        double ux, uy, vx, vy, wx, wy;
        ux = e1.end.getX() - e1.start.getX();
        uy = e1.end.getY() - e1.start.getY();
        vx = e2.end.getX() - e2.start.getX();
        vy = e2.end.getY() - e2.start.getY();
        wx = e1.start.getX() - e2.start.getX();
        wy = e1.start.getY() - e2.start.getY();

        double k1, k2, k3;
        k1 = vx * wy - vy * wx;
        k2 = ux * wy - uy * wx;
        k3 = vy * ux - vx * uy;

        double alpha, beta;
        alpha = k1 / k3;
        beta = k2 / k3;

        if (alpha > 1 || beta > 1 || alpha < 0 || beta < 0) {
            return null;
        }

        Point2D pt;
        pt = new Point2D.Double(e1.start.getX() + alpha * ux, e1.start.getY() + alpha * uy);

        return new IntersectionPoint(pt, e1, e2, alpha, beta);
    }

    public static List<IntersectionPoint> allIntersections(Polygon polyA, Polygon polyB) {
        List<IntersectionPoint> ints;
        ints = new LinkedList<>();

        for (Edge ea : polyA.edges) {
            for (Edge eb : polyB.edges) {
                IntersectionPoint pt;
                pt = calcIntersection(ea, eb);
                if (pt != null) {
                    ints.add(pt);
                }
            }
        }

        return ints;
    }

    public static List<Edge> divideEdge(Edge e, List<IntersectionPoint> points) {
        List<IntersectionPoint> myInts;
        myInts = new LinkedList<>();

        for (IntersectionPoint pt : points) {
            if (pt.e1 == e || pt.e2 == e) {
                myInts.add(pt);
            }
        }

        Collections.sort(myInts, (IntersectionPoint t, IntersectionPoint t1) -> {
            double c1, c2;
            if (t.e1 == e) {
                c1 = t.alpha;
            } else {
                c1 = t.beta;
            }
            if (t1.e1 == e) {
                c2 = t1.alpha;
            } else {
                c2 = t1.beta;
            }
            if (c1 < c2) {
                return -1;
            } else if (c1 > c2) {
                return 1;
            }
            return 0;
        });

        List<Edge> divided;
        divided = new LinkedList<>();

        if (myInts.size() == 0) {
            divided.add(e);
            return divided;
        }

        divided.add(new Edge(e.start, myInts.get(0).point));
        for (int i = 0; i < myInts.size() - 1; i++) {
            divided.add(new Edge(myInts.get(i).point, myInts.get(i + 1).point));
        }
        divided.add(new Edge(myInts.get(myInts.size() - 1).point, e.end));

        return divided;
    }

    public static Polygon divideAll(Polygon poly, List<IntersectionPoint> points) {
        Polygon out;
        out = new Polygon();
        for (Edge e : poly.edges) {
            out.edges.addAll(divideEdge(e, points));
        }
        return out;
    }

    public enum PositionEnum {
        INSIDE, OUTSIDE, BOUNDARY
    }

    public static double dotProd(double ux, double uy, double vx, double vy) {
        return ux * vx + uy * vy;
    }

    public static double len(double ux, double uy) {
        return sqrt(ux * ux + uy * uy);
    }

    public static double angle(double ux, double uy, double vx, double vy) {
        double prod = dotProd(ux, uy, vx, vy);
        double ulen = len(ux, uy);
        double vlen = len(vx, vy);
        return acos(prod / (ulen * vlen));
    }

    public static int deter(double x1, double y1, double x2, double y2,
            double x, double y) {
        double det = (x2 - x1) * (y - y1) - (y2 - y1) * (x - x1);
        //if determinant is close to 0 (10) it takes it as boundary
        // this interval is really small - need to use decimal numbers, to get there
        if (det <= -0.001 && det >= 0.001) {
            if (x >= x2 && x <= x1 || x <= x2 && x >= x1) {
                if (y >= y2 && y <= y1 || y <= y2 && y >= y1) {
                    return 0;
                }
                return 3;
            }
            return 3;
        } else if (det > 0) {
            return 1;
        } else {
            return 2;
        }
    }
    protected static final double EPSILON = 0.001;

    protected enum OrientationEnum {
        CW, COLINEAR, CCW
    }

    // method, which does the determinant test
    // and determinates orientation of third point based on a test
    public static OrientationEnum getOrientation(Point2D p1, Point2D p2, Point2D p3) {
        double val = (p2.getY() - p1.getY()) * (p3.getX() - p2.getX())
                - (p2.getX() - p1.getX()) * (p3.getY() - p2.getY());
        if (abs(val) < EPSILON) {
            return OrientationEnum.COLINEAR;
        } else if (val > 0) {
            return OrientationEnum.CW;
        } else {
            return OrientationEnum.CCW;
        }
    }

    public static PositionEnum pointPolygonWinding(Point2D pt, Polygon poly) {
        double sumAngle = 0;
        final double eps = 0.01;

        for (Edge e : poly.edges) {
            double ux = e.start.getX() - pt.getX();
            double uy = e.start.getY() - pt.getY();

            double vx = e.end.getX() - pt.getX();
            double vy = e.end.getY() - pt.getY();

            if (getOrientation(pt, e.start, e.end) == OrientationEnum.CCW) {
                // count angle between two vectors and deduct it to sumAngle
                sumAngle = sumAngle - angle(ux, uy, vx, vy);
            } else {
                // count angle between two vectors and add it up to sumAngle
                sumAngle = sumAngle + angle(ux, uy, vx, vy);
            }
        }

        if (abs(sumAngle) >= 2 * Math.PI - eps) {
            return PositionEnum.INSIDE;
        }
        return PositionEnum.OUTSIDE;
    }

    public static OrientationEnum getPolygonOrientation(Polygon poly) {
        double sumAngle = 0;
        final double eps = 0.01;

        Point2D pt;
        pt = poly.edges.get(0).start;

        List<Edge> otherEdges;
        otherEdges = poly.edges.subList(1, poly.edges.size() - 1);

        for (Edge e : otherEdges) {
            double ux = e.start.getX() - pt.getX();
            double uy = e.start.getY() - pt.getY();

            double vx = e.end.getX() - pt.getX();
            double vy = e.end.getY() - pt.getY();

            if (getOrientation(pt, e.start, e.end) == OrientationEnum.CCW) {
                // count angle between two vectors and deduct it to sumAngle
                sumAngle = sumAngle - angle(ux, uy, vx, vy);
            } else {
                // count angle between two vectors and add it up to sumAngle
                sumAngle = sumAngle + angle(ux, uy, vx, vy);
            }
        }

        if (sumAngle > 0) {
            return OrientationEnum.CCW;
        }
        return OrientationEnum.CW;
    }

    public static void setInside(Edge e, Polygon poly) {
        Point2D middle;
        middle = new Point2D.Double(
                (e.start.getX() + e.end.getX()) / 2,
                (e.start.getY() + e.end.getY()) / 2);

        e.side = pointPolygonWinding(middle, poly);
    }

    public static void setInside(Polygon polyA, Polygon polyB) {
        for (Edge e : polyA.edges) {
            setInside(e, polyB);
        }
        for (Edge e : polyB.edges) {
            setInside(e, polyA);
        }
    }

    public static List<Edge> filterEdges(Polygon poly, PositionEnum side) {
        List<Edge> out = new LinkedList<>();
        for (Edge e : poly.edges) {
            if (e.side == side) {
                out.add(e);
            }
        }
        return out;
    }

    public static List<Polygon> buildRings(List<Edge> edges) {
        List<Polygon> out;
        out = new LinkedList<>();

        Polygon poly;
        poly = new Polygon();

        Edge cure;
        cure = edges.get(0);
        edges.remove(0);
        poly.edges.add(cure);

        while (!edges.isEmpty()) {
            Edge nexte;
            nexte = null;
            for (Edge e : edges) {
                if (cure.end.equals(e.start)) {
                    nexte = e;
                    break;
                }
            }
            if (nexte == null) {
                // mam hotovy ring
                System.out.format("Edges: %d %h %h\n", poly.edges.size(),
                        poly.edges.get(poly.edges.size() - 1).end,
                        poly.edges.get(0).start);
                if (poly.edges.get(poly.edges.size() - 1).end.equals(poly.edges.get(0).start)) {
                    if (getPolygonOrientation(poly) == OrientationEnum.CCW) {
                        poly.side = PositionEnum.OUTSIDE;
                    } else {
                        poly.side = PositionEnum.INSIDE;
                    }
                    out.add(poly);

                    poly = new Polygon();
                    cure = edges.get(0);
                    edges.remove(0);
                    poly.edges.add(cure);
                } else {
                    System.err.println("Nexte null");
                }    
            }else {
                edges.remove(nexte);
                poly.edges.add(nexte);
                cure = nexte;
            }
        }

        // pridame posledni ring
        if (getPolygonOrientation(poly) == OrientationEnum.CCW) {
            poly.side = PositionEnum.OUTSIDE;
        } else {
            poly.side = PositionEnum.INSIDE;
        }
        out.add(poly);

        return out;
    }

    public static List<Polygon> polyUnion(Polygon polyA, Polygon polyB) {
        List<Edge> edges;
        edges = new LinkedList<>();

        edges.addAll(filterEdges(polyA, PositionEnum.OUTSIDE));
        edges.addAll(filterEdges(polyB, PositionEnum.OUTSIDE));

        return buildRings(edges);
    }

    public static List<Polygon> polyIntersect(Polygon polyA, Polygon polyB) {
        List<Edge> edges;
        edges = new LinkedList<>();

        edges.addAll(filterEdges(polyA, PositionEnum.INSIDE));
        edges.addAll(filterEdges(polyB, PositionEnum.INSIDE));

        return buildRings(edges);
    }

    public static List<Polygon> polyDiff(Polygon polyA, Polygon polyB) {
        List<Edge> edges;
        edges = new LinkedList<>();

        edges.addAll(filterEdges(polyA, PositionEnum.INSIDE));
        edges.addAll(filterEdges(polyB, PositionEnum.OUTSIDE));

        return buildRings(edges);
    }

}
