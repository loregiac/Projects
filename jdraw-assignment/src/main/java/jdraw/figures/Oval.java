/*
 * Copyright (c) 2017 Fachhochschule Nordwestschweiz (FHNW)
 * All Rights Reserved. 
 */

package jdraw.figures;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.Rectangle;

//import java.awt.geom.Ellipse2D;


/**
 * Represents rectangles in JDraw.
 * 
 * @author Christoph Denzler
 *
 */
public class Oval extends AbstractRectangularFigure {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	/**
	 * Use the java.awt.Rectangle in order to save/reuse code.
	 */
	//private Ellipse2D ellipse;
	/**
	 * Create a new rectangle of the given dimension.
	 * @param x the x-coordinate of the upper left corner of the rectangle
	 * @param y the y-coordinate of the upper left corner of the rectangle
	 * @param w the rectangle's width
	 * @param h the rectangle's height
	 */
	//public Oval(int x, int y, int w, int h) {
		//ellipse = new Ellipse2D.Float(x, y, w, h);
	//}

	public Oval(Point p) {
		super(p);	
	}

	/**
	 * Draw the rectangle to the given graphics context.
	 * @param g the graphics context to use for drawing.
	 */
	@Override
	public void draw(Graphics g) {
		Rectangle r = getBounds();
		g.setColor(Color.WHITE);
		g.fillOval(r.x, r.y, r.width, r.height);
		g.setColor(Color.BLACK);
		g.drawOval(r.x, r.y, r.width, r.height);
	}

}
