/*
 * Copyright (c) 2017 Fachhochschule Nordwestschweiz (FHNW)
 * All Rights Reserved. 
 */

package jdraw.figures;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.geom.Line2D;
import java.util.*;

import jdraw.framework.FigureHandle;
import jdraw.states.LineHandle1;
import jdraw.states.LineHandle2;



/**
 * Represents rectangles in JDraw.
 * 
 * @author Christoph Denzler
 *
 */
public class Line extends AbstractFigure {
	/**
	 * 
	 */
	private int x1;
	private int x2;
	private int y1;
	private int y2;
	private int dx,dy;
	
	private static final long serialVersionUID = 1L;
	/**
	 * Use the java.awt.Rectangle in order to save/reuse code.
	 */
	private Line2D line;
	ArrayList<FigureHandle> handles = new ArrayList<FigureHandle>(2);

	/**
	 * Create a new rectangle of the given dimension.
	 * @param x the x-coordinate of the upper left corner of the rectangle
	 * @param y the y-coordinate of the upper left corner of the rectangle
	 * @param w the rectangle's width
	 * @param h the rectangle's height
	 */
	public Line(int x1, int y1, int x2, int y2) {
		this.x1=x1;
		this.x2=x2;
		this.y1=y1;
		this.y2=y2;
		line = new Line2D.Float(x1, y1, x2, y2);
		setHandles();
	}

	public Line(Point p){
		this((int) p.getX(),(int) p.getY(), (int) p.getX(),(int) p.getY());
		setHandles();
	}
	
	/**
	 * Draw the rectangle to the given graphics context.
	 * @param g the graphics context to use for drawing.
	 */
	@Override
	public void draw(Graphics g) {
		g.setColor(Color.BLACK);
		g.drawLine((int) line.getX1(),(int) line.getY1(),(int) line.getX2(),(int) line.getY2());
	}
	
	void setHandles(){
		handles.add(new Handle((int) line.getX1(),(int) line.getY1(),this, new LineHandle1(this))); 
		handles.add(new Handle((int) line.getX2(),(int) line.getY2(),this, new LineHandle2(this))); 		
	}
	
	@Override
	public void setBounds(Point origin, Point corner) {
		line.setLine(origin, corner);
		this.x1 = (int) origin.getX();
		this.x2 = (int) corner.getX();
		this.y1 = (int) origin.getY();
		this.y2 = (int) corner.getY();
		propagateFigureEvent();
	}

	@Override
	public void move(int dx, int dy) {
		if(dx!=0 || dy!=0){
			this.dx = dx;
			this.dy = dy;
			line.setLine(line.getX1()+dx, line.getY1()+dy, line.getX2()+dx, line.getY2()+dy);
			this.x1 = this.x1 + dx;
			this.x2 = this.x2 + dx;
			this.y1 = this.y1 + dy;
			this.y2 = this.y2 + dy;
			propagateFigureEvent();
		}
	}

	private static final int TOL = 6;
	
	public boolean contains(int x, int y) { 
		return line.ptSegDist(x, y) <= TOL;
	}

	//@Override
	public Rectangle getBounds() {
		return line.getBounds();
	}

	/**
	 * Returns a list of 8 handles for this Rectangle.
	 * @return all handles that are attached to the targeted figure.
	 * @see jdraw.framework.Figure#getHandles()
	 */	
	@Override
	public List<FigureHandle> getHandles() {
		if(handles == null){
			handles = new ArrayList<FigureHandle>(2);
			handles.add(new Handle((int) line.getX1(),(int) line.getY1(),this, new LineHandle1(this))); 
			handles.add(new Handle((int) line.getX2(),(int) line.getY2(),this, new LineHandle2(this))); 
		}
		return handles;
	}


	@Override
	public Line clone() {
		Line copy= (Line) super.clone();
		copy.line = (Line2D) copy.line.clone();
		copy.x1 = x1;
		copy.x2 = x2;
		copy.y1 = y1;
		copy.y2 = y2;
		copy.handles = null;
		return copy;
	}

	@Override
	public void swapHorizontal() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void swapVertical() {
		// TODO Auto-generated method stub
		
	}
	
	public int getX1(){
		return this.x1;
	}
	
	public int getY1(){
		return this.y1;
	}
	
	public int getX2(){
		return this.x2;
	}
	
	public int getY2(){
		return this.y2;
	}

	@Override
	public int getDx() {
		return dx;
	}

	@Override
	public int getDy() {
		return dy;
	}

}
