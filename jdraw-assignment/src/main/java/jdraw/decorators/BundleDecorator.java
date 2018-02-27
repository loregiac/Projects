package jdraw.decorators;

import java.awt.Point;

import jdraw.framework.Figure;

public class BundleDecorator extends DecoratorFigure {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//private static final int OFFSET = 3;
	
	public BundleDecorator(Figure figure){
		super(figure);
	}
	
	@Override
	public void move(int dx,int dy){
	}
	
	@Override
	public void setBounds(Point origin, Point corner){
	}
}