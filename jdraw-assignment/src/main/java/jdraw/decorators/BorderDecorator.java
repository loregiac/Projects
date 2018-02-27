package jdraw.decorators;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Rectangle;

import jdraw.framework.Figure;

public class BorderDecorator extends DecoratorFigure {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private static final int OFFSET = 3;
	//private static int count = 0;
	
	public BorderDecorator(Figure figure){
		super(figure);
		//BorderDecorator.count = BorderDecorator.count+1;
	}
	
	@Override
	public Rectangle getBounds(){
		Rectangle r = super.getInner().getBounds();
		r.grow(OFFSET, OFFSET);
		return r;
	}
	
	@Override 
	public boolean contains(int x,int y){
		return getBounds().contains(x,y);
	}
	@Override
	public void draw(Graphics g){
		super.getInner().draw(g);
		Rectangle r = this.getBounds();
		int x1 = r.x;
		int y1 = r.y;
		int x2 = (int) (r.x+r.getWidth());
		int y2 = r.y;
		int x3 = (int) (r.x+r.getWidth());
		int y3 = (int) (r.y+r.getHeight());
		int x4 = r.x;
		int y4 = (int) (r.y +r.getHeight());
		g.setColor(Color.GRAY);
		g.drawLine(x1, y1, x2, y2);
		g.drawLine(x1, y1, x4, y4);
		g.setColor(Color.BLACK);
		g.drawLine(x2, y2,x3,y3);
		g.drawLine(x3, y3, x4, y4);
	}
}
