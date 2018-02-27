package jdraw.states;

import java.awt.Cursor;
import java.awt.Point;
import java.awt.event.MouseEvent;

import java.awt.Rectangle;

import jdraw.figures.HandleState;
import jdraw.framework.DrawView;
import jdraw.framework.Figure;

public class E extends HandleState {
	
	public E(Figure owner) {
		super(owner);
	}

	@Override
	public Cursor getCursor() {
		return new Cursor(Cursor.E_RESIZE_CURSOR);
	}

	@Override
	public Point getAnchor() {
		Rectangle r = getOwner().getBounds();
		return new Point(r.x+r.width,r.y+r.height/2);    //could be wrong
	}

	@Override
	public void dragInteraction(int x, int y, MouseEvent e, DrawView v) {
		Rectangle r = getOwner().getBounds();
		getOwner().setBounds(new Point(r.x,r.y), new Point(x,r.y+r.height));
		
		if(x <= r.x-r.getWidth()){
			owner.swapHorizontal();
		}
		
	}

}
