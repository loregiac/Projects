package jdraw.states;

import java.awt.Cursor;
import java.awt.Point;
import java.awt.event.MouseEvent;

import java.awt.Rectangle;

import jdraw.figures.HandleState;
import jdraw.framework.DrawView;
import jdraw.framework.Figure;

public class SWgroup extends HandleState {
	
	public SWgroup(Figure owner) {
		super(owner);
	}

	@Override
	public Cursor getCursor() {
		return new Cursor(Cursor.SW_RESIZE_CURSOR);
	}

	@Override
	public Point getAnchor() {
		Rectangle r = getOwner().getBounds();
		return new Point(r.x,r.y+r.height);    //could be wrong
	}

	@Override
	public void dragInteraction(int x, int y, MouseEvent e, DrawView v) {
		Rectangle r = getOwner().getBounds();
		getOwner().move(x-r.x,y-(r.y+r.height));
	}

}
