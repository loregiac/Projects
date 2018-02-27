package jdraw.states;

import java.awt.Cursor;
import java.awt.Point;
import java.awt.event.MouseEvent;

import java.awt.Rectangle;

import jdraw.figures.HandleState;
import jdraw.framework.DrawView;
import jdraw.framework.Figure;

public class NEgroup extends HandleState {
	
	public NEgroup(Figure owner) {
		super(owner);
	}

	@Override
	public Cursor getCursor() {
		return new Cursor(Cursor.NE_RESIZE_CURSOR);
	}

	@Override
	public Point getAnchor() {
		Rectangle r = getOwner().getBounds();
		return new Point(r.x+r.width,r.y);    //could be wrong
	}

	@Override
	public void dragInteraction(int x, int y, MouseEvent e, DrawView v) {
		Rectangle r = getOwner().getBounds();
		getOwner().move(x-r.x-r.width , y-r.y);
	}

}

