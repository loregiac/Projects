package jdraw.states;

import java.awt.Cursor;
import java.awt.Point;
import java.awt.event.MouseEvent;

import jdraw.figures.HandleState;
import jdraw.figures.Line;
import jdraw.framework.DrawView;
import jdraw.framework.Figure;

public class LineHandle2 extends HandleState{
	public LineHandle2(Figure owner) {
		super(owner);
	}

	@Override
	public Cursor getCursor() {
		return new Cursor(Cursor.HAND_CURSOR);
	}

	@Override
	public Point getAnchor() {
		Line l = (Line) getOwner();
		return new Point((int) l.getX2(), (int) l.getY2());    //could be wrong
	}

	@Override
	public void dragInteraction(int x, int y, MouseEvent e, DrawView v) {
		Line l = (Line) getOwner();
		getOwner().setBounds(new Point((int) l.getX1(), (int) l.getY1()),new Point(x,y));
	}
}
